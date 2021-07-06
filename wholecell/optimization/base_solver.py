"""

"""

import os
import pickle
from typing import Any, Dict, Optional

from wholecell.fireworks.firetasks import FitSimDataTask, InitRawDataTask, SimulationTask, SimulationDaughterTask, VariantSimDataTask
from wholecell.sim.simulation import ALTERNATE_KWARG_NAMES
from wholecell.utils import constants, data, scriptBase
import wholecell.utils.filepath as fp


def run_sim(args):
    for alt, original in ALTERNATE_KWARG_NAMES.items():
        if alt not in args:
            args[alt] = args[original]

    sim_args = data.select_keys(args, scriptBase.SIM_KEYS)

    variant_type = args['variant type']
    variant_directory = os.path.join(args['sim dir'], '{}_{:06n}'
        .format(variant_type, args['variant']))
    variant_sim_data_directory = os.path.join(variant_directory,
        VariantSimDataTask.OUTPUT_SUBDIR_KB)

    variant_sim_data_modified_file = os.path.join(
        variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

    for j in range(args['seed'], args['seed'] + args['init_sims']):  # init sim seeds, TODO: allow for multiple seeds with index
        seed_directory = fp.makedirs(variant_directory, "%06d" % args['index'])

        for k in range(args['generations']):  # generation number k
            gen_directory = fp.makedirs(seed_directory,
                "generation_%06d" % k)

            # l is the daughter number among all of this generation's cells,
            # which is 0 for single-daughters but would span range(2**k) if
            # each parent had 2 daughters.
            l = 0
            cell_directory = fp.makedirs(gen_directory, "%06d" % l)
            cell_sim_out_directory = fp.makedirs(cell_directory, "simOut")

            options = dict(sim_args,
                input_sim_data=variant_sim_data_modified_file,
                output_directory=cell_sim_out_directory,
                )

            if k == 0:
                task = SimulationTask(seed=j, **options)
            else:
                parent_gen_directory = os.path.join(seed_directory,
                    "generation_%06d" % (k - 1))
                parent_cell_directory = os.path.join(parent_gen_directory,
                    "%06d" % (l // 2))
                parent_cell_sim_out_directory = os.path.join(
                    parent_cell_directory, "simOut")
                daughter_state_path = os.path.join(
                    parent_cell_sim_out_directory,
                    constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))
                task = SimulationDaughterTask(
                    seed=(j + 1) * ((2 ** k - 1) + l),
                    inherited_state_path=daughter_state_path,
                    **options
                    )
            task.run_task({})

    return variant_directory

def run_parca(args, raw_data_file, sim_data_file, metrics_data_file):
    task = FitSimDataTask(
        input_data=raw_data_file,
        output_data=sim_data_file,
        output_metrics_data=metrics_data_file,
        cached=False,
        load_intermediate=None,
        save_intermediates=False,
        intermediates_directory=os.path.dirname(sim_data_file),
        cpus=args.get('cpus', 1),
        debug=False,
        variable_elongation_transcription=args.get('variable_elongation_transcription', False),
        variable_elongation_translation=args.get('variable_elongation_translation', False),
        disable_ribosome_capacity_fitting=not args.get('ribosome_fitting', False),
        disable_rnapoly_capacity_fitting=not args.get('rnapoly_fitting', False))

    task.run_task({})


class BaseSolver():
    def __init__(self, method, sim_dir, **options):
        self._method = method
        self._sim_dir = sim_dir

    def parameter_update(self, param, objectives, paths, reference_path):
        raise NotImplementedError('Need to implement in a subclass.')

    def get_parameter_perturbations(self, iteration, index):
        raise NotImplementedError('Need to implement in a subclass.')

    def n_variants_per_iteration(self):
        raise NotImplementedError('Need to implement in a subclass.')

    def perturb_parameters(self, iteration, variants, raw_data_file, sim_data_file):
        sim_data_files = []
        for variant in variants:
            index = variant - variants[0]
            raw_updates, sim_updates = self.get_parameter_perturbations(iteration, index)
            new_raw_data_file, new_sim_data_file, metrics_file = self.data_paths(variant)
            if raw_updates:
                self.apply_updates(raw_data_file, raw_updates, new_raw_data_file)
                run_parca(self._method.parca_args, raw_data_file, sim_data_file, metrics_file)
                sim_data_file = new_sim_data_file

            # TODO: apply variant modifications before updating sim_data
            self.apply_updates(sim_data_file, sim_updates, new_sim_data_file)

            sim_data_files.append(new_sim_data_file)

        return sim_data_files

    def update_parameters(self, variants, objectives):
        raw_data_paths = []
        sim_data_paths = []
        for variant in variants:
            raw_data, sim_data, _ = self.data_paths(variant)
            raw_data_paths.append(raw_data)
            sim_data_paths.append(sim_data)

        self._method.update_raw_data(objectives, raw_data_paths)
        self._method.update_sim_data(objectives, sim_data_paths)

    def run_sims(self, sim_params):
        # TODO: run sims in parallel
        sim_dirs = []
        for params in sim_params:
            sim_dir = run_sim(params)
            sim_dirs.append(sim_dir)

        return sim_dirs

    def run(self, iteration, variant):
        variants = list(range(variant, variant+self.n_variants_per_iteration()))
        raw_data_file, sim_data_file = self.get_reference_parameters(variant)
        if not self._method.initialized:
            self._method.initialize(raw_data_file, sim_data_file)
        sim_data_files = self.perturb_parameters(iteration, variants, raw_data_file, sim_data_file)
        sim_params = self._method.get_sim_params(self._sim_dir, variants)
        sim_out_dirs = self.run_sims(sim_params)
        objectives = self._method.get_objective(sim_out_dirs, sim_data_files)
        self.update_parameters(variants, objectives)
        new_variant = variants[-1] + 1

        return new_variant, objectives

    def print_update(self, objective):
        print(objective)
        self._method.print_update()

    def get_param(self, param, path):
        with open(path, 'rb') as f:
            obj = pickle.load(f)

        return self._method.get_attrs(obj, param)

    def apply_updates(self, old_path: str, updates: Dict[str, Any], new_path: str):
        """
        Apply parameter updates to attributes in a pickle object and save the new object.

        Args:
            old_path: path to the old pickle object that will be modified
            updates: updates to apply to attributes in the old pickle object,
                nested attributes should be separated by '.'
            new_path: path to the new pickle object to store the modifications

        TODO:
            way to index into raw_data
        """

        with open(old_path, 'rb') as f:
            obj = pickle.load(f)

        for param, val in updates.items():
            parent = obj
            attrs = param.split('.')
            for a in attrs[:-1]:
                parent = getattr(parent, a)
            setattr(parent, attrs[-1], val)

        with open(new_path, 'wb') as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)

    def get_reference_parameters(self, variant: int) -> (str, str):
        """
        Get paths to raw_data and sim_data files for the reference variant.

        Args:
            variant: current variant number

        Returns:
            raw_data: path to the reference raw_data pickle object for this variant
            sim_data: path to the reference sim_data pickle object for this variant
        """

        if variant == 0:
            raw_data, sim_data, metrics = self.data_paths()
            if not os.path.exists(raw_data):
                InitRawDataTask(output=raw_data).run_task({})
            if not os.path.exists(sim_data):
                run_parca(self._method.parca_args, raw_data, sim_data, metrics)
        else:
            raw_data, sim_data, _ = self.data_paths(variant - 1)

        return raw_data, sim_data

    def data_paths(self, variant: Optional[int] = None) -> (str, str, str):
        """
        Get paths to raw_data and sim_data files for a given variant.

        Args:
            variant: variant number to get the path for, if variant is None,
                gets the default, unmodified files

        Returns:
            raw_data: path to the raw_data pickle object for this variant
            sim_data: path to the sim_data pickle object for this variant
            metrics: path to the output metrics pickle object for this variant
        """

        if variant is None:
            kb_dir = fp.makedirs(self._sim_dir, constants.KB_DIR)
            sim_data_filename = constants.SERIALIZED_SIM_DATA_FILENAME
        else:
            kb_dir = fp.makedirs(self.get_variant_dir(variant), VariantSimDataTask.OUTPUT_SUBDIR_KB)
            sim_data_filename = constants.SERIALIZED_SIM_DATA_MODIFIED

        raw_data = os.path.join(kb_dir, constants.SERIALIZED_RAW_DATA)
        sim_data = os.path.join(kb_dir, sim_data_filename)
        metrics = os.path.join(kb_dir, constants.SERIALIZED_METRICS_DATA_FILENAME)

        return raw_data, sim_data, metrics

    def get_variant_dir(self, variant: int) -> str:
        return os.path.join(self._sim_dir, f'{self._method.variant_name}_{variant:06n}')
