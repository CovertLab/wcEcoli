"""
Make one or more sim_data variants via VariantSimDataTask, writing e.g.
`wildtype_000000/kb/simData_Modified.cPickle`, in preparation for running
cell simulations via `runSim --require_variants`. (Without --require_variants,
runSim will make the sim_data variants, which is handy until you want to launch
multiple first-gen runSim runs in parallel without collisions writing the same
`simData_Modified.cPickle` file.)

See models/ecoli/sim/variants/*.py for the variant choices.

Prerequisite: Run the parameter calculator (runParca.py).

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import os

from wholecell.fireworks.firetasks import VariantSimDataTask
from wholecell.utils import constants, scriptBase
import wholecell.utils.filepath as fp

DEFAULT_VARIANT = ['wildtype', '0', '0']


class MakeVariants(scriptBase.ScriptBase):
	"""Make sim_data variants like wildtype_000000/kb/simData_Modified.cPickle."""

	def help(self):
		"""Return help text for the Command Line Interface."""
		return '''Run {}. Given a variant type name and an index range,
			this makes the variant subdirectories and their
			simData_Modified.cPickle files.'''.format(self.description())

	def define_parameters(self, parser):
		super(MakeVariants, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)
		self.define_make_variants_option(parser)

	def run(self, args):
		sim_data1, sim_data2 = scriptBase.sim_data_paths(args.sim_path)

		variant_type = args.variant[0]
		variant_spec = (variant_type, int(args.variant[1]), int(args.variant[2]))

		# args.sim_path is called INDIV_OUT_DIRECTORY in fw_queue.
		for base_index, index, subdir in fp.iter_variants3(
				*variant_spec, both_operons=bool(sim_data2)):
			sim_data_file = sim_data1 if base_index == index else sim_data2
			variant_sim_data_directory = os.path.join(args.sim_path, subdir,
				constants.VKB_DIR)
			variant_metadata_directory = os.path.join(args.sim_path, subdir,
				constants.METADATA_DIR)

			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

			task = VariantSimDataTask(
				variant_function=variant_type,
				variant_index=base_index,
				input_sim_data=sim_data_file,
				output_sim_data=variant_sim_data_modified_file,
				variant_metadata_directory=variant_metadata_directory,
				)
			task.run_task({})


if __name__ == '__main__':
	script = MakeVariants()
	script.cli()
