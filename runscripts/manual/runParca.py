"""
Run the Parameter Calculator (Parca).
The output goes into the named subdirectory of wcEcoli/out/, defaulting to "manual".

TODO: Share lots of code with fw_queue.py and AnalysisPaths.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

import os
import sys

from wholecell.fireworks.firetasks import ParcaTask
from wholecell.utils import constants, data, parallelization, scriptBase
from wholecell.utils import filepath as fp


class RunParca(scriptBase.ScriptBase):
	"""Runs the Parca, aka simulation parameter calculator."""

	def define_parameters(self, parser):
		super(RunParca, self).define_parameters(parser)

		# NOTE: Don't name this arg sim_dir since that makes parse_args() look
		# for an existing sim_dir directory while here we aim to create one.
		parser.add_argument('sim_outdir', nargs='?', default='manual',
			help='The simulation root directory to hold the "kb/" and'
				 ' eventually simulation variants and generations.'
				 ' Default = "manual".')

		parser.add_argument('--timestamp', action='store_true',
			help='Timestamp the given `sim_outdir`, transforming e.g.'
				 ' "Fast run" to "20190514.135600__Fast_run".')
		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use. Default = 1.')

		parser.add_argument('--load-intermediate',
			help='The function in the parca to load (skips functions that would'
				 ' have run before the function). Must run with --save-intermediates first.')
		parser.add_argument('--save-intermediates', action='store_true',
			help='If set, saves sim_data and cell_specs at intermediate'
			     ' function calls in the parca.')
		parser.add_argument('--intermediates-directory',
			help='Directory to save or load intermediate sim_data and cell_specs'
			     ' results from if --load-intermediate or --save-intermediates'
			     ' are set.')

		self.define_parca_options(parser)
		self.define_elongation_options(parser)

	def parse_args(self):
		args = super(RunParca, self).parse_args()
		args.cpus = parallelization.cpus(args.cpus)

		args.time = fp.timestamp()
		args.description = args.sim_outdir.replace(' ', '_')

		if args.timestamp:
			args.sim_outdir = args.time + '__' + args.description

		args.sim_path = fp.makedirs(fp.ROOT_PATH, "out", args.sim_outdir)
		return args

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, constants.KB_DIR)

		# Write the metadata file.
		metadata = {
			'git_hash': fp.git_hash(),
			'git_branch': fp.git_branch(),
			'description': args.description,
			'operons': args.operons,
			'new_genes': args.new_genes,
			'remove_rrna_operons': args.remove_rrna_operons,
			'remove_rrff': args.remove_rrff,
			'stable_rrna': args.stable_rrna,
			'time': args.time,
			'python': sys.version.splitlines()[0],
		}
		metadata_dir = fp.makedirs(args.sim_path, constants.METADATA_DIR)
		metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
		fp.write_json_file(metadata_path, metadata)

		if args.debug_parca:
			print('DEBUG Parca')

		python_args = data.select_keys(
			vars(args),
			scriptBase.PARCA_KEYS,
			debug=args.debug_parca,
			load_intermediate=args.load_intermediate,
			save_intermediates=args.save_intermediates,
			intermediates_directory=args.intermediates_directory,
			output_directory=kb_directory)

		task = ParcaTask(**python_args)
		task.run_task({})


if __name__ == '__main__':
	script = RunParca()
	script.cli()
