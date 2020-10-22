'''Generate expression figure for operon grant.

For usage information, run:
	python make_figures_operon_grant.py -h
'''

import argparse
import os
import sys

from vivarium_cell.analysis.analyze import Analyzer
from vivarium_cell.plots.multibody_physics import plot_tags

from colony.constants import OUT_DIR
import wholecell.utils.filepath as fp


PROTEIN_KEYS_MAP = {
    'EG11250-MONOMER[c]': 'chpS',
    'EG11222­MONOMER[c]': 'alkA',
    'G6504­MONOMER[o]': 'gfcE',
    'G7263­MONOMER[c]': 'murQ',
}
TAG_PATH_NAME_MAP = {
	('boundary', 'bulk_molecules_report', key): name
    for key, name in PROTEIN_KEYS_MAP.items()
}
FIG_OUT_DIR = os.path.join(OUT_DIR, 'figs')
FILE_EXTENSION = 'pdf'
EXPERIMENT_ID = '20201020.174038'
METADATA_FILE = 'metadata.json'


def get_metadata():
	'''Get information on which experiments and code were used.'''
	metadata = {
		'git_hash': fp.run_cmdline('git rev-parse HEAD'),
		'git_branch': fp.run_cmdline('git symbolic-ref --short HEAD'),
		'time': fp.timestamp(),
		'python': sys.version.splitlines()[0],
        'experiment_id': EXPERIMENT_ID,
	}
	return metadata


def main():
	'''Generate figure.'''
	if not os.path.exists(FIG_OUT_DIR):
		os.makedirs(FIG_OUT_DIR)
	fp.write_json_file(os.path.join(
		FIG_OUT_DIR, METADATA_FILE), get_metadata())
	parser = argparse.ArgumentParser()
	Analyzer.add_connection_args(parser)
	args = parser.parse_args()

	data, environment_config = Analyzer.get_data(
		args, EXPERIMENT_ID)
	tags_data = Analyzer.format_data_for_tags(data, environment_config)
	plot_config = {
		'out_dir': FIG_OUT_DIR,
		'tagged_molecules': TAG_PATH_NAME_MAP.keys(),
		'filename': 'operon_grant_fig6_expression.{}'.format(FILE_EXTENSION),
		'tag_path_name_map': TAG_PATH_NAME_MAP,
		'tag_label_size': 48,
		'default_font_size': 48,
        'n_snapshots': 5,
        'background_color': 'gray',
	}
	plot_tags(tags_data, plot_config)


if __name__ == '__main__':
	main()
