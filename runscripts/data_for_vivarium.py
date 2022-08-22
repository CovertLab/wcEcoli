"""Transform simulation data for comparison with Vivarium E. coli.

This script reads from wcEcoli output tables and produces JSON files
with the wcEcoli simulation data in the format expected by
vivarium-ecoli. These JSON files can then be used for analyses that
compare vivarium-ecoli and wcEcoli.
"""

import argparse
import json
import os
import pickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units


def _make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--sim_data',
        default='out/manual/kb/simData.cPickle',
        help='Path to simData.'
    )
    parser.add_argument(
        '--sim_path',
        default='out/manual/',
        help='Path to root simulation output directory.'
    )
    parser.add_argument(
        '--seed',
        nargs='+',
        default=[0],
        type=int,
        help='Path to seed for random number generators.'
    )
    parser.add_argument(
        '--out',
        default='out/',
        help='Path to folder where data files should be written.'
    )
    parser.add_argument(
        '--output_prefix',
        default='monomer_counts',
        help='Prefix for output file names.'
    )
    return parser


def main(cli_args=None):
    parser = _make_parser()
    args = parser.parse_args(cli_args)

    # Get the names of proteins from the KB.
    with open(args.sim_data, 'rb') as f:
        sim_data = pickle.load(f)
    monomer_ids = sim_data.process.translation.monomer_data["id"].tolist(),

    # Get the protein counts, averaged over time, from the output data.
    for seed in args.seed:
        sim_out_dir = os.path.join(
            args.sim_path, 'wildtype_000000', f'{seed:06d}',
            'generation_000000', '000000', 'simOut')
        monomer_counts_table = TableReader(os.path.join(sim_out_dir, "MonomerCounts"))
        mass_table = TableReader(os.path.join(sim_out_dir, "Mass"))

        monomer_counts = monomer_counts_table.readColumn("monomerCounts")
        masses = {
            submass: mass_table.readColumn(submass)
            for submass in (
                'dryMass', 'proteinMass', 'tRnaMass', 'rRnaMass',
                'mRnaMass', 'dnaMass', 'smallMoleculeMass',
            )
        }

        # Vivarium uses dry_mass instead of dryMass.
        masses['dry_mass'] = masses['dryMass']
        del masses['dryMass']

        data = {
            i * 2: {
                'listeners': {
                    'monomer_counts': counts.tolist(),
                    'mass': {
                        submass_type: timeseries[i]
                        for submass_type, timeseries in masses.items()
                    },
                },
            }
            for i, counts in enumerate(monomer_counts)
        }

        output = {
            'data': data,
            'environment_config': {},
        }
        output_path = os.path.join(
            args.out, f'{args.output_prefix}_{seed}.json')
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(output, f)


if __name__ == '__main__':
    main()
