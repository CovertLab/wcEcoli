#! /usr/bin/env python

"""
Pulls the lastest versions of the raw data files that are sourced from
BioCyc/EcoCyc from their webservices server and updates the local versions of
the files under reconstruction/ecoli/flat.
"""

import argparse
import os
import requests
import time

from wholecell.utils.filepath import ROOT_PATH

BASE_API_URL = 'https://websvc.biocyc.org/wc-get?type='
BASE_API_URL_PREVIEW = 'https://brg-preview.ai.sri.com/wc-get?type='
BIOCYC_FILE_IDS = [
    'complexation_reactions',
    'dna_sites',
    'equilibrium_reactions',
    'genes',
    'metabolic_reactions',
    'metabolites',
    'proteins',
    'rnas',
    'transcription_units',
    'trna_charging_reactions',
    ]
FLAT_DIR = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli', 'flat')
OK_STATUS_CODE = 200

def update_biocyc_files(preview=False):
    file_id_to_text = {}

    if preview:
        base_url = BASE_API_URL_PREVIEW
    else:
        base_url = BASE_API_URL

    # Get files from BioCyc webservices
    for file_id in BIOCYC_FILE_IDS:
        response = requests.get(base_url + file_id)
        response.raise_for_status()

        file_id_to_text[file_id] = response.text
        time.sleep(1)

    # Write each file to FLAT_DIR
    for (file_id, text) in file_id_to_text.items():
        with open(os.path.join(FLAT_DIR, file_id + '.tsv'), 'w') as f:
            f.write(text)

def parse_args():
    """
	Parses command line arguments.

	Returns:
		argparse.Namespace object: parsed arguments and values
	"""
    parser = argparse.ArgumentParser(description='Script to import latest'
        'versions of EcoCyc-sourced flat files from EcoCyc servers')

    parser.add_argument('--preview', action='store_true',
        help='Import files from the preview server.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    update_biocyc_files(preview=args.preview)
