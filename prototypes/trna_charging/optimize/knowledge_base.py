from __future__ import division

import os
import cPickle
import numpy as np

from wholecell.utils import units
import wholecell

# The wcEcoli project root path
ROOT_PATH = os.path.dirname(
    os.path.dirname(os.path.abspath(wholecell.__file__)))

# Path to sim_data
SIM_DATA_PATH = os.path.join(
    ROOT_PATH, "out", "manual", "kb", "simData_Fit_1.cPickle")

# Constants
CELL_DENSITY = 1100. * units.g / units.L
AVOGADROS_NUM = 6.02e23 * 1 / units.mol
INITIAL_CELL_TO_AVG_CELL =  1.35
FRACTION_RIBOSOMES_ACTIVE = 0.85
CELL_WATER_MASS_FRACTION = 0.7

class knowledge_base(object):
    def __init__(self):
        # Load sim data
        sim_data = cPickle.load(open(SIM_DATA_PATH, "rb"))
        transcription = sim_data.process.transcription
        translation = sim_data.process.translation
        cell_specs = sim_data.cell_specs_basal
        bulk_container = cell_specs["bulkAverageContainer"]

        # Compute cell volume
        cell_dry_mass_init = cell_specs["avgCellDryMassInit"]
        cell_dry_mass_avg = cell_dry_mass_init * INITIAL_CELL_TO_AVG_CELL
        cell_mass_avg = cell_dry_mass_avg / (1. - CELL_WATER_MASS_FRACTION)
        cell_volume = cell_mass_avg / CELL_DENSITY

        """
        Ribosomes
        Notes:
            Resembles computation in initial_conditions.py for acquiring number
            of active ribosomes:
            - number of total (ie. active + inactive) ribosomes is minimum of
              30s and 50s
            - 85% of total ribosomes are active
        """
        s30_id = "CPLX0-3953[c]"
        s50_id = "CPLX0-3962[c]"
        n_ribosome_total = np.min(bulk_container.counts([s30_id, s50_id]))
        n_ribosome_active = np.round(
            FRACTION_RIBOSOMES_ACTIVE * n_ribosome_total)
        self.ribosome = n_ribosome_active / AVOGADROS_NUM / cell_volume

        """
        tRNA synthetases
        Note:
            Ordering of synthetase_ids matches that of aa_ids
        """
        synthetase_ids = [
            "ALAS-CPLX[c]",
            "ARGS-MONOMER[c]",
            "ASNS-CPLX[c]",
            "ASPS-CPLX[c]",
            "CYSS-MONOMER[c]",
            "GLURS-MONOMER[c]",
            "GLNS-MONOMER[c]",
            "GLYS-CPLX[c]",
            "HISS-CPLX[c]",
            "ILES-MONOMER[c]",
            "LEUS-MONOMER[c]",
            "LYSS-CPLX[c]",
            "METG-CPLX[c]",
            "PHES-CPLX[c]",
            "PROS-CPLX[c]",
            "SERS-CPLX[c]",
            "THRS-CPLX[c]",
            "TRPS-CPLX[c]",
            "TYRS-CPLX[c]",
            "CPLX0-1141[c]",
            "VALS-MONOMER[c]",
        ]
        n_synthetase = bulk_container.counts(synthetase_ids)
        self.synthetase = n_synthetase / AVOGADROS_NUM / cell_volume

        """tRNA isoacceptors"""
        trna_ids = transcription.rnaData["id"][transcription.rnaData["isTRna"]]
        n_trna = bulk_container.counts(trna_ids)
        n_trna_rep = np.dot(transcription.aa_from_trna, n_trna)
        self.trna_total = n_trna_rep / AVOGADROS_NUM / cell_volume

        """amino acids"""
        aa_ids = sim_data.moleculeGroups.aaIDs
        self.aa_ids = aa_ids
        n_aa = bulk_container.counts(aa_ids)
        self.aa = n_aa / AVOGADROS_NUM / cell_volume

        """atp"""
        n_atp = bulk_container.count("ATP[c]")
        self.atp = n_atp / AVOGADROS_NUM / cell_volume

        """
        fraction of proteome that is composed of amino acid
        Notes:
            Sums to 1
            Assumes proteome is composed of 1 count of all translatable genes
        """
        polypeptide_seq = translation.translationSequences
        aa_indices, counts = np.unique(polypeptide_seq, return_counts = True)
        count_each_aa = counts[np.where(aa_indices != -1)]
        f = count_each_aa.astype(float) / count_each_aa.sum()
        self.f = f

        # Targets
        elongation_rate = translation.ribosomeElongationRateDict["minimal"].asNumber(units.aa / units.s) # aa / s / rib
        self.target_v_ribosome = elongation_rate * self.ribosome / units.s

kb = knowledge_base()
