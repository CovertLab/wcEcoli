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
# todo: get these constants from sim_data if possible
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
        cell_specs = sim_data.cell_specs

        # Initialize attributes
        self.ribosome = []
        self.synthetase = []
        self.trna_total = []
        self.aa = []
        self.atp = []
        self.target_v_translation = []
        self.f = []

        # For each condition, gather molecule counts
        self.conditions = ["basal", "with_aa"]
        condition_ids_translation = ["minimal", "minimal_plus_amino_acids"]
        for i, condition in enumerate(self.conditions):

            # Load cell specs
            cell_specs_condition = cell_specs[condition]
            bulk_container = cell_specs_condition["bulkAverageContainer"]

            # Compute cell volume
            cell_dry_mass_init = cell_specs_condition["avgCellDryMassInit"]
            cell_dry_mass_avg = cell_dry_mass_init * INITIAL_CELL_TO_AVG_CELL
            cell_mass_avg = cell_dry_mass_avg / (1. - CELL_WATER_MASS_FRACTION)
            cell_volume = cell_mass_avg / CELL_DENSITY

            """
            Ribosomes
            Notes:
                Resembles computation in initial_conditions.py for acquiring
                number of active ribosomes:
                - number of total (ie. active + inactive) ribosomes is minimum
                  of 30s and 50s
                - 85% of total ribosomes are active
            """
            s30_id = "CPLX0-3953[c]"
            s50_id = "CPLX0-3962[c]"
            n_ribosome_total = np.min(bulk_container.counts([s30_id, s50_id]))
            n_ribosome_active = np.round(
                FRACTION_RIBOSOMES_ACTIVE * n_ribosome_total)
            self.ribosome.append(
                (n_ribosome_active / AVOGADROS_NUM / cell_volume) \
                .asNumber(units.mol / units.L))

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
            synthetase_count = bulk_container.counts(synthetase_ids)
            synthetase_concentration = synthetase_count / AVOGADROS_NUM / cell_volume
            self.synthetase.append(self._remove_units(synthetase_concentration))

            """tRNA isoacceptors"""
            trna_ids = transcription.rnaData["id"][transcription.rnaData["isTRna"]]
            trna_count = bulk_container.counts(trna_ids)
            rep_trna_count = np.dot(transcription.aa_from_trna, trna_count)
            trna_total_concentration = rep_trna_count / AVOGADROS_NUM / cell_volume
            self.trna_total.append(self._remove_units(trna_total_concentration))

            """amino acids"""
            aa_ids = sim_data.moleculeGroups.aaIDs
            aa_count = bulk_container.counts(aa_ids)
            aa_concentration = aa_count / AVOGADROS_NUM / cell_volume
            self.aa.append(self._remove_units(aa_concentration))

            """atp"""
            atp_count = bulk_container.count("ATP[c]")
            self.atp.append(
                (atp_count / AVOGADROS_NUM / cell_volume) \
                .asNumber(units.mol / units.L))

            """
            fraction of proteome that is composed of each amino acid
            Note: sums to 1
            """
            # Get counts of protein
            counts_proteome = bulk_container.counts(translation.monomerData["id"])

            # Get counts of each amino acid per protein
            polypeptide_seq = translation.translationSequences
            proteome_to_aa = np.zeros((len(aa_ids), counts_proteome.shape[0]))
            for j, seq in enumerate(polypeptide_seq):
                indices, counts = np.unique(seq, return_counts = True)
                for index, count in zip(indices, counts):
                    if index == -1:
                        continue

                    proteome_to_aa[index, j] = count

            # Compute amino acid representation across proteome
            counts_aa = np.dot(proteome_to_aa, counts_proteome)
            f = counts_aa.astype(float) / counts_aa.sum()
            self.f.append(f)

            # Targets
            elongation_rate = translation.ribosomeElongationRateDict \
                [condition_ids_translation[i]].asNumber(units.aa / units.s)
            self.target_v_translation.append(elongation_rate * self.ribosome[i])

    def _remove_units(self, array):
        return np.array([x.asNumber() for x in array])
