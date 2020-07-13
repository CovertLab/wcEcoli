from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader

simOutDir = 'out/no-operon-0629/wildtype_000000/000000/generation_000000/000000/simOut'
bulkMolecules = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
mRNACounts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
main_reader = TableReader(os.path.join(simOutDir, 'Main'))
# monomerCounts_reader = TableReader(os.path.join(simOutDir, 'MonomerCounts'))

mRNA_ids = mRNACounts_reader.readAttribute('mRNA_ids')
mRNA_counts = mRNACounts_reader.readColumn('mRNA_counts')
bulkMolecule_ids = bulkMolecules.readAttribute('objectNames')
bulkMolecule_counts = bulkMolecules.readColumn('counts')
time = main_reader.readColumn('time')

coExpTime_count = 0
time_total = time.shape[0]

for i in range(time_total):
    mRNA1_curr_count = mRNA_counts[i, mRNA_ids.index('EG12197_RNA[c]')]
    mRNA2_curr_count = mRNA_counts[i, mRNA_ids.index('EG12144_RNA[c]')]
    if mRNA1_curr_count > 0 and mRNA2_curr_count > 0:
        coExpTime_count = coExpTime_count + 1

print("fraction of time when mRNAs are co-expressed: " + str(float(coExpTime_count/time_total)))

plt.figure()
plt.plot(mRNA_counts[:,mRNA_ids.index('EG12197_RNA[c]')], label='EG12197_RNA')
plt.plot(mRNA_counts[:,mRNA_ids.index('EG12144_RNA[c]')], label='EG12144_RNA')
plt.xlabel('t')
plt.ylabel('count')
plt.legend()
plt.title('EG12197, EG12144 RNA')
plt.savefig('out/mRNA_no_operon_test')
