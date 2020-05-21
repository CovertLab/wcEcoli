from wholecell.io.tablereader import TableReader
import os
import matplotlib.pyplot as plt

simOutDir = '/Users/mialydefelice/Documents/code_repositories/wcEcoli_2/wcEcoli/out/manual/wildtype_000000/000000/generation_000000/000000/simOut'

monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
ids = monomerCounts.readAttribute("monomerIds")
counts = monomerCounts.readColumn("monomerCounts")

bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
objectNames = bulkMolecules.readAttribute("objectNames")
bulkMoleculeCounts = bulkMolecules.readColumn("counts")


mRNA_names = mRNA_counts_reader.readAttribute('mRNA_ids')
mRNACounts = mRNA_counts_reader.readColumn('mRNA_counts')

plt.figure()
plt.subplot(211)
plt.plot(bulkMoleculeCounts[:,objectNames.index('EG12197-MONOMER[c]')], color='r', label='seqA')
plt.plot(bulkMoleculeCounts[:,objectNames.index('PHOSPHOGLUCMUT-MONOMER[c]')], color='g', label='pgm')
plt.legend()

plt.title('protein')

plt.subplot(212)
plt.plot(mRNACounts[:,mRNA_names.index('EG12197_RNA[c]')], color='r', label='EG12197')
plt.plot(mRNACounts[:,mRNA_names.index('EG12144_RNA[c]')], color='g', label='EG12144')
plt.plot(mRNACounts[:,mRNA_names.index('EG12197_EG12144_RNA[c]')], color='b', label='EG12197_EG12144')
plt.legend()
plt.title('rna')
plt.tight_layout()
plt.savefig('out/pgm-seqA')

#--------------------------------------------

plt.figure()
plt.subplot(211)
plt.plot(bulkMoleculeCounts[:,objectNames.index('BETAGALACTOSID-MONOMER[c]')], color='r', label='lacZ')
plt.plot(bulkMoleculeCounts[:,objectNames.index('GALACTOACETYLTRAN-MONOMER[c]')], color='g', label='lacA')
plt.plot(bulkMoleculeCounts[:,objectNames.index('LACY-MONOMER[i]')], color='b', label='lacY')
plt.legend()


plt.title('protein')

plt.subplot(212)
plt.plot(mRNACounts[:,mRNA_names.index('EG10527_EG10526_EG10524_RNA[c]')], color='r', label='ZYA_RNA')
plt.plot(mRNACounts[:,mRNA_names.index('EG10526_EG10524_RNA[c]')], color='g', label='YA_RNA')
plt.legend()
plt.title('rna')
plt.savefig('out/lacZYA_plots')


#--------------------------------------------

plt.figure()
plt.subplot(211)
plt.plot(bulkMoleculeCounts[:,objectNames.index('RIBULOKIN-MONOMER[c]')], color='r', label='araB')
plt.plot(bulkMoleculeCounts[:,objectNames.index('ARABISOM-MONOMER[c]')], color='g', label='araA')
plt.plot(bulkMoleculeCounts[:,objectNames.index('RIBULPEPIM-MONOMER[c]')], color='b', label='araD')
plt.legend()


plt.title('protein')

plt.subplot(212)
plt.plot(mRNACounts[:,mRNA_names.index('EG10053_EG10052_EG10055_RNA[c]')], color='r', label='araBAD_RNA')
plt.legend()
plt.title('rna')
plt.savefig('out/araBAD_plots')


