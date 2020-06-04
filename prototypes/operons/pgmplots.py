from wholecell.io.tablereader import TableReader
import os
import matplotlib.pyplot as plt

simOutDir = '/Users/mialydefelice/Documents/code_repositories/wcEcoli_2/wcEcoli/out/operon_1/wildtype_000000/000000/generation_000000/000000/simOut'

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
plt.plot(bulkMoleculeCounts[:,objectNames.index('EG12197-MONOMER[c]')], color='#00A2FF', label='SeqA')
plt.plot(bulkMoleculeCounts[:,objectNames.index('PHOSPHOGLUCMUT-MONOMER[c]')], color='#CB297B', label='PGM')
plt.legend()

plt.title('protein')

plt.subplot(212)
plt.plot(mRNACounts[:,mRNA_names.index('EG12144_RNA[c]')], color='#00A2FF', label='seqA')
plt.plot(mRNACounts[:,mRNA_names.index('EG12197_EG12144_RNA[c]')], color='#CB297B', label='pgm_seqA_RNA')
plt.legend()
plt.title('rna')
plt.tight_layout()
plt.savefig('out/pgm-seqA')

#--------------------------------------------

plt.figure()
plt.subplot(211)
plt.plot(bulkMoleculeCounts[:,objectNames.index('BETAGALACTOSID-MONOMER[c]')], color='#00A2FF', label='lacZ')
plt.plot(bulkMoleculeCounts[:,objectNames.index('GALACTOACETYLTRAN-MONOMER[c]')], color='#CB297B', label='lacA')
plt.plot(bulkMoleculeCounts[:,objectNames.index('LACY-MONOMER[i]')], color='#F8BA00', label='lacY')
plt.legend()


plt.title('protein')

plt.subplot(212)
plt.plot(mRNACounts[:,mRNA_names.index('EG10527_EG10526_EG10524_RNA[c]')], color='k', label='ZYA_RNA')
plt.plot(mRNACounts[:,mRNA_names.index('EG10526_EG10524_RNA[c]')], color='#00A89D', label='YA_RNA')
plt.legend()
plt.title('rna')
plt.savefig('out/lacZYA_plots')


#--------------------------------------------

plt.figure()
plt.subplot(211)
plt.plot(bulkMoleculeCounts[:,objectNames.index('RIBULOKIN-MONOMER[c]')], color='#00A2FF', label='araB')
plt.plot(bulkMoleculeCounts[:,objectNames.index('ARABISOM-MONOMER[c]')], color='#CB297B', label='araA')
plt.plot(bulkMoleculeCounts[:,objectNames.index('RIBULPEPIM-MONOMER[c]')], color='#F8BA00', label='araD')
plt.legend()


plt.title('protein')

plt.subplot(212)
plt.plot(mRNACounts[:,mRNA_names.index('EG10053_EG10052_EG10055_RNA[c]')], color='k', label='araBAD_RNA')
plt.legend()
plt.title('rna')
plt.savefig('out/araBAD_plots')

#--------------------------------------------

plt.figure()
plt.subplot(211)
plt.plot(bulkMoleculeCounts[:,objectNames.index('HYAA-MONOMER[i]')], color='#00A2FF', label='hyaA')
plt.plot(bulkMoleculeCounts[:,objectNames.index('HYAB-MONOMER[i]')], color='#CB297B', label='hyaB')
plt.plot(bulkMoleculeCounts[:,objectNames.index('HYAC-MONOMER[i]')], color='#F8BA00', label='hyaC')
plt.plot(bulkMoleculeCounts[:,objectNames.index('EG10471-MONOMER[c]')], color='#F8BA00', label='hyaD')
plt.plot(bulkMoleculeCounts[:,objectNames.index('EG10472-MONOMER[c]')], color='#F8BA00', label='hyaE')
plt.plot(bulkMoleculeCounts[:,objectNames.index('EG10473-MONOMER[c]')], color='#F8BA00', label='hyaF')
plt.legend()


plt.title('protein')

plt.subplot(212)
plt.plot(mRNACounts[:,mRNA_names.index('EG10468_EG10469_EG10470_EG10471_EG10472_EG10473_RNA[c]')], color='k', label='hyaABCDEF_RNA')
plt.legend()
plt.title('rna')
plt.savefig('out/hyaABCDEF_plots')


