"""
Test new gene insertion into the chromosome.

Covers the coordinate conventions that the insertion relies on, which are easy
to break because they are only checked implicitly elsewhere:

	* new gene coordinates in new_gene_data/<dir>/genes.tsv are 1-indexed and
	  agree with the actual sequence lengths,
	* the inserted sequence reads back through the model's own
	  genome[left - 1:right] convention,
	* the genome grows by exactly the insertion length,
	* every original feature either keeps its coordinates or shifts by exactly
	  the insertion length, and none is split by the insertion,
	* the conflict search resolves to a position that splits nothing.

Building the knowledge base is slow, so the whole module shares one
``new_genes_option='off'`` build and one ``'gfp'`` build.
"""

import csv
import os
import unittest

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

FLAT_DIR = os.path.join(os.path.dirname(os.path.dirname(
	os.path.abspath(__file__))), 'ecoli', 'flat')
NEW_GENE_DIR = os.path.join(FLAT_DIR, 'new_gene_data')
NEW_GENE_OPTION = 'gfp'


def _rows(path):
	"""Rows of a flat TSV as dicts, with the quoting the flat files use."""
	with open(path) as f:
		reader = csv.reader(
			(line for line in f if not line.startswith('#')), delimiter='\t')
		header = [c.strip('"') for c in next(reader)]
		return [dict(zip(header, [c.strip('"') for c in row]))
			for row in reader if row]


def _coordinates(knowledge_base):
	"""{id: (left, right)} for every gene with usable coordinates."""
	out = {}
	for row in knowledge_base.genes:
		left, right = row['left_end_pos'], row['right_end_pos']
		if left not in (None, '') and right not in (None, ''):
			out[row['id']] = (int(left), int(right))
	return out


class TestNewGeneFlatFiles(unittest.TestCase):
	"""Checks on the input files alone; these need no knowledge base."""

	def test_coordinates_are_one_indexed(self):
		"""Every new gene set starts at relative coordinate 1, not 0."""
		for subdir in sorted(os.listdir(NEW_GENE_DIR)):
			genes_path = os.path.join(NEW_GENE_DIR, subdir, 'genes.tsv')
			if not os.path.isfile(genes_path):
				continue
			genes = _rows(genes_path)
			if not genes:
				continue        # the template directory carries no genes
			first = min(genes, key=lambda row: int(row['left_end_pos']))
			self.assertEqual(int(first['left_end_pos']), 1,
				'%s: new gene coordinates must be 1-indexed. A 0-indexed start'
				' shifts the whole insert by one base, so the modelled RNA'
				' begins one base before the insertion and loses its last'
				' base.' % subdir)

	def test_coordinates_match_sequence_lengths(self):
		"""right - left + 1 equals the length of the gene's sequence."""
		for subdir in sorted(os.listdir(NEW_GENE_DIR)):
			genes_path = os.path.join(NEW_GENE_DIR, subdir, 'genes.tsv')
			seq_path = os.path.join(NEW_GENE_DIR, subdir, 'gene_sequences.tsv')
			if not os.path.isfile(genes_path):
				continue
			genes = _rows(genes_path)
			if not genes:
				continue
			lengths = {row['id']: len(row['gene_seq'])
				for row in _rows(seq_path)}
			for gene in genes:
				span = (int(gene['right_end_pos'])
					- int(gene['left_end_pos']) + 1)
				self.assertEqual(span, lengths[gene['id']],
					'%s/%s: coordinate span %d does not match its %d nt'
					' sequence' % (subdir, gene['id'], span,
						lengths[gene['id']]))


class TestNewGeneInsertion(unittest.TestCase):
	"""Checks on the genome the knowledge base actually builds."""

	@classmethod
	def setUpClass(cls):
		cls.plain = KnowledgeBaseEcoli(
			operons_on=True, remove_rrna_operons=False, remove_rrff=False,
			stable_rrna=False, new_genes_option='off')
		cls.inserted = KnowledgeBaseEcoli(
			operons_on=True, remove_rrna_operons=False, remove_rrff=False,
			stable_rrna=False, new_genes_option=NEW_GENE_OPTION)
		new_genes = [row for row in cls.inserted.genes
			if row['id'].startswith('NG')]
		cls.new_gene = min(new_genes, key=lambda row: row['left_end_pos'])
		cls.insert_length = (len(cls.inserted.genome_sequence)
			- len(cls.plain.genome_sequence))

	def test_genome_grows_by_the_insertion_length(self):
		sequences = getattr(self.inserted.new_gene_data, NEW_GENE_OPTION)
		total = sum(len(row['gene_seq'])
			for row in sequences.gene_sequences)
		self.assertEqual(self.insert_length, total)

	def test_inserted_sequence_reads_back(self):
		"""
		The model reads any feature as genome[left - 1:right], so the new
		gene's coordinates must recover its own sequence exactly. This is the
		check that fails when the coordinates are off by one.
		"""
		sequences = getattr(self.inserted.new_gene_data, NEW_GENE_OPTION)
		expected = {row['id']: row['gene_seq']
			for row in sequences.gene_sequences}
		for row in self.inserted.genes:
			if not row['id'].startswith('NG'):
				continue
			read_back = self.inserted.genome_sequence[
				row['left_end_pos'] - 1:row['right_end_pos']]
			self.assertEqual(str(read_back), expected[row['id']],
				'%s does not read back from its own coordinates' % row['id'])

	def test_original_features_shift_by_zero_or_the_insert_length(self):
		before = _coordinates(self.plain)
		after = _coordinates(self.inserted)
		for gene_id, (old_left, old_right) in before.items():
			if gene_id not in after:
				continue
			new_left, new_right = after[gene_id]
			shift = new_left - old_left
			self.assertIn(shift, (0, self.insert_length),
				'%s moved by %d, which is neither 0 nor the insert length'
				% (gene_id, shift))
			self.assertEqual(new_right - old_right, shift,
				'%s was stretched rather than shifted' % gene_id)

	def test_no_feature_is_split_by_the_insertion(self):
		"""
		The insert goes immediately after 1-indexed base insert_pos, so no
		surviving feature may span that point.
		"""
		insert_pos = self.new_gene['left_end_pos'] - 1
		for gene_id, (left, right) in _coordinates(self.inserted).items():
			if gene_id.startswith('NG'):
				continue
			self.assertFalse(left <= insert_pos < right,
				'%s at %d-%d is split by an insertion after base %d'
				% (gene_id, left, right, insert_pos))

	def test_replication_sites_shift_with_the_genome(self):
		"""oriC and TerC set gene dosage, so a missed shift would be silent."""
		def sites(knowledge_base):
			return {row['common_name']:
				(row['left_end_pos'], row['right_end_pos'])
				for row in knowledge_base.dna_sites
				if row['common_name'] in ('oriC', 'TerC')}
		before, after = sites(self.plain), sites(self.inserted)
		insert_pos = self.new_gene['left_end_pos'] - 1
		for name, (old_left, old_right) in before.items():
			new_left, new_right = after[name]
			expected = self.insert_length if old_left > insert_pos else 0
			self.assertEqual(new_left - old_left, expected,
				'%s shifted by %d, expected %d'
				% (name, new_left - old_left, expected))
			self.assertEqual(new_right - old_right, expected)


if __name__ == '__main__':
	unittest.main()
