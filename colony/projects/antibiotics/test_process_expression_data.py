import pandas as pd

from colony.projects.antibiotics.process_expression_data import (
    raw_data_to_end_expression_table,
)


class TestRawDataToEndExpressionTable:

	@staticmethod
	def make_agent_data(volume, counts_dict):
		agent_data = {
			'boundary': {
				'volume': volume,
			},
		}
		agent_data['counts'] = dict()
		for variable, count in counts_dict.items():
			agent_data['counts'][variable] = count
		return agent_data

	def test_simple(self):
		data = {
			1: {
				'agents': {
					'agent1': self.make_agent_data(2, {'protein': 1}),
					'agent2': self.make_agent_data(4, {'protein': 4}),
				},
			},
		}
		name_to_path_map = {
			'protein': ('counts', 'protein'),
		}
		table = raw_data_to_end_expression_table(data, name_to_path_map)
		assert set(table['protein']) == set([1 / 2, 4 / 4])

	def test_get_end_time(self):
		data = {
			2: {
				'agents': {
					'agent1': self.make_agent_data(2, {'protein': 0}),
					'agent2': self.make_agent_data(4, {'protein': 0}),
				},
			},
			3: {
				'agents': {
					'agent1': self.make_agent_data(2, {'protein': 1}),
					'agent2': self.make_agent_data(4, {'protein': 4}),
				},
			},
			1: {
				'agents': {
					'agent1': self.make_agent_data(2, {'protein': 0}),
					'agent2': self.make_agent_data(4, {'protein': 0}),
				},
			},
		}
		name_to_path_map = {
			'protein': ('counts', 'protein'),
		}
		table = raw_data_to_end_expression_table(data, name_to_path_map)
		assert set(table['protein']) == set([1 / 2, 4 / 4])

	def test_multiple_proteins(self):
		data = {
			1: {
				'agents': {
					'agent1': self.make_agent_data(
						2,
						{'protein1': 1, 'protein2': 2, 'protein3': 0}
					),
					'agent2': self.make_agent_data(
						4,
						{'protein2': 3, 'protein1': 8, 'protein3': 0}
					),
				},
			},
		}
		name_to_path_map = {
			'protein1': ('counts', 'protein1'),
			'protein2': ('counts', 'protein2'),
		}
		table = raw_data_to_end_expression_table(data, name_to_path_map)
		print(table['protein1'][0])
		if table['protein1'][0] == 1 / 2:
			assert table['protein1'].to_list() == [1 / 2, 8 / 4]
			assert table['protein2'].to_list() == [2 / 2, 3 / 4]
		else:
			assert table['protein1'].to_list() == [8 / 4, 1 / 2]
			assert table['protein2'].to_list() == [3 / 4, 2 / 2]
		assert 'protein3' not in table.columns

	def test_zeros(self):
		data = {
			1: {
				'agents': {
					'agent1': self.make_agent_data(0, {'protein': 1}),
					'agent2': self.make_agent_data(4, {'protein': 0}),
					'agent3': self.make_agent_data(0, {'protein': 0}),
				},
			},
		}
		name_to_path_map = {
			'protein': ('counts', 'protein'),
		}
		table = raw_data_to_end_expression_table(data, name_to_path_map)
		assert set(table['protein']) == set([0, 0, 0])
