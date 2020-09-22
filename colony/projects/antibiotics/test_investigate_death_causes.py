from colony.projects.antibiotics.investigate_death_causes import (
	filter_raw_data_by_time,
	split_raw_data_by_survival,
)

class TestSplitRawDataBySurvival:

	dead_agent = {
		'boundary': {
			'dead': True,
		},
	}
	alive_agent = {
		'boundary': {
			'dead': False,
		},
	}

	def test_simple(self):
		data = {
			1: {
				'agents': {
					'survives': self.alive_agent,
					'dies': self.alive_agent,
				},
			},
			2: {
				'agents': {
					'survives': self.alive_agent,
					'dies': self.dead_agent,
				},
			},
		}
		survive, die = split_raw_data_by_survival(data)

		survive_expected = {
			1: {
				'agents': {
					'survives': self.alive_agent,
				},
			},
			2: {
				'agents': {
					'survives': self.alive_agent,
				},
			},
		}
		die_expected = {
			1: {
				'agents': {
					'dies': self.alive_agent,
				},
			},
			2: {
				'agents': {
					'dies': self.dead_agent,
				},
			},
		}
		assert survive_expected == survive
		assert die_expected == die

	def test_all_survive(self):
		data = {
			1: {
				'agents': {
					'a': self.alive_agent,
					'b': self.alive_agent,
				},
			},
			2: {
				'agents': {
					'a': self.alive_agent,
					'b': self.alive_agent,
				},
			},
		}
		survive, die = split_raw_data_by_survival(data)

		die_expected = {
			1: {
				'agents': {},
			},
			2: {
				'agents': {},
			},
		}
		assert data == survive
		assert die_expected == die

	def test_all_die(self):
		data = {
			1: {
				'agents': {
					'survives': self.alive_agent,
					'dies': self.alive_agent,
				},
			},
			2: {
				'agents': {
					'survives': self.dead_agent,
					'dies': self.dead_agent,
				},
			},
		}
		survive, die = split_raw_data_by_survival(data)

		survive_expected = {
			1: {
				'agents': {},
			},
			2: {
				'agents': {},
			},
		}
		assert survive_expected == survive
		assert data == die


class TestFilterRawDataByTime:

	@staticmethod
	def _gen_data(times):
		data = {
			time: {
				'agents': {},
			}
			for time in times
		}
		return data

	def test_simple(self):
		data = self._gen_data([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
		filtered = filter_raw_data_by_time(data, (0.2, 0.8))
		expected = self._gen_data([2, 3, 4, 5, 6, 7, 8])
		assert expected == filtered

	def test_filter_none(self):
		data = self._gen_data(range(7))
		filtered = filter_raw_data_by_time(data, (0, 1))
		assert data == filtered

	def test_filter_all(self):
		data = self._gen_data([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
		filtered = filter_raw_data_by_time(data, (0.5, 0.5))
		expected = self._gen_data([5])
		assert expected == filtered
