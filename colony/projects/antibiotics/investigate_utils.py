'''Utilities for investigation scripts.
'''

import copy

from vivarium.core.experiment import get_in, assoc_path


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')


def filter_raw_data_by_time(raw_data, time_range):
	floor, ceil = time_range
	end = max(raw_data.keys())
	filtered = {
		time: time_data
		for time, time_data in raw_data.items()
		if floor * end <= time <= ceil * end
	}
	return filtered


def split_raw_data_by_survival(raw_data):
	# Establish which agents die
	agents_die = set()
	for time_data in raw_data.values():
		agents_data = get_in(time_data, PATH_TO_AGENTS)
		for agent, agent_data in agents_data.items():
			dead = get_in(agent_data, PATH_TO_DEAD, False)
			if dead:
				agents_die.add(agent)

	# Split the data
	survive_data = dict()
	for time in raw_data:
		path = (time,) + PATH_TO_AGENTS
		assoc_path(survive_data, path, dict())
	die_data = copy.deepcopy(survive_data)

	for time, time_data in raw_data.items():
		agents_data = get_in(time_data, PATH_TO_AGENTS)
		for agent, agent_data in agents_data.items():
			dest = die_data if agent in agents_die else survive_data
			path = (time,) + PATH_TO_AGENTS + (agent,)
			assoc_path(dest, path, agent_data)

	return survive_data, die_data
