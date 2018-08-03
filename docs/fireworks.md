Fireworks - patching to avoid slow runtimes
=======================================

FireWorks can be slowed down when refreshing a firework that depends on a lot of other fireworks. This can be alleviated by applying the following fix after pip install. For FireWorks (1.6.0) the lines of interest are from 977-999 in $PI_HOME/pyenv/versions/wcEcoli2/lib/python2.7/site-packages/fireworks/core/firework.py:

```python
# what are the parent states?
parent_states = [self.id_fw[p].state for p in self.links.parent_links.get(fw_id, [])]

completed_parent_states = ['COMPLETED']
if fw.spec.get('_allow_fizzled_parents'):
    completed_parent_states.append('FIZZLED')

if len(parent_states) != 0 and not all(s in completed_parent_states for s in parent_states):
    m_state = 'WAITING'

else:  # not DEFUSED/ARCHIVED, and all parents are done running. Now the state depends on the launch status
    # my state depends on launch whose state has the highest 'score' in STATE_RANKS
    m_launch = self._get_representative_launch(fw)
    m_state = m_launch.state if m_launch else 'READY'
    m_action = m_launch.action if (m_launch and m_launch.state == "COMPLETED") else None

    # report any FIZZLED parents if allow_fizzed allows us to handle FIZZLED jobs
    if fw.spec.get('_allow_fizzled_parents') and "_fizzled_parents" not in fw.spec:
        parent_fws = [self.id_fw[p].to_dict() for p in self.links.parent_links.get(fw_id, [])
                      if self.id_fw[p].state == 'FIZZLED']
        if len(parent_fws) > 0:
            fw.spec['_fizzled_parents'] = parent_fws
            updated_ids.add(fw_id)
```

Should be changed to:

```python
completed_parent_states = ['COMPLETED']
if fw.spec.get('_allow_fizzled_parents'):
    completed_parent_states.append('FIZZLED')

# what are the parent states?
for p in self.links.parent_links.get(fw_id, []):
    if self.id_fw[p].state not in completed_parent_states:
        m_state = 'WAITING'
        break
else:  # not DEFUSED/ARCHIVED, and all parents are done running. Now the state depends on the launch status
    # my state depends on launch whose state has the highest 'score' in STATE_RANKS
    m_launch = self._get_representative_launch(fw)
    m_state = m_launch.state if m_launch else 'READY'
    m_action = m_launch.action if (m_launch and m_launch.state == "COMPLETED") else None

    # report any FIZZLED parents if allow_fizzed allows us to handle FIZZLED jobs
    if fw.spec.get('_allow_fizzled_parents') and "_fizzled_parents" not in fw.spec:
        parent_fws = [self.id_fw[p].to_dict() for p in self.links.parent_links.get(fw_id, [])
                      if self.id_fw[p].state == 'FIZZLED']
        if len(parent_fws) > 0:
            fw.spec['_fizzled_parents'] = parent_fws
            updated_ids.add(fw_id)
```

Compile the project's Cython code.

    make clean compile

Yes, it prints deprecation warnings.
