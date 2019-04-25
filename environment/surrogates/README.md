# Surrogates

Light weight cell models that can be added to an environmental simulation.

# Uses

Surrogates are used for developing some feature of either the environment or the cell model that is focused on cell-environment interactions. The whole cell model can also be compared with these much simpler models to see what added predictive value it provides.

# Execution

Add surrogates into an environmental simulation using ```experiment``` from agent framework. This is the same way the whole cell model is added to an environmental simulation: 

        python -m environment.boot experiment --number N --type T
        
Here, ```T``` specifies the surrogate type, such as ```chemotaxis```, or ```transport``` and ```N``` specifies the number of cells.

They can also be added to an already-running experiment with ```add```:

    python -m environment.boot add --type T
    
# Making New Surrogates

In order to make a new surrogate, create a new python file in this directory. It requires some interface with ```agent.inner```, including the functions: 
- ```run_incremental(run_until)```, for running the cell until the time specified by run_until. 
- ```generate_inner_update()```, for passing values such as environment_change to the environment.
- ```apply_outer_update(update)```, for passing updates from the environment to the cell.

A new surrogate also requires the addition of functions to ```environment.boot```:
- ```boot_SURROGATE(agent_id, agent_type, agent_config)```, for setting up the interface with ```agent.inner```
- ```initialize_SURROGATE(agent_id, agent_type, agent_config)```, for passing configuration data to the surrogate, calling the boot function, and starting the simulation.

Finally, it also requires an initializer:
- ```initializers['SURROGATE'] = initialize_SURROGATE```

# The Transport Surrogate

The transport surrogate (laid out in surrogates/transport.py) utilize a lookup table derived from the whole-cell model to simulate the flux of molecules through the cell membrane in different media conditions. The transport surrogate utilizes a timeline function (condition/timelines) to change the media environment at different time points. At each time point, the transport surrogate updates its local media condition, then parses a lookup table associated with that media condition (condition/tables/*.tsv).

## The Look-up Table
The lookup tables record each molecular transport reaction, along with its molecules (substrates) and their ratios (stoichiometry). These ratios are multiplied by the flux, which is sampled by the transport surrogate from a list of values from multiple whole-cell model experiments.

## Pre- and Post-Processing the data
The lookup tables require some pre-processing to read the keys correctly in the nested stoichiometry and flux distribution dictionaries. This is done outside of the source code, as it requires ```ast.literal_eval()``` which is potentially very dangerous. Additionally, a CSV file is generated in out/manual that tracks some information about each cell (see ```transport_listener()```), but as everyone's use for these data will be slightly different, the post-processing scripts will not be included until the surrogate model is being prepared for wide dissemination. 
