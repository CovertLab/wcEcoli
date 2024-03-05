# Running the model

See the top level [README](../README.md) for general instructions and [docs/README](README.md) for details on setup and running
and tradeoffs between different ways to run the model.


## FireWorks on Google Cloud

See [How to run the Whole Cell Model on the Google Cloud Platform](google-cloud.md)
for instructions to run a FireWorks workflow of cell simulations and analysis plots in Google Cloud or your local computer.

## Fireworks on Sherlock

See [Setting up to run FireWorks](wholecell/fireworks/README.md) for instructions to run a FireWorks workflow of cell simulations and analysis plots on Sherlock or your local computer.

**NOTE:** If you get this error message from the pymongo library trying to connect to a MongoDB server:

> pymongo.errors.OperationFailure: This MongoDB deployment does not support retryable writes. Please add retryWrites=false to your connection string.

the fix is to add the query parameter `retryWrites=true` to the `host: mongodb+srv://...` URL in your `my_launchpad.yaml` file or add these lines to `my_launchpad.yaml`:

```yaml
mongoclient_kwargs:
  retryWrites: false
```

## Using the Manual Runscripts

Here's a summary. See the top level [README](../README.md) for more info.

You can accomplish largely the same thing that running via Fireworks does by running:

1. `runParca.py` -- creates the `kb/*` sim data files for the control case
2. `runSim.py` -- creates the modified modified sim data files for requested variants and runs the sims
3. `analysisParca.py`, `analysisCohort.py`, `analysisMultigen.py`, `analysisSingle.py`, `analysisVariant.py`, `analysisComparison.py` -- generate output plots

When using the manual runscripts:

* You can run it in a debugger.
* It runs directly without FireWorks or MongoDB, but it doesn't distribute multiple sims and analysis plots across multiple computers.
* You can run just the parts you want and rerun them as needed.
* It doesn't do automate dependency management such as simulation outputs as analysis plot inputs. It's on you to rerun code if things change. (Some analysis scripts get confused if the sim runs are more varied than expected. See [https://github.com/CovertLab/wcEcoli/issues/199](#199).)
* The scripts have command line help and are more rigorous at arg parsing than `fw_queue.py`. Argparse provides helpful features like the ability to use any unambiguous option name prefix so, e.g. you can use `--var` as short for `--variant_index`.
* The manual runscripts do smart defaulting. runParca defaults to creating the sim dir `out/manual/` rather than a timestamped directory name. The others default to finding the most likely sim dir, i.e. the latest timestamped subdirectory of `out/` or else the alphabetically first subdirectory of `out/`. You can pass in a directory name that's absolute, or relative to `wcEcoli/`, or relative to `wcEcoli/out/`.
* Analysis scripts can accept a parameter like `--variant_index 1` to look for a subdirectory like `condition_000001/`.
* You can run a single analysis script via

       python models/ecoli/analysis/cohort/growthDynamics.py

   or run one or more of them like this:

       runscripts/manual/analysisCohort.py growthDynamics transcriptFrequency
* Manual analysis scripts accept range arguments like  
  ```
  python runscripts/manual/analysisVariant.py --variant-path-range START_VARIANT END_VARIANT --seed-path-range START_SEED END_SEED --generation-path-range START_GENERATION END_GENERATION
  ```  
  where the all-caps words are placeholders for your chosen values.
* The manual scripts don't compress the output files into tar files.
