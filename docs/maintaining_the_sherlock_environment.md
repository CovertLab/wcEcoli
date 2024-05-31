# Maintaining the shared Sherlock runtime environment

_Sherlock_ is Stanford's scientific computing cluster. The Covert lab runs development workflows on Sherlock, also Continuous Integration (CI) builds, including nightly builds and Pull Request (PR) builds.

For consistent results and debugging, it's important for our local builds, development builds on Sherlock, and the CI builds to run with the same versions of Python, native libraries, and Python libraries (pips).

When we update to newer versions, it requires skimming release notes for required code changes, working out a compatible set of library versions, and coordinating to update all of our development environments.

This document describes best practices for updating the versions of Python, native libraries, environment modules, and Python libraries used in the shared Sherlock runtime environment to reduce the chances of breaking things for the CI builds and our development runs.

(BTW, since Sherlock uses an NFS file server, each directory has a hidden `.snapshot/` subdirectory containing previous snapshots of files. They can be very helpful to recover from goofs.)


Also see:
* [Required development tools](dev-tools.md) -- installation and tips
* [Creating the "pyenv" runtime environment](create-pyenv.md) -- build your local Python runtime environment
* [Setting up to run FireWorks](wholecell/fireworks/README.md) -- instructions to run a FireWorks workflow of cell simulations and analysis plots


## User profile setup

Put these lines in your `.bash_profile` on Sherlock:

```bash
export PI_HOME=$GROUP_HOME
export PI_SCRATCH=$GROUP_SCRATCH

##### Add group-wide path settings #####
if [ -f "${PI_HOME}/etc/bash_profile" ]; then
    . "${PI_HOME}/etc/bash_profile"
fi

##### For pyenv #####
export PYENV_ROOT="${PI_HOME}/pyenv"

# Environment modules used by wcEcoli
module load git/2.39.1 git-lfs/2.11.0 system parallel
module load wcEcoli/python3

# for wcEcoli
export PYTHONPATH="$HOME/wcEcoli"

if [ -d "${PYENV_ROOT}" ]; then
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
    eval "$(pyenv virtualenv-init -)"
fi
```

(The `parallel` module is needed only for rarely used shell scripts.)

The `$PYENV_ROOT` directory holds the lab's shared collection of pyenv virtual environments.

The above `PYTHONPATH` assumes your wcEcoli repo will be in `$HOME/wcEcoli`,
so it's handy to have an alias like this in your `.bashrc`:

```bash
alias cdw='cd ~/wcEcoli'
```

### `umask 0002` is needed before modifying $GROUP_HOME or $GROUP_SCRATCH

**NOTE:** Before installing pips in the virtual environments within `${PI_HOME}/pyenv`, or otherwise changing files or directories in `$GROUP_HOME` or `$GROUP_SCRATCH`, be sure that the `umask` command prints `0002`, or `umask -S` prints `u=rwx,g=rwx,o=rx`. This way, new and edited files won't be locked for the rest of the group.

The shared `"${PI_HOME}/etc/bash_profile"` does set `umask 0002` and aliases `mkdir="mkdir -m u=rwx,g=rwx,o=rx"`. It also adds `${PI_HOME}/modules` to the `MODULEPATH` to make the lab's shared modules available for loading.

### Install git hooks

Follow the instructions in wcEcoli's `runscripts/git_hooks/README.md` to set up git hooks in each of your wcEcoli development repos.

In particular, the `post-checkout` hook will diff the list of pips in the checked-out code (`requirements.txt`) against your installed pips to show if there are pip updates to install into your development environment. See the `pip` commands in `requirements.txt`.


## Environment Modules

The purpose of an _environment module_ is to make specific versions of executables and native libraries available by adding paths to environment variables such as `$PATH`, `$C_INCLUDE_PATH`, and `$LD_LIBRARY_PATH`. The module software can cleanly unload any environment module.

Example commands:
* `module list` -- list the loaded modules
* `module spider git` -- search for loadable versions of git
* `module load git/2.39.1` -- load a particular git version
* `ml` -- a shortcut for `module list`
* `ml git/2.39.1` -- a shortcut for `module load git/2.39.1`

After installing a native library or executable, make a module for it in the `${PI_HOME}/modules/` directory. A module can be written in lua, e.g. `openblas/0.3.26.lua`, or in the `Module1.0` text format, e.g. `jenkins/1.596.2`.


## Installing native libraries and executables

**Note:** Always use a **Sherlock compute node** to compile native libraries or install Python libraries so the compiler will detect the right CPU model.

Look in `$GROUP_HOME/downloads-sherlock2/` and `$GROUP_HOME/installation_notes/` for downloaded software packages and notes on recompiling them as needed to install new packages, new libraries, and new Python releases for the team.

When adding or updating libraries or executables (such as Python), please add the installation steps to `$GROUP_HOME/installation_notes/`. These range from simple like `mongodb` to not-so-simple like `openblas.txt` and `python-3-11-3.txt`.

Then add an environment module for the new library or library version.


## Installing new Python versions

**Note:** Always use a compute node to build Python so the compiler will detect the right CPU model.

Look in `$GROUP_HOME/downloads/` and `$GROUP_HOME/installation_notes/` for downloaded software packages and notes on recompiling them as needed to install new packages, new libraries, and new Python releases for the team.

We use `pyenv` to install versions of Python. Begin by updating `pyenv` itself:

```bash
cd $GROUP_HOME/pyenv
git fetch --all && git checkout master && git merge origin/master
```

Then list relevant versions of Python that are available, e.g.:

```bash
pyenv install -l | grep ' 3\.11\.'
```

Then follow steps similar to the previously installed version of Python.

The new version of Python might need updated/added native libraries. For a needed library, use `module spider <library>` to check if it's already available on Sherlock and just has to be added to the existing module (try `module show wcEcoli/python3`).


## Updating `wcEcoli3` or another shared pyenv virtualenv

When updating a shared pyenv virtualenv, plan to test it out both on your development computer and in a Jenkins PR build.

Once everything is working on your development computer and in Sherlock CI builds, plan to cut over to the new virtualenv in a short time window so it's less likely to break other builds happening at the same time.

* On Sherlock, remember to build and install all software on a **compute node** so it'll compile for the correct CPU model.
* Within the Pull Request, temporarily edit the file:  
  `runscripts/jenkins/setup-environment.sh`  
  changing the line:  
  `pyenv local wcEcoli3`  
  to instead select a staging virtualenv:  
  `pyenv local wcEcoli3-staging`
* Update this staging virtualenv (`wcEcoli3-staging`) on Sherlock with the updated pip versions. In a case like changing to a new Python version, you'll need to delete and rebuild the virtualenv. Use `pip freeze` to make a complete list of pips with "frozen" versions, merge that into `requirements.txt`, add it to the PR, then run `pip install -r requirements.txt` on Sherlock to update to that list.
* Watch for PR build success. Iterate as needed to fix problems and respond to code review feedback.
* When ready to merge the PR into the `master` branch:
  1. Update Python and libraries in virtualenv `wcEcoli3` to match `wcEcoli3-staging`.
  2. Revert the `pyenv local ...` edit in `setup-environment.sh`.
  3. Merge the PR into `master`.
