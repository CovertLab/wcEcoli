Create a pyenv
================

Install the relevant version of Python in pyenv if it's not already there, enabling shared (dynamic) libraries.

    PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.15

Create a pyenv and use it for the local directory.

    cd ~/wcEcoli  # or wherever your `wcEcoli` development directory is
    pyenv local 2.7.15
    pyenv virtualenv wcEcoli2  # <-- picking a new virtualenv name
    pyenv local wcEcoli2

Upgrade the installers.

     pip install --upgrade pip setuptools virtualenv wheel

Pip needs a local numpy before installing the other packages:

    pip install numpy==1.14.5

Then install all the other packages listed in requirements.txt:

    CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary cvxopt
    pyenv rehash

NOTE: Per Issue #161, the backend for matplotlib might need to be changed from `TkAgg` to the default (or specifically to `agg`) in `$PYENV_ROOT/versions/wcEcoli2/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc`.

   * This is only a problem when running under Fireworks.
   * All the wcEcoli code expects to run in the `wcEcoli` directory and with that directory on the `PYTHONPATH`, but Fireworks `rlaunch rapidfire` sets the current working directory to its launcher_... subdirectory instead. `rlaunch singleshot` does not have that problem, manual runscripts do not have that problem, and we didn't test `qlaunch`.
   * After installing or updating matplotlib, edit the .../site-packages/matplotlib/mpl-data/matplotlibrc file, commenting out the line `backend : TkAgg` _OR_ first test if it can import pyplot:
      1. cd to a directory that does not have a matplotlibrc file, e.g. `wcEcoli/docs/`.
      2. Use the `pyenv version` command to verify that the new pyenv is active. (This should be fine if the current directory is a subdirectory of wcEcoli and if the latter has a `.python-version` file set via the `pyenv local` command.)
      3. Start a python shell and type `import matplotlib.pyplot`
      4. If it raised `ImportError: No module named _tkinter`, then you need to edit the .../site-packages/matplotlib/mpl-data/matplotlibrc file, commenting out the line `backend : TkAgg`
      5. If it didn't raise an error, run `matplotlib.get_backend()` and check that it returns `'agg'` or similar.
      6. Retest after adjusting the matplotlibrc configuration.

Run the "small" unit tests.

    nosetests -a smalltest
