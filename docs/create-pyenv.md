Create a pyenv
================

Note: When using pyenv and pip to install Python packages, note that some options may be needed like PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.15 and CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary cvxopt. It's not completely clear when to use them.

Install the relevant version of Python in pyenv if it's not already there, enabling shared (dynamic) libraries.

Create a pyenv and use it for the local directory.

    cd ~/wcEcoli  # or wherever your `wcEcoli` development directory is
    pyenv local 2.7.15
    pyenv virtualenv wcEcoli2  # <-- picking a new virtualenv name
    pyenv local wcEcoli2

Upgrade the installers.

     pip install --upgrade pip setuptools virtualenv wheel

NOTE: A plain pip install might now suffice for scipy.
TODO: Try using conda to install packages, esp. scipy, Intel MKL, and other binary packages. Also try the Intel Distribution for Python.
Then install all the other packages listed in requirements.txt:

    CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary cvxopt
    pyenv rehash

NOTE: Per #161, the backend for matplotlib might need to be changed from TkAgg to agg in $PYENV_ROOT/versions/wcEcoli2/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc.

This is only a problem when running under Fireworks. All the wcEcoli code expects to run in the wcEcoli directory and with that directory on the PYTHONPATH, but Fireworks rlaunch rapidfire sets the current working directory to its launcher_... subdirectory instead. rlaunch singleshot does not have that problem, manual runscripts do not have that problem, and we didn't test qlaunch.

After installing or updating matplotlib, test if it can import pyplot:

* cd to a directory that does not have a matplotlibrc file, e.g. wcEcoli/docs/.
* Use the pyenv version command to verify that the pyenv wcEcoli2 is active. (This should be fine if the current directory is a subdirectory of wcEcoli and if the latter has a .python-version file set via the command pyenv local wcEcoli2.)
* Start a python shell and type import matplotlib.pyplot
* If it raised ImportError: No module named _tkinter, then edit the .../site-packages/matplotlib/mpl-data/matplotlibrc file, commenting out the line backend : TkAgg
* If it didn't raise an error, run matplotlib.get_backend() and check that it returns 'agg' or similar.
* Retest after adjusting the matplotlibrc configuration.

Run the "small" unit tests.
    nosetests -a smalltest
