# Create a pyenv

The [requirements.txt](https://github.com/CovertLab/wcEcoli/blob/master/requirements.txt) file contains terse setup instructions for a Python environment, mainly the package list to use with `pip install`. This page contains additional details.

1. On Sherlock, see `$PI_HOME/downloads/` for notes on downloading and compiling native code packages like `glpk`, as needed to install new packages or new versions.
2. See [Development-Tools](Development-Tools) for instructions to install the development tools including a C compiler and `pyenv`. Those steps include using your package manager to install `pyenv`, `pyenv-virtualenv`, and `pyenv-virtualenvwrapper`.
3. Use your package manager (e.g. **homebrew** on macOS) to install the needed packages [see the `requirements.txt` file for the current list], e.g.:

   ```bash
   brew install glpk openssl swig suite-sparse xz
   ```

## Setup mods to fix Non-determinism Issues, Jan. 8, 2019

**Problem:** Per [Issue #196](https://github.com/CovertLab/wcEcoli/issues/196), the Fitter does not always produce consistent results from run to run. The following modifications to the setup instructions seem to fix that by building NumPy and SciPy on OpenBLAS v0.3.5+.

**Background:** You can install different implementations of BLAS (Basic Linear Algebra Subprograms) and link NumPy and SimPy to call them. OpenBLAS 0.3.5 fixes several bugs with non-determinism and hanging.

   * Linux might come with an unoptimized version of BLAS.
   * MacOS comes with Apple's Accelerate. As of Dec. 5, 2018 OS 10.14.1 Mojave, it exhibits non-determinism bugs. See [#196](https://github.com/CovertLab/wcEcoli/issues/196).
   * Intel provides [Math Kernel Library (MKL) and related software](https://software.intel.com/en-us/articles/installing-the-intel-distribution-for-python-and-intel-performance-libraries-with-pip-and) optimized for their CPUs. It seems to run a little faster than Accelerate and OpenBLAS, and their `tbb` tool should fix the "oversubscription" problem where multiple processes collectively request too many CPU cores. **But** as of Dec. 5, 2018, it exhibits the non-determinism bugs. Eventually it'll become a good choice for speed and easier installation, in brief, `pip uninstall numpy scipy scikit-learn -y; pip install intel-numpy intel-scipy tbb4py`.
   * OpenBLAS [release](https://github.com/xianyi/OpenBLAS/releases) `0.3.5` fixed several multi-threading bugs.

**Installation Steps:** For OpenBLAS, NumPy, and SciPy.

   1. Fetch OpenBLAS source, compile it, and install it in a local directory [it does not usually belong on the compiler include path and the linker library path] or install 0.3.5+ using a package manager. [On Sherlock, it's in `$PI_HOME/downloads-sherlock2/compiled/openblas/lib`.]:

      ```bash
      brew install gcc  # if you don't have a gfortran compiler
      cd <the directory containing your git projects>
      git clone https://github.com/xianyi/OpenBLAS
      cd OpenBLAS
      git checkout v0.3.5
      make FC=gfortran
      sudo make PREFIX=/opt/OpenBLAS install  # <-- pick another PREFIX dir if you don't/can't sudo
      ```

   2. Download https://github.com/numpy/numpy/blob/master/site.cfg.example to your local file `~/.numpy-site.cfg`
   3. Uncomment and edit this part of the openblas section (_naming the directory you installed it into_) (or just save these lines; the rest of the downloaded `site.cfg.example` file is just instructions and commented-out configuration choices):

      ```
      [openblas]
      libraries = openblas
      library_dirs = /opt/OpenBLAS/lib
      include_dirs = /opt/OpenBLAS/include
      ```

   4. Install NumPy and SciPy from source (without manually downloading them) in your pyenv using your openblas build:

      ```bash
      cd wcEcoli
      pip install numpy==1.14.5 scipy==1.0.1 --no-binary numpy,scipy --force-reinstall
      ```

      _Or_ combine that step with installing all the other required pips by running:

      ```bash
      CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary numpy,scipy,cvxopt
      ```

      and if you had a previously-installed theano:

      ```bash
      theano-cache clear
      ```

   5. Check the installation

      ```bash
      python runscripts/debug/summarize_environment.py
      ```
      It should print something like this for numpy and scipy,  referring to the
      `library_dirs` that you set above:
      ```
      lapack_opt_info:
          libraries = ['openblas', 'openblas']
          library_dirs = ['/opt/OpenBLAS/lib']
          define_macros = [('HAVE_CBLAS', None)]
          language = c
      ```

      ```bash
      python
      import theano
      theano.config.blas.ldflags
      ```
      It should print something like `-L/opt/OpenBLAS/lib -lopenblas -lopenblas`, referring to the
      library_dirs that you set above.


## Create a pyenv

1. **Sherlock only:** Before installing another version of Python on Sherlock, you might need to do these steps [TODO: Put them in module `wcEcoli/sherlock2`?]:

   ```bash
   module load readline/7.0
   module load sqlite/3.18.0
   ```

   or:

   ```bash
   export CPPFLAGS=-I/share/software/user/open/sqlite/3.18.0/include:/share/software/user/open/readline/7.0/include
   ```

   to avoid:

   WARNING: The Python readline extension was not compiled. Missing the GNU readline lib?
   WARNING: The Python sqlite3 extension was not compiled. Missing the SQLite3 lib?

2. Install the relevant version of Python via `pyenv`, enabling shared (dynamic) libraries.

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.15
   ```

   [If running the unit tests (below) fails with an error message saying the loader can't load `.../pyenv/versions/.../lib/libpython2.7.a`, that means you didn't successfully `--enable-shared` when installing python.]

3. Create a pyenv and use it for the local directory.

   ```bash
   cd ~/wcEcoli  # or ~/dev/wcEcoli or wherever you cloned `wcEcoli` to
   pyenv local 2.7.15
   pyenv virtualenv wcEcoli2  # <-- picking a suitable virtualenv name
   pyenv local wcEcoli2
   ```

4. Upgrade the installers.

   ```bash
   pip install --upgrade pip setuptools virtualenv wheel
   ```

5. Then install all the other packages listed in `requirements.txt`:

   ```bash
   CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary cvxopt
   pyenv rehash
   ```
   - **NOTE:** Per [#161](https://github.com/CovertLab/wcEcoli/issues/161), the backend for `matplotlib` might need to be changed from `TkAgg` to the default (or specifically to `agg`) in `$PYENV_ROOT/versions/wcEcoli2/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc`.
      1. This is only a problem when running under Fireworks. All the wcEcoli code expects to run in the `wcEcoli` directory and with that directory on the `PYTHONPATH`, but Fireworks `rlaunch rapidfire` sets the current working directory to its `launcher_...` subdirectory instead. `rlaunch singleshot` does not have that problem, manual runscripts do not have that problem, and I don't know about `qlaunch`.
      1. After installing or updating matplotlib, test if it can import pyplot:
      1. `cd` to a directory that does not have a `matplotlibrc` file, e.g. `wcEcoli/docs/`.
      1. Use the `pyenv version` command to verify that the pyenv `wcEcoli2` is active. (This should be fine if the current directory is a subdirectory of `wcEcoli` and if the latter has a `.python-version` file set via the command `pyenv local wcEcoli2`.)
      1. Start a python shell and type `import matplotlib.pyplot`
         * If it raised `ImportError: No module named _tkinter`, then edit the `.../site-packages/matplotlib/mpl-data/matplotlibrc` file, commenting out the line `backend : TkAgg`, then retest.
         * If it didn't raise an error, run `matplotlib.get_backend()` and check that it returns `'agg'` or similar.

6. Compile the project's native code.

   ```bash
   make clean compile
   ```

(Yes, that prints deprecation warnings.)

7. Run the unit tests.

   ```bash
   nosetests
   ```
