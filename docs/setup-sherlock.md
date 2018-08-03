Sherlock - how to set up the runtime environment for the model
===================================================

The requirements.txt file contains terse setup instructions for a Python environment. See the additional details, below.

See $PI_HOME/downloads/ for notes on downloading and compiling native code packages like NumPy and glpk. This is useful to add new packages or new versions.

Setup Details
-------------

Example .bash_profile lines for Sherlock 2.0:

    module load wcEcoli/sherlock2

    export PYENV_ROOT="${PI_HOME}/pyenv"

    if [ -d "${PYENV_ROOT}" ]; then
        export PATH="${PYENV_ROOT}/bin:${PATH}"
        eval "$(pyenv init -)"
        eval "$(pyenv virtualenv-init -)"
    fi

Also consider:

    export PYTHONPATH="$HOME/wcEcoli:$PYTHONPATH"

pyenv, pyenv-virtualenv, pyenv-virtualenvwrapper, glpk, swig (needed to install pip CVXcanon), and suite-sparse (if needed for cvxopt to call glpk) are already available on Sherlock so you don't need to install these yourself.

Note: Before installing a new Python version on Sherlock, you might need to do [TODO: Should these be in module wcEcoli/sherlock2?]

    module load readline/7.0
    module load sqlite/3.18.0

or:

    export CPPFLAGS=-I/share/software/user/open/sqlite/3.18.0/include:/share/software/user/open/readline/7.0/include

to avoid

    WARNING: The Python readline extension was not compiled. Missing the GNU readline lib?
    WARNING: The Python sqlite3 extension was not compiled. Missing the SQLite3 lib?
    PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.15

[If running the "small" unit tests (below) fails with an error that the loader can't load .../pyenv/versions/.../lib/libpython2.7.a, that means you didn't successfully enable the shared library.]

