# Mac - how to set up the runtime environment for the model

The `requirements.txt` file contains terse setup instructions for a Python environment, mainly the package list to use with `pip install`. See the additional details, below.

**NOTE:** The following instructions don't include installing OpenBLAS and then compiling NumPy and SciPy to use it. See [Create a pyenv](create-pyenv.md) for those instructions.

## Prerequisites

* See [Development tools](dev-tools.md) for instructions to set up the development tools.
* You need to have Xcode's command line tools to run the C compiler. Unless brew did that for you, run:

    xcode-select --install


## pyenv tools

Use homebrew to install:

    brew install pyenv pyenv-virtualenv pyenv-virtualenvwrapper
    brew install glpk openssl readline swig suite-sparse xz

When using Mojave or higher (10.14+) you will also need to install the additional SDK headers:

    sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /

Set `pyenv` and optionally `pyenv-virtualenv` to initialize in your shell login script (either `.bash_profile` or `.profile` on macOS):

    export PYENV_ROOT=/usr/local/var/pyenv
    if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi
    if which pyenv-virtualenv-init > /dev/null; then eval "$(pyenv virtualenv-init -)"; fi
    ## -- Do this before sourcing iterm2_shell_integration
    
Open a new shell so it runs the updated profile.

Also consider (when you're not working on other Python projects):

    export PYTHONPATH="$HOME/wcEcoli:$PYTHONPATH"

(NOTE: if you have a `~/.local/` directory, paths might not work properly with `pyenv` and you might receive error messages. [TODO] In that case, delete the directory?)

Now that you have pyenv and related libraries installed, you can install python and set up your local environment:

    PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.15
    cd ~/dev/wcEcoli  # or wherever you cloned `wcEcoli` to
    pyenv local 2.7.15
    pyenv virtualenv wcEcoli2
    pyenv local wcEcoli2
    pip install --upgrade pip setuptools virtualenv wheel

Finally, run the `requirements.txt`:

    CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary cvxopt
    pyenv rehash

Post-finally, compile the native code:

    make clean compile

To make sure everything is working, run the tests:

    nosetests

If these pass, you are good to go.

