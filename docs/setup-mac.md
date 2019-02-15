Mac - how to set up the runtime environment for the model
===================================================

The `requirements.txt` file contains terse setup instructions for a Python environment, mainly the package list to use with `pip install`. See the additional details, below.

You need to have Xcode's command line tools to run the C compiler. Unless brew did that for you, run

    xcode-select --install

Tip: On macOS, iTerm2 is much nicer than the built-in Terminal app.

pyenv tools
------------

Use homebrew to install

    brew install pyenv pyenv-virtualenv pyenv-virtualenvwrapper glpk swig suite-sparse zlib readline xz

swig is needed to install pip CVXcanon.
suite-sparse is needed for cvxopt to call glpk.

When using Mojave or higher you will aslo need to install the additional SKD headers.

    sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /g

Set pyenv and optionally pyenv-virtualenv to initialize in your shell login script (`.profile` on macOS).

Example `.bash_profile` lines for macOS:

    export PYENV_ROOT=/usr/local/var/pyenv
    if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi
    if which pyenv-virtualenv-init > /dev/null; then eval "$(pyenv virtualenv-init -)"; fi
    
Run your .bash_profile

    source ~/.bash_profile

Also consider (when you're not working on other Python projects):

    export PYTHONPATH="$HOME/wcEcoli:$PYTHONPATH"
    
Now that you have pyenv and related libraries installed, you can install python and set up your local environment (NOTE: if you have a `~/.local/` directory, paths might not work properly with `pyenv` and you might receive error messages):

    PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.15
    pyenv local 2.7.15
    pyenv virtualenv wcEcoli2
    pyenv local wcEcoli2
    pip install --upgrade pip setuptools virtualenv wheel
    pip install numpy==1.14.3

Finally, run the `requirements.txt`:

    CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary cvxopt
    pyenv rehash

Post-finally, compile the CYTHON code:

    make clean compile

To make sure everything is working, run the tests:

    nosetests -a smalltest

If these pass, you are good to go.

