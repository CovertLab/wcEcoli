Mac - how to set up the runtime environment for the model
===================================================

The `requirements.txt` file contains terse setup instructions for a Python environment, mainly the package list to use with `pip install`. See the additional details, below.

You need to have Xcode's command line tools to run the C compiler. Unless brew did that for you, run

    xcode-select --install

Tip: On macOS, iTerm2 is much nicer than the built-in Terminal app.

pyenv tools
------------

Use homebrew's `brew install` command to install

    pyenv
    pyenv-virtualenv
    pyenv-virtualenvwrapper
    glpk
    swig
    suite-sparse

swig is needed to install pip CVXcanon.
suite-sparse is needed for cvxopt to call glpk.

Set pyenv and optionally pyenv-virtualenv to initialize in your shell login script (`.profile` on macOS).

Example `.profile` lines for macOS:

    export PYENV_ROOT=/usr/local/var/pyenv
    if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi
    if which pyenv-virtualenv-init > /dev/null; then eval "$(pyenv virtualenv-init -)"; fi

Also consider (when you're not working on other Python projects):

    export PYTHONPATH="$HOME/wcEcoli:$PYTHONPATH"

