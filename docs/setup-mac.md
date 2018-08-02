Mac - how to set up the runtime environment for the model
===================================================

The requirements.txt file contains terse setup instructions for a Python environment. See the additional details, below.

You need to have Xcode installed, then run

    xcode-select --install

(install its command line tools) in order to run the C compiler.

On macOS, iTerm2 is much nicer than the built-in Terminal app.

pyenv tools
------------

Install pyenv, pyenv-virtualenv, pyenv-virtualenvwrapper, glpk, swig (needed to install pip CVXcanon), and suite-sparse (if needed for cvxopt to call glpk). These are available through homebrew.

Set pyenv and optionally pyenv-virtualenv to initialize in your shell login script (.profile on macOS).

Example .profile lines for macOS:

    export PYENV_ROOT=/usr/local/var/pyenv
    if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi

Also consider:

    export PYTHONPATH="$HOME/wcEcoli:$PYTHONPATH"

