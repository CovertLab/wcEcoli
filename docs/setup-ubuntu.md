Ubuntu - how to set up the runtime environment for the model
===================================================

There are a number of dependencies that need to be set up before the model will run on Ubuntu.

    sudo apt install -y gcc make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev glpk-utils libglpk-dev glpk-doc python-cvxopt

You may also have to find the proprietary package `python-glpk` and install that as well. 

Also, set your `PYTHONPATH` to point to the WcEcoli model directory:

    export PYTHONPATH=$CODE/WcEcoli

pyenv
-----

pyenv has a related set of packages that need to be installed, which on ubuntu is currently through cloning git repositories and adding things to your `.bash_profile`.

* pyenv
* pyenv-virtualenv
* pyenv-virtualenvwrapper

To set these up there are a number of commands to run:

    git clone https://github.com/pyenv/pyenv.git ~/.pyenv
    echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile
    echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile
    echo -e 'if command -v pyenv 1>/dev/null 2>&1; then\n  eval "$(pyenv init -)"\nfi' >> ~/.bash_profile

    git clone https://github.com/pyenv/pyenv-virtualenv.git $(pyenv root)/plugins/pyenv-virtualenv
    echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.bash_profile

    git clone https://github.com/pyenv/pyenv-virtualenvwrapper.git $(pyenv root)/plugins/pyenv-virtualenvwrapper
    source ~/.bash_profile

Now that you have pyenv and related libraries installed, you can install python and set up your local environment:

    pyenv install 2.7.15
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

    nosetests -a smalltests

If these pass, you are good to go.