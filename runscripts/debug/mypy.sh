#!/bin/sh
# Run the mypy type checker on the wcEcoli sources (except prototypes, etc.)
# The mypy.ini file sets some configuration options. So far it only runs a
# Python 2.7 check with Python 3 to come. A key goal is to find unicode/bytes
# mismatches in Python 3.
#
# ASSUMES: The current working directory is your wcEcoli git directory.
#
# ASSUMES: You created a `mypy` virtualenv like this:
#   pyenv install 3.8.2
#   pyenv shell 3.8.2
#   pyenv virtualenv mypy
#   pyenv shell mypy
#   pip install mypy
#   pyenv shell --unset

# Run in the Python 3 virtualenv `mypy` with access to the current virtualenv,
# which is typically the Python 2.7 virtualenv `wcEcoli2`.
export PYENV_VERSION="mypy:$(pyenv version-name)"

mypy --py2
