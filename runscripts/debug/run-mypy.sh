#!/bin/sh
# Run the mypy type checker on the wcEcoli sources (except prototypes, etc.)
#
# The mypy.ini config file sets `python_version = 2.7` etc. so this script
# doesn't need the --py2 option.
#
# This script will get more interesting when it needs to type-check the code
# for Python 3 or both.
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

pyenv shell mypy wcEcoli2
mypy models reconstruction runscripts validation wholecell
pyenv shell -
