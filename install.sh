#!/bin/bash

# this literally just ensures that FoBiS is installed
# and then clones json-fortran (required dependency) and builds it

JSONDIR='json-fortran'
REPO="https://github.com/jacobwilliams/${JSONDIR}"
PACKAGE='FoBiS.py'

pip install $PACKAGE
git clone $REPO
cd $JSONDIR
./build.sh
