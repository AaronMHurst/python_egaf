#!/bin/bash

# Script to remove folders created during package-building process.
# Uninstalls pyEGAF package, then re-installs it.
# Run this script after editing source modules in the pyEGAF package.

rm -rf build
rm -rf dist
rm -rf pyEGAF.egg-info

rm -f *~
rm -f pyEGAF/*~
rm -f tests/*~

pip uninstall pyEGAF<<EOF
y
EOF

python setup.py install
