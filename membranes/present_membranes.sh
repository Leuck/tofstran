#!/usr/bin/bash
# runs the three membrane examples from tofstran
PYBIN=/usr/bin/python3
$PYBIN testmembrane.py &
$PYBIN testmembrane2.py &
$PYBIN testmembrane3.py &

