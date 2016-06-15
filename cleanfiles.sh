#!/bin/bash
for i in $(find horton tools scripts | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
(cd doc; make clean)
rm -v MANIFEST
rm -vr dist
rm -vr build
rm -vr doctrees
rm -v horton/cext.cpp
rm -v horton/*/cext.cpp
rm -v data/examples/hamiltonian/ring.xyz
rm -v data/examples/hamiltonian/*.h5
rm -v data/examples/hamiltonian/*FCIDUMP*
rm -v data/examples/hf_dft/*.h5
rm -v data/examples/hf_dft/*.molden
rm -v data/examples/hf_dft/*FCIDUMP*
rm -v data/examples/hf_compare/*.png
rm -v data/examples/ap1rog/*.h5
rm -v data/examples/ap1rog/*FCIDUMP*
rm -v data/examples/orbital_entanglement/*.h5
rm -v data/examples/orbital_entanglement/*.dat
rm -v data/examples/perturbation_theory/*.h5
rm -v data/examples/wpart/charges.txt
rm -v .coverage
exit 0
