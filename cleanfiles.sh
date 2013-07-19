#!/bin/bash
for i in $(find horton tools scripts | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
(cd doc; make clean)
rm -v MANIFEST
rm -vr dist
rm -vr build
rm -vr doctrees
rm -v horton/cext.cpp
rm -v horton/gbasis/cext.cpp
rm -v horton/grid/cext.cpp
rm -v horton/meanfield/cext.cpp
rm -v horton/espfit/cext.cpp
