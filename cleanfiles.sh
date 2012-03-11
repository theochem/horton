#!/bin/bash
for i in $(find horton | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
rm -vr doc/_build/
rm -v MANIFEST
rm -vr dist
rm -vr build
rm -v horton/cext.c
