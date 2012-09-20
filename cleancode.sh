#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in $(find data doc horton tools *.py | egrep "(\.rst$)|(README)|(\.bib$)|(\.py$)|(\.c$)|(\.h$)|(\.nwchem)|(\.pyx$)|(\.pxd$)|(\.cpp)"); do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /g' ${file}
  sed -i -e $'s/[ \t]\+$//g' ${file}
  sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' ${file}
done
