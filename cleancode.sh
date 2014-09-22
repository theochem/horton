#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in $(find data doc horton tools scripts examples *.py *.sh | egrep "(\.rst$)|(README)|(\.bib$)|(\.py$)|(\.c$)|(\.h$)|(\.nwchem)|(\.pyx$)|(\.pxd$)|(\.cpp)|(\.sh)"); do
  echo Cleaning ${file}
  dos2unix -q ${file}
  sed -i -e $'s/\t/    /g' ${file}
  sed -i -e $'s/[ \t]\+$//g' ${file}
  sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' ${file}
done
