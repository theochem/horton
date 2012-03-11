#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in $(find data doc horton *.py | egrep "(\.rst$)|(\.py$)|(\.c$)|(\.h$)|(\.nwchem)|(\.pyx$)|(\.pxd$)"); do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /g' ${file}
  sed -i -e $'s/[ \t]\+$//g' ${file}
  sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' ${file}
done
