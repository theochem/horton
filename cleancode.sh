#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
for file in $(find doc horton *.py | egrep "(\.rst$)|(\.py$)|(\.c$)|(\.h$)"); do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
  sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' ${file}
done
