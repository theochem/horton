#!/bin/bash
git checkout master && (
  ./cleanfiles.sh
  (
    cd doc
    make html
    cd _build/html
    mv _static static
    for f in $(find); do sed -e 's/_static/static/' -i $f; done
    mv _sources sources
    for f in $(find); do sed -e 's/_sources/sources/' -i $f; done
  )
  git checkout gh-pages && (
    git rm $(git ls-tree $(git log --color=never | head -n1 | cut -f2 -d' ') -r --name-only)
    cp -rv doc/_build/html/* .
    for f in $(find doc/_build/html | cut -c17-); do
      
      git add $f
    done
    git commit -a -m 'Automatic documentation update'
    git push official gh-pages:gh-pages
  )
) && git checkout master
