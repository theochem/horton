Citing scientific work
======================

Whenever you add a feature based on a scientific publication, it should be cited
properly:

1. Add an item to the file ``data/references.bib``. Include the ``doi`` if
   possible. (The ``url`` field can be used as an alternative if the ``doi`` is
   not available.) Maintain the chronological order. Do not use reference
   manager software to edit this file. Use a plain text editor instead and don't
   touch existing references if not needed.

2. Add ``log.cite('someref', 'a reason')`` to the code based on the publication, e.g.\
   ``log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')``.
   This guarantees that paper is properly cited at the end of the Horton screen
   output.
