Citing scientific work
======================

Whenever you add a feature based on a scientific publication, it should be cited
properly:

1. Add an item to the file ``data/references.bib``. Include the ``doi`` if
   possible. (The ``url`` field can be used as an alternative if the ``doi`` is
   not available.) Maintain the chronological order.

2. Add ``log.cite('someref', 'a reason')`` to the code based on the publication, e.g.\
   ``log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')``.

3. Update the references in ``doc/ref_literature.rst`` by running the script
   ``updateliterature.sh``.
