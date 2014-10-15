Writing examples scripts
########################

The examples are located in the ``data/examples`` directory. Each example is located
in a subdirectory ``XXX_some_name``, where ``XXX`` is a counter. In each
subdirectory there is an executable script ``run.py`` that demonstrates a
feature. Additional input files may be provided here for the example to work.
For each example, there is also a corresponding test in
``horton/test/test_examples.py`` Please, follow these conventions when adding a
new example.

The file ``run.py`` should not contain more than 50 lines of functional Python
code, preferably even less. Put excessive comments before each line to explain
how the example work. Also write the test such that it does not take more than a
second to complete.
