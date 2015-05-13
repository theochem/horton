Writing examples scripts
########################

The examples are located in the ``data/examples`` directory. The first level of
subdirectories is grouping the examples by category. Additional (small) input
files may be provided next to the example scripts. Make sure that each script
properly executes in a few seconds and that it does not require input that can't
be included in the examples directory (because it is too big or binary). Keep
the examples short and put excessive comments before each line to explain how
the example work. For each example, one should also add a corresponding test in
``horton/test/test_examples.py`` that just executes the example.
