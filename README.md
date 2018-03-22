Event Reweighting Tests
=======================
Authors: A. Lister, A. Mastbaum

These scripts check that reweighting of events in an analysis tree, via the
TreeReader input source, produces exactly the same as running EventWeight
directly on art ROOT files.

The `TestModule` module produces an example analysis tree with all the truth
information required for reweighting. Feel free to adapt this code for
writing analysis-specific output trees. Note that the branch names can be
adjusted to your preferences (the mapping is specified when the analysis trees
are read back in for reweighting), to minimize need to duplicate information
and/or adapt existing analysis code.

N.B. This code assumes uboonecode `v06_26_01_12`.

Steps
-----

1. Reweight some MC of your choice.

   For example:

    ```
    lar -c run_eventweight_microboone.fcl /pnfs/uboone/mc/uboone/detector-simulated/prodgenie_bnb_nu_uboone/sim/prod_v06_26_00/prodgenie_bnb_nu_uboone_0_20170228T015121_gen2_0aa93a49-868a-41dd-a795-e3ef4643f5b9.root
    mv <that output> mcc8_test.root
    ```

2. Turn that into an analysis tree with `TestModule`.

   Place `TestModule` into the uboonecode source, add it to the CMakeLists,
   and rebuild uboonecode. Then:

    ```
    lar -c run_test_module.fcl mcc8_test.root
    ```

3. Read the tree back in and re-run eventweight (same seeds).

    ````
    lar -c test_treereader.fcl testOutput.root
    ```

4. Check that all the weights match.

    ```
    make  # Build the test program
    ./check_weights mcc8_test.root treereader.root
    ```

   if any weights mismatch, it will print an error and abort immediately.

