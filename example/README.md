Steps to run conta example (scripts need to be run in the example directory):

1) ./run_small.sh -- if it works (if conta R library and dependencies are correctly installed)
Check the output files under test directory to get familiar with conta output.

2) ./run.sh -- takes a while. This is an experiment of cross contamination titrations of 0.8% to
0.01%. Each sample is separately run through conta.

3) ./combine.sh -- combines conta.tsv files under each output folder to a single results.tsv 

The sample with "_1" extension is actually pure (not contaminated on purpose). This should
be the only sample not being conta called. 

4) ./clean.sh -- clean up once you are done
