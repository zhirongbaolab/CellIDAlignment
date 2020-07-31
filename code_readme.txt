high level overview

test_em_embryo[1,2,9]_alignment_new.m are drivers that perform alignment on an individual embryo performing a pre-alignment of the unknown sample to the atlas model based on a subset of manually named cells (specific to the sample and specified in the script)


matchPointCloundsConsensusMultiembOptimizedN.m
matchPointCloundsConsensusMultiembOptimizedIteratedN.m

are the main meta alignment functions. These are provided with a set of coaligned named embryo data, expected neighbor constraint matricies, and in the case of iterated driver the set of confident matches from the first non iterated alignment.

WG_prealignment_model_driver.m provided for reference loads a set of acetree files representing fully edited example data, coaligns them to eachother using names and prepares lists of nearby divisions, expected neighbor constraint graphs etc. This script generates the general model over time in Global_Model_Fully_Loaded.mat 