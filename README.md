# RANN
# TODO
WEEK 1  
Write pseudo code for the random tranformation  
Pseudo-Code for P Q and F  
Explain the steps in the pseudo code and explain the math  
  
  
WEEK 2  
Implement the pseudo code! Do it in small chuncks and test that everything works  
Start by finding out how you create and run a futhark program  
  
WEEK 3  
Write issue, install cuda.  
Get access to A100. NEED HELP...  
Plan out the K-D trees section of the algorithm  
Start on related work  
Extra: describe the work/span of code from week 2.  
  
WEEK 4  
Construct the k-d tree  
call it-support for a100  

WEEK 5&6  
Make the functions for:  
    finding the leaf a point naturally falls into and its path  
    finding the paths that differ by 1.  
    bruteforceing the found paths.  
Ask troels to help with my installation issues for futhark 0.23/ look at it again myself  
  
WEEK 7  
Look at the minor adjustments to treeProcess cosmin talked about (Not too important)  
Fix the tree-building such that leaves does not nessasarily have the same size.  
    These trees are build from only knowing the height. computeTreeShape now unnesessary.  
    I need to have a shape array that says where the leaves starts/stops (See cosmins comments)  
    Adjust treeProcess so it uses this new tree  
Alternatively continue on the algorithm.  