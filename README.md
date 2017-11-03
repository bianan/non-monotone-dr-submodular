# Non-monotone Continuous DR-submodular Maximization: Structure and Algorithms


This repository collects source code for the paper:

"Non-monotone Continuous DR-submodular Maximization: Structure and Algorithms"

NIPS 2017. An Bian, Kfir Y. Levy, Andreas Krause and Joachim M. Buhmann

### File Structure:

 - src/:  contians source files
 - main.m: the main file to run different experiments

## Usage:

See the setup guide inside the main file "main.m" to run different
experiments.


### Dependencies:

 - This implementation uses the quadprogIP solver provided by  [quadprogIP](https://github.com/xiawei918/quadprogIP), the MATLAB LP solver [linprog](https://ch.mathworks.com/help/optim/ug/linprog.html) and the
 MATLAB quadratic programming solver [quadprog](https://ch.mathworks.com/help/optim/ug/quadprog.html).

 - Additionally, the quadprogIP solver requires the Matlab interface to
 CPLEX 12.2 or later.

 - The code has been tested on Ubuntu 16.04 LTS, 64 bits with MATLAB R2016a. It should work with other OS with little change.

### Copyright:

 Copyright (2017) [An Bian | ETH Zurich | http://neocortex.ch/].
 Please cite the above paper if you are using this code in your work.
