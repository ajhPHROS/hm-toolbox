# Hierarchical matrix toolbox

[![Build status](https://api.travis-ci.org/numpi/hm-toolbox.svg?branch=master)](https://travis-ci.org/numpi/hm-toolbox)

The <code>hm-toolbox</code> is a toolbox implementing the arithmetic of HODLR and HSS matrices in MATLAB. 

The HODLR case is handled in the @hodlr class, and correspond to H-matrices with partitioning
recursively done in 2 x 2 blocks, where the off-diagonal blocks are of low-rank. The HSS 
arithmetic uses the same partitioning (with nested bases), and is available through the
@hss class. 

Routines to compute matrix functions [1] and to solve matrix equations are included [1,2]. 

Chebfun2 is required for the construction of HODLR / HSS matrices which sample 
(piecewise) regular functions on a grid.

# Installation instructions

To install the toolbox download the latest revision from Github by running
```
 git clone https://github.com/numpi/hm-toolbox.git
```

or downloading the ZIP file from the webpage [github.com/numpi/hm-toolbox](https://github.com/numpi/hm-toolbox). 
Rename the folder to <code>hm-toolbox</code> if needed. Then, add it to your MATLAB path by running
```Matlab
 >> addpath /path/to/hm-toolbox
```

You are now ready to create new @hodlr and @hss objects. Check some examples in the
<code>examples/</code> folder. 

# References

1. Massei, S., Palitta, D., & Robol, L. (2018). Solving Rank-Structured Sylvester and Lyapunov Equations. SIAM Journal on Matrix Analysis and Applications, 39(4), 1564-1590.
2. Kressner, D., Massei, S., & Robol, L. (2017). Low-rank updates and a divide-and-conquer method for linear matrix equations. arXiv preprint arXiv:1712.04349.
