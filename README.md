# A function to estimate the total stopping time distribution of the Collatz problem based on the stochastic approach.

## Details
This repository contains code implementing a function that approximate the total stopping time distribution for the Collatz function [1] based on the Brownian motion.
The Collatz function is 
[tex: f(n)=\frac n2 \quad if n 0]

We derived this distribution function in a similar approach to works [2][3].
The max of the total stopping time can be estimated by using the distribution function.
##File structure
|folder name|explanation                         |
|:--        |:--                          |
|doc        | Mathematical details on this approach.  |
|results    |Comparison results of Collatz sequences and estimation by using the distribution function|
|src        | A Code implementing a distribution function.|

The results folder also containts
-comparison of the total stopping time distribution of the Collatz sequences (3-point moving average) and estimation results for numbers less than 10^6.
-comparison result the max of the total stopping time of the Collatz sequences [4] and the prediction results less than 10^40.
## License
This project is licensed under the terms of the [MIT license](LICENSE.md).

## References
-  [1] Lagarias, Jeffrey C., ed. *The ultimate challenge: The 3x+ 1 problem*. American Mathematical Soc., 2010.
-  [2] Kontorovich, Alex V., and Jeffrey C. Lagarias. "Stochastic Models for the 3x+ 1 and 5x+ 1 Problems." *arXiv preprint arXiv:0910.1944 (2009).*
-  [3] Borovkov, Konstantin Aleksandrovich, and Dietmar Pfeifer. "Estimates for the Syracuse problem via a probabilistic model." *Theory of Probability & Its Applications* 45.2 (2001): 300-310.
-  [4] Roosendaal, Eric. "3x+1 Delay Records". http://www.ericr.nl/wondrous/delrecs.html

