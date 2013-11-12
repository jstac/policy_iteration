The main program is main.m which solves for the optimal policy
using policy function iteration. It calls upon the following subroutines

1. ComputePolicy: find the fixed point of the FOC give a guess of the policy function

2. CalculateRHSEuler: evaluates the expected marginal utility 

3. interpne: linear interpolation on the grid and extrapolate outside the grid using the value at the
	     closest grid point 

4. u1: computes marginal utility 

5. InverseMargU: computes the consumption level from the price function

6. makegrid: creates a grid with user-specified degree of density at the low end of the grid

ComputeStatistics.m: produce plots and simulates shocks to compute the insurability measure 