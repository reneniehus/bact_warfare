# bact_warfare

This folder contains MATLAB and R code to reproduce results from Niehus et al. "The Evolution of Strategy in Bacterial Warfare". 

In particular, the file invasion_analysis.m implements an invasion analysis according to classic game theory (see for example McElreath, R. & Boyd, R. Mathematical Models of Social Evolution.

Mathematical Models of Social Evolution (2013)). It assumes a monomorphic metapopulation where repeatedly rare mutants appear and it is tested whether these mutants can or cannot invade. This pairwise invasibility algorithm is efficient in finding the optimal consitutive strategy. It follows the pseudocode shown in Supplementary Methods and Results of Niehus et al. 

The function mass_battle.m implements the co-evoluationary tournament described in Niehus et al. and that leads to Figure 4. It allows for a polymorphic metapopulation where all three sensing strategies meet and compete in 1-on-1 competitions. Mutation, migration and selection is at work to establish the long-term winner of this evolutionary race. Note that this script can run up to hours, and that the control parameters can be tuned to optimise the time until convergence. 

