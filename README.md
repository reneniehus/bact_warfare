This folder contains MATLAB and R code to reproduce results from Niehus et al. "The Evolution of Strategy in Bacterial Warfare".

In particular, the file invasion_analysis.m implements an invasion analysis according to classic game theory (see for example McElreath, R. & Boyd, R. Mathematical Models of Social Evolution. Mathematical Models of Social Evolution (2013)). It assumes a monomorphic metapopulation where repeatedly rare mutants appear, and for each mutant it is tested whether it can or cannot invade. This pairwise invasibility algorithm is efficient in finding the optimal constitutive strategy. It follows the pseudocode shown in Supplementary Methods and Results of Niehus et al.

The function mass_battle.m implements the co-evolutionary tournament described in Niehus et al. (see Methods) and that leads to Figure 4. It allows for a polymorphic metapopulation where all three sensing strategies meet and compete in 1-on-1 competitions. Mutation, migration and selection at the meta-population level is at work to establish the long-term winner of this evolutionary race. Note that this script can run up to hours, and that the control parameters can be tuned to improve convergence speed.

The R file vis_mass_batt_fig4c.R loads the results of mass_battle.m and plots them just like in Figure 4 c, allowing to inspect convergence and dominating strategies.
