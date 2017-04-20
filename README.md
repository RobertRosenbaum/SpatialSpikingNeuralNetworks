# SpatialSpikingNeuralNetworks
C and Matlab code for simulatung large, spatially extended spiking neuronal networks

This repository contains Matlab and C code to generate Figures 3 and 4 from the manuscript:

Rosenbaum, R., Smith, M. A., Kohn, A., Rubin, J. E., & Doiron, B. (2017). The spatial structure of correlated neuronal variability. Nature neuroscience, 20(1), 107.

Please cite this paper if you use any of the included software or any modifications of it.

Instructions for use:
- Before running any of the Matlab code, first you must compile all of the C code using the mex compiler in Matlab.  Specifically, in the Matlab command line, run the following commands:
mex EIF2DSpatialNetworkNoJitter.c
mex CalcRheoBaseEIF.c
This will produce compiled code for the networks sims in the next step.

- After compiling the C code, you can run the network simulations by running the scripts NetworkSimForFigureX.m where X=3,4 is the Figure number.  This will produce large *.mat data files used to make the plots in the next step.

- After running each simulation, you can generate the figures by running the scripts MakeFigureX.m.  You need to run the network simulation only once before you can generate the associated figure.

NOTE: The network simulations take on the order of 10 minutes to run and they each produce a .mat data files that is on the order of 1GB in size.

- For questions about the code, please e-mail Robert Rosenbaum (currently at Robert.Rosenbaum@nd.edu)
