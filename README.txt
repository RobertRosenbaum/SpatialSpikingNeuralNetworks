
This folder contains Matlab and C code to generate Figures from the manuscript:

Rosenbaum, R., Smith, M. A., Kohn, A., Rubin, J. E., & Doiron, B. (2017). The spatial structure of correlated neuronal variability. Nature neuroscience, 20(1), 107.

Please cite this paper if you use any of the included software or any modifications of it.


The GPFA code in the subfolder gpfa-master was downloaded directly from: 
http://toliaslab.org/publications/ecker-et-al-2014/
and comes from the algorithms describes in the following papers:

Yu, B, et al. "Gaussian Process Factor Analysis for Low-Dimensional Single-Trial Analysis of Neural Population Activity." Journal of Neurophysiology 102.3 (2009).

Ecker, Alexander S., et al. "State dependence of noise correlations in macaque primary visual cortex." Neuron 82.1 (2014): 235-248.

Please cite these papers if you use the GPFA algorithms.

------

Instructions for use:

- Before running any of the Matlab code, first you must compile all of the C code using the mex compiler in Matlab.  This can be achieved by running the script CompileMexFiles.m  This will produce compiled code for the networks sims in the next step.  You only need to compile the C code once on each computer you run the code on.

- After compiling the C code, you can run the network simulations by running the scripts NetworkSimForFigureX.m where X is the Figure number.  This will produce large *.mat data files used to make the plots in the next step.  The *.mat files for Figures 3 and 4 are especially large (about 1GB) since they contain synaptic currents for 100 neurons (200 were used in the manuscript, but I reduced this number to make the files smaller).

- If you want to run all of the network simulations instead of just some of them, run the script RunAllNetworkSims.m.  This will take some time to complete.

- For Figure 4, you will also need to run TheoreticalCorrsForFigure4.m which computes the theoretical values of the correlations (the dashed curve).

- After running each simulation, you can generate the figures by running the scripts MakeFigureX.m.  You need to run the network simulation only once before you can generate the associated figure.


NOTES: 

- You must run NetworkSimForFigure6LGNtoL4 before NetworkSimForFigure6L4toL23 since NetworkSimForFigure6L4toL23 uses the spike trains generated in NetworkSimForFigure6LGNtoL4.

- The network simulations take some time to run, particularly Figures 3, 4 and 6 which take on the order of 10 minutes to run. 

- The network simulations use on the order of 6GB of memory while running.

- Figures 3 and 4 produce .mat data files that is on the order of 1GB in size.

- For questions about the code, please e-mail Robert Rosenbaum (currently at Robert.Rosenbaum@nd.edu)
