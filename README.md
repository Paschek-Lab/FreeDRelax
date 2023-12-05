# DrelaxD: Computing the Frequency-Dependent NMR Relaxation of <sup>1</sup>H Nuclei in Liquid Water

This repository contains a collection of input files and source code as described in the paper "Computing the Frequency-Dependent NMR Relaxation of <sup>1</sup>H Nuclei in Liquid Water", presenting a computational framework for reliably determining the frequency dependent intermolecular NMR dipole-dipole relaxation rate of spin 1/2 nuclei from MD simulations.

A preprint of the paper is available on XXX.

The repsitory is structured as follows:
* `MD`: contains the input parameter and topology files for the MD simulations with system sizes of 512, 1024, 2048, 4096 and 8192 molecules. The simulations were performed with the GROMACS 5.0.6 software package. All parameters different from default can be found in `SIMXX.mdp` while `sim1out.mdp` include all used parameters. For naming convention, please refer to the Gromacs-5.0.6-manual. The start configurations are stored in `START.gro`, the corresponding force field parameters can be found in `topol.top`. 
  
* `MC`: contains the code for performing randm walker Monte Carlo simulations
  
* `relaxation`: contains the code for computing the intra- and intermolecular relaxation rates 


For questions, please contact [the authors](mailto:dietmar.paschek@uni-rostock.de)