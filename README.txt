This repository contains the Granger Causality (GC) and Synchrony analysis scripts and codes for the results reported in: 
Time-dependent homeostatic mechanisms underlie Brain-Derived Neurotrophic Factor action on neural circuitry

K. M. O'Neill, E. D. Anderson, S. Mukherjee, S. Gandu, S. A. McEwan, A. Omelchenko, A. R. Rodriguez, W. Losert, D. F. Meaney, B. Babadi, B. L. Firestein*

*Corresponding Author:
Bonnie L. Firestein: firestein@biology.rutgers.edu

Contact information regarding these codes:
Shoutik Mukherjee: smukher2@umd.edu
Behtash Babadi: behtash@umd.edu

MEA spike time recordings are available at Figshare:
https://doi.org/10.6084/m9.figshare.24596121

BDNF Dose Response Scripts:
	preprocessing.m:					Conversion of spike times to burstlet trains
	GCAnalysis_script.m:				Script performing GC network analysis of MEA burstlet trains
	SynchAnalysis_script.m:				Script performing analysis of higher-order burstlet synchrony between electrodes


Injury+Recovery Scripts:
	preprocessing.m:					Conversion of spike times to burstlet trains
	GCAnalysis_script.m:				Script performing GC network analysis of MEA burstlet trains
	SynchAnalysis_script.m:				Script performing analysis of higher-order burstlet synchrony between electrodes
	
	
Codes:
	omp.m:							Implementation of OMP for discretized point process models of burstlet trains
	omp_cv.m:						Cross-validate for sparsity level of discretized point process model
	getDesMat.m:					Construct set of history covariates for GC analysis
	
	FDRcontrolBY.m:					Statistical inference of GC links with Benjamini-Hochbery-Yakutelli procedure for False Discovery Rate control 
	SynchTest_static.m:				Statistical inference of rth-order coordinated spiking
	
										
Codes were developed and tested using MATLAB R2017b.
