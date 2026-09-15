Codes for numerical experiments in

Bosco García-Archilla, Alicia García-Mascaraque and Julia Novo
	Using BDF schemes in the temporal integration of POD-ROM methods
	Journal of Scientific Computing (to appear).
 Please check volume, year and pages with the journal for proper citation.
 Preprint: arXiv:2506.14543 [math.NA],   	
           https://doi.org/10.48550/arXiv.2506.14543

Feel free to use data and codes in this repository, but please cite the above-mentioned paper when doing so.

Codes written by Alicia García-Mascaraque. Codes come with no guarantee or warranty of any kind.

Run with Matlab R2025b or earlier. 


This directory is used to compute the POD basis using finite differences for different stepsizes M for the brusselator equation with beta=3, alpha=1 and nu=0.002.

Initial conditions of the problem are known from the data file: 'initial_conditions_nu0p002_alpha1_beta3.mat'. 

To get the snapshots at 2049 time instants, run  the code 'step1_compute_snapshots_diff. The results are saved in 'snapshots_alpha1_beta3_nu0p002'. 
Then, to get the finite differences for different M values (time step sizes), run 'step2_compute_snapshots_diff'. Results are saved in 'snapshots_diff_M' files.
To get the POD basis for the corresponding M, run 'step3_base_pod_diff'. Results are saved in 'base_pod_diff_M' files.

To make figure 1 in the article above, run 'make_figure1'.