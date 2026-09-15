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


This directory is used to compute the POD approximation for different choices on M and r for the brusselator equation with beta=3, alpha=1 and nu=0.002.

To get the POD approximation, run 'todo' and define the number of the stepsizes in the POD basis (M) and the number of modes (r). Results are saved in 'results_POD_'M'_r' files. POD basis are saved in '00_compute_pod_basis' directory. A tensor structure is used to compute the POD basis at 2049 time instants.

To get the results in Table 1, 2 and 3, run 'tabla_de_errores' file and select the corresponding columns in the error matrices. 