## Simulations 

This code is largely from Jusuf et al., NSMB (2026). The only additions are (1) incorporation of increased cohesin loading probability at discrete sets of monomers and (2) some refactoring to make parameter sweeps easier to track and record.

The main 3D simulation files is volume_density_sweep.py, which can run either or both of chromatin compaction and targeted loading simulations. This runs the 1D cohesin position simulation DSB_SMCTranslocator_v2.pyx, written by Hugo Brandao for Yang et al., Nat Comms (2022).

The auxiliary files are:

process_sim_bulk.py makes contact maps and gets loop strengths. Can specify balancing, loop params, cutoff radii, etc.
1D_sim_trajectories.ipynb runs the cohesin accumulation simulation shown in Fig. S12A of the manuscript.
sim_param_check.py runs a "dummy" of the simulation that just shows where all the elements are (for verification purposes). This generates the schematic in Fig. S11A.
