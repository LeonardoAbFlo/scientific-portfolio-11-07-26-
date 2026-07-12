# RNA hairpin molecular dynamics

A staged GROMACS workflow for preparing RNA hairpin systems, running energy
minimization/equilibration/production simulations, and quantifying structural and
thermodynamic behavior across replicates.

## Stages

1. Download the CHARMM force field.
2. Build the solvated and ionized RNA system.
3. Prepare molecular-dynamics input files.
4. Run minimization, NVT, and NPT equilibration.
5. Run production replicates with checkpoint continuation.
6. Extract thermodynamic observables.
7. Prepare centered, fitted trajectories.
8. Calculate RMSD, RMSF, radius of gyration, and hydrogen-bond metrics.
9. Calculate the hairpin index.
10. Summarize hairpin-specific measurements across replicates.

Scripts are numbered in execution order. Set `GMX` to the GROMACS executable
(`gmx_mpi` by default) and `GMXLIB` to the force-field parent directory.

```bash
FORCE_FIELD_DIR="$PWD/force-fields" bash 01_download_force_field.sh
export GMXLIB="$PWD/force-fields"
```

Simulation parameters and input structures remain study-specific and should be
reviewed before use on a new system.

