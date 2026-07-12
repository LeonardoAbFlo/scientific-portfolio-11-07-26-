#!/usr/bin/env bash
set -Eeuo pipefail

: "${GMXLIB:?Set GMXLIB to the directory containing charmm36-mar2019.ff}"
GMX="${GMX:-gmx_mpi}"

printf "2\n3\n" | "$GMX" pdb2gmx -f RNA_charmm_fix.pdb -o processed.gro -p topol.top \
  -ff charmm36-mar2019 -water tip3p -ignh -ter

"$GMX" editconf -f processed.gro -o boxed.gro -bt dodecahedron -d 1.2

"$GMX" solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top

"$GMX" grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
"$GMX" genion -s ions.tpr -o solv_ions.gro -p topol.top \
  -pname K -nname CL -neutral -conc 0.15

"$GMX" make_ndx -f solv_ions.gro -o index.ndx
