gmx_mpi grompp -f em.mdp  -c solv_ions.gro -p topol.top -n index.ndx -o em.tpr
CUDA_VISIBLE_DEVICES=0 gmx_mpi mdrun -deffnm em -nb gpu -pme cpu -ntomp 12 -pin on

gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
CUDA_VISIBLE_DEVICES=0 gmx_mpi mdrun -deffnm nvt -nb gpu -pme gpu -ntomp 12 -pin on

gmx_mpi grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
CUDA_VISIBLE_DEVICES=0 gmx_mpi mdrun -deffnm npt -nb gpu -pme gpu -ntomp 12 -pin on
