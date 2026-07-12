#!/usr/bin/env bash
set -euo pipefail

NREPS=5
NTOMP=12
TEMP=310
GPUFLAGS="-nb gpu -pme gpu -pin on"
CPUS_GPU0="0-15"
CPUS_GPU1="16-31"

GMX=${GMX:-gmx_mpi}
BASE_DIR="$(pwd)"

for f in npt.gro topol.top index.ndx md.mdp; do
  [[ -f "$f" ]] || { echo "Falta $f en $BASE_DIR"; exit 1; }
done

echo "[OK] Archivos base encontrados en: $BASE_DIR"
echo "GMX=$GMX | NREPS=$NREPS | NTOMP=$NTOMP | TEMP=$TEMP"

fix_mdp_comments () {
  local mdp="$1"
  sed -i 's/^[[:space:]]*#[[:space:]]*/; /' "$mdp"
  sed -i '/^[[:space:]]*title[[:space:]]*=/d' "$mdp"
}

ensure_replica_vel_init () {
  local mdp="$1"
  local seed="$2"
  fix_mdp_comments "$mdp"

  if grep -q '^[[:space:]]*continuation[[:space:]]*=' "$mdp"; then
    sed -i 's/^[[:space:]]*continuation[[:space:]]*=.*/continuation = no/' "$mdp"
  else
    echo "continuation = no" >> "$mdp"
  fi

  if grep -q '^[[:space:]]*gen_vel[[:space:]]*=' "$mdp"; then
    sed -i 's/^[[:space:]]*gen_vel[[:space:]]*=.*/gen_vel      = yes/' "$mdp"
  else
    echo "gen_vel      = yes" >> "$mdp"
  fi

  if grep -q '^[[:space:]]*gen_temp[[:space:]]*=' "$mdp"; then
    sed -i "s/^[[:space:]]*gen_temp[[:space:]]*=.*/gen_temp     = $TEMP/" "$mdp"
  else
    echo "gen_temp     = $TEMP" >> "$mdp"
  fi

  if grep -q '^[[:space:]]*gen_seed[[:space:]]*=' "$mdp"; then
    sed -i "s/^[[:space:]]*gen_seed[[:space:]]*=.*/gen_seed     = $seed/" "$mdp"
  else
    echo "gen_seed     = $seed" >> "$mdp"
  fi
}

run_rep () {
  local rep="$1"
  local gpu="$2"
  local cpus="$3"
  (
    cd "$rep"
    echo "[RUN] $rep | GPU=$gpu | CPUs=$cpus | NTOMP=$NTOMP"
    taskset -c "$cpus" env CUDA_VISIBLE_DEVICES="$gpu" \
      $GMX mdrun -deffnm md $GPUFLAGS -ntomp "$NTOMP"
  )
}

echo "[1/3] Preparando réplicas y generando md.tpr..."
for i in $(seq 1 "$NREPS"); do
  rep="rep$i"
  seed=$((1000+i))

  mkdir -p "$rep"
  cp -f npt.gro topol.top index.ndx "$rep/"
  cp -f md.mdp "$rep/md.mdp"

  ensure_replica_vel_init "$rep/md.mdp" "$seed"

  (
    cd "$rep"
    $GMX grompp -f md.mdp -c npt.gro -p topol.top -n index.ndx -o md.tpr -maxwarn 1
  )
done

echo "[OK] md.tpr generado en rep1..rep$NREPS"
echo

echo "[2/3] Corriendo MD (2 réplicas en paralelo con 2 GPUs)..."

i=1
while [[ $i -le $NREPS ]]; do
  j=$((i+1))

  if [[ $j -le $NREPS ]]; then
    run_rep "rep$i" 0 "$CPUS_GPU0" & pid1=$!
    run_rep "rep$j" 1 "$CPUS_GPU1" & pid2=$!
    wait "$pid1" "$pid2"
  else
    run_rep "rep$i" 0 "$CPUS_GPU0"
  fi

  i=$((i+2))
done

echo
echo "[OK] Producción terminada."
echo

cat <<'EOF'
Monitoreo:
  tail -f rep1/md.log
  nvidia-smi dmon -s u

Reanudar si se corta (en una réplica):
  (cd rep1 && gmx_mpi mdrun -deffnm md -cpi md.cpt -nb gpu -pme gpu -ntomp 12 -pin on)

EOF
