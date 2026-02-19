#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 23 ]; then
  echo "Usage: $0 <outDir> <inSuffix> <setId> <epsblob> <epsij_scale> <surfself> <sigdeli> <sigdelf> <sigactdeli> <sigactdelf> <Temi> <Temf> <lmni> <lmnf> <surfDiffu> <sgactau> <actscalei> <actscalef> <pref> <start> <start_anneal> <simlength> <Tp>" >&2
  exit 1
fi

out_dir="$1"
in_suffix="$2"
set_id="$3"
epsblob="$4"
epsij_scale="$5"
surfself="$6"
sigdeli="$7"
sigdelf="$8"
sigactdeli="$9"
sigactdelf="${10}"
Temi="${11}"
Temf="${12}"
lmni="${13}"
lmnf="${14}"
surf_diffu="${15}"
sgactau="${16}"
actscalei="${17}"
actscalef="${18}"
pref="${19}"
start="${20}"
start_anneal="${21}"
simlength="${22}"
Tp="${23}"

calc() {
  echo "scale=6; $*" | bc
}

# Derived simulation scales used in the parameter file.
kspring="$(calc "$epsblob / 8")"
epshertz="$(calc "$epsblob * 4.0")"
surftauinv="$(calc "$surfself * 1.0")"
surfcross="$(calc "$surfself * 0.125")"
lmn_scalei="$(calc "$epsij_scale * $lmni")"
lmn_scalef="$(calc "$epsij_scale * $lmnf")"
seed="$(calc "$set_id * 1000 + $actscalei * 100 + $sigactdeli * 100 + 5 * $sgactau")"

# Select ellipsoid axes from geometry mode.
case "$pref" in
  ell)
    ea=9
    eb=6
    ec=2.5
    ;;
  sph)
    ea=5.5
    eb=5.5
    ec=5.5
    ;;
  *)
    echo "Error: pref must be 'ell' or 'sph' (got '$pref')." >&2
    exit 1
    ;;
esac

# Build canonical Set_<set>_in<inSuffix>_<pref>_Tp<Tp> folders.
main_dir="Set_${set_id}_in${in_suffix}_${pref}_Tp${Tp}"
starts_dir="${in_suffix}_${pref}_Tp${Tp}_starts"
surfcyc=1000
if [ "$Tp" = "0" ]; then
  main_dir="Set_${set_id}_in${in_suffix}_${pref}"
  starts_dir="${in_suffix}_${pref}_starts"
  surfcyc=500
fi

mkdir -p "$main_dir/$out_dir" "ParaFiles"
par_file="ParaFiles/Para_test_46c_whole_${set_id}.par"

cat > "$par_file" <<PAR
Nchrom            46
nPatchTot         1058
Temi              ${Temi}
Temf              ${Temf}
dt                0.001
gamma             10
rotgamma          100
chromAct          0
actscalei         ${actscalei}
actscalef         ${actscalef}
rho               0.1
pressure          0
init_read         ${start}
startGeom         0
seed              ${seed}
Label             Conf
TSfact            2
eqRun             ${start_anneal}
totRun            ${simlength}
thermo            50
dump              1000
ellai             ${ea}
ellbi             ${eb}
ellci             ${ec}
ellaf             9
ellcf             2.5
dellwin           5000
eps_blob          ${epsblob}
springscale       ${kspring}
espij_scale       ${epsij_scale}
eps_hertz         ${epshertz}
nSurfBasicPoints  2500
surfNeighfact     1.75
nSurfCycPoints    ${surfcyc}
refineRounds      5
delta_Ar_th       0.1
surfPhi           0.35
surfSelf          ${surfself}
surfCrossi        ${surfcross}
surfCrossf        ${surfcross}
surfDiffu         ${surf_diffu}
lamin_scalei      ${lmn_scalei}
lamin_scalef      ${lmn_scalef}
taup              100
sigdeli           ${sigdeli}
sigdelf           ${sigdelf}
nGrid             1000
tau_theta         ${sgactau}
sigactdeli        ${sigactdeli}
sigactdelf        ${sigactdelf}
surftauinv        ${surftauinv}
surfrestore       1.0
hic_cut           0.5
hic_decay         0.25
PAR

cp "$par_file" "$main_dir/$out_dir/"

surf_dir="prepSurf_perturb"
cp "$main_dir/$surf_dir/Conf-finalfixedSurfFile_ea${ea}_eb${eb}_ec${ec}_5.dat" "$main_dir/$out_dir/"
cp "$main_dir/$surf_dir/Conf-SurfCyc_ea${ea}_eb${eb}_ec${ec}_5.dat" "$main_dir/$out_dir/"

# Run in foreground so Slurm task accounting tracks completion correctly.
./mdChrom_batch_gr.exe \
  "$par_file" \
  "$starts_dir/testWhole_46c-SimStart_${in_suffix}_${set_id}.dat" \
  "$starts_dir/testWhole_46c-IntMat_${in_suffix}_${set_id}.dat" \
  "$main_dir/$out_dir/"
