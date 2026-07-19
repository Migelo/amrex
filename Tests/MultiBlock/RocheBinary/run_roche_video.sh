#!/usr/bin/env bash
# run_roche_video.sh OUTNAME [roche.key=value ...]
#
# End-to-end RocheBinary pipeline: run the sim with the given roche.* params
# in an isolated workdir (video_runs/OUTNAME/), render one rho frame per
# plotfile step, assemble an h264 MP4, and publish it to the static site.
#
# Extra args are passed verbatim to the executable; roche.plot_int=10 is
# added unless already supplied. Star geometry for the renderer's masks is
# derived from roche.m1/m2/sep/a (defaults 1 1 4 1):
#   x1 = -sep*m2/(m1+m2),  x2 = +sep*m1/(m1+m2)
# Env knobs: FPS (default 12), CRF (default 15), DRIFT=1 (render relative
# (rho-rho0)/rho0 instead of rho -- for subtle-transfer runs).
#
# Examples:
#   ./run_roche_video.sh baseline roche.perturb=0.05 roche.stop_time=4.0
#   ./run_roche_video.sh m2_05 roche.m2=0.5 roche.stop_time=2 roche.perturb=0.05
set -euo pipefail

TEST_DIR=/home/cernetic/amrex/Tests/MultiBlock/RocheBinary
EXE="$TEST_DIR/main2d.gnu.MPI.ex"
SITE=/home/cernetic/stardisk-site
RENDER="$TEST_DIR/render_frames.py"
NIX_PY='with import (builtins.fetchTarball { url = "https://github.com/NixOS/nixpkgs/archive/nixos-unstable.tar.gz"; }) {}; mkShell { packages = [ (python3.withPackages (ps: [ ps.numpy ps.matplotlib ])) ]; }'

if [ $# -lt 1 ]; then
    sed -n '2,20p' "$0"
    exit 1
fi
OUTNAME=$1; shift

# Pass-through exe args; add plot_int=10 unless the user supplied one.
EXE_ARGS=("$@")
have_plot_int=0
M1=1; M2=1; SEP=4; A=1
for arg in "$@"; do
    case "$arg" in
        roche.plot_int=*) have_plot_int=1 ;;
        roche.m1=*)  M1=${arg#roche.m1=} ;;
        roche.m2=*)  M2=${arg#roche.m2=} ;;
        roche.sep=*) SEP=${arg#roche.sep=} ;;
        roche.a=*)   A=${arg#roche.a=} ;;
    esac
done
[ "$have_plot_int" -eq 1 ] || EXE_ARGS+=(roche.plot_int=10)

X1=$(awk -v s="$SEP" -v m1="$M1" -v m2="$M2" 'BEGIN { printf "%.6f", -s*m2/(m1+m2) }')
X2=$(awk -v s="$SEP" -v m1="$M1" -v m2="$M2" 'BEGIN { printf "%.6f",  s*m1/(m1+m2) }')

WORKDIR="$TEST_DIR/video_runs/$OUTNAME"
FRAMES="$WORKDIR/frames"
mkdir -p "$WORKDIR"
rm -rf "$WORKDIR/RocheBinary" "$FRAMES"   # drop stale steps from prior runs

echo "== [$OUTNAME] sim: ${EXE_ARGS[*]}"
echo "== [$OUTNAME] workdir: $WORKDIR; stars at x1=$X1 x2=$X2 a=$A"
(cd "$WORKDIR" && direnv exec . mpirun -n 2 "$EXE" "${EXE_ARGS[@]}")

echo "== [$OUTNAME] rendering frames"
nix-shell -E "$NIX_PY" --run "python3 $RENDER --base $WORKDIR --outdir $FRAMES --star-x1=$X1 --star-x2=$X2 --star-a=$A ${DRIFT:+--drift}"

echo "== [$OUTNAME] encoding ${FPS:-12} fps, crf ${CRF:-15}"
ffmpeg -y -framerate "${FPS:-12}" -i "$FRAMES/frame_%04d.png" \
    -c:v libx264 -crf "${CRF:-15}" -pix_fmt yuv420p \
    -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -movflags +faststart \
    "$WORKDIR/$OUTNAME.mp4"

cp "$WORKDIR/$OUTNAME.mp4" "$SITE/$OUTNAME.mp4"
echo "== [$OUTNAME] published: $SITE/$OUTNAME.mp4"
ffprobe -v error -show_entries stream=codec_name,pix_fmt \
    -show_entries format=duration,size -of default=nw=1 "$SITE/$OUTNAME.mp4"
