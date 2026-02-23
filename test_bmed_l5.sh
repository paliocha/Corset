#!/bin/bash
#SBATCH -J corset-test-BMED-l5
#SBATCH -p orion
#SBATCH -c 32
#SBATCH --mem 131072M
#SBATCH -t 48:00:00
#SBATCH -o corset-test-BMED-l5_%j.log
#SBATCH -e corset-test-BMED-l5_%j.err

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────
CORSET=/net/fs-2/scale/OrionStore/Home/martpali/AnnualPerennial/Corset/corset
WORKBASE=/net/fs-2/scale/OrionStore/Projects/FjellheimLab/martpali/AnnualPerennial/nf-denovoslim/BMED/work
OUTDIR=/net/fs-2/scale/OrionStore/Home/martpali/AnnualPerennial/Corset/test_output
CONDA_LIB=/mnt/users/martpali/micromamba/envs/samtools/lib

export OMP_NUM_THREADS=32
export LD_LIBRARY_PATH="${CONDA_LIB}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

mkdir -p "$OUTDIR"
cd "$TMPDIR"

# ── Stage eq_classes symlinks ────────────────────────────────────────
declare -A QUANT_DIRS=(
    [BMED02_T1_L]=81/2b9ef8c215dc17724caa4a9d6ba067
    [BMED03_T1_L]=e6/8d550efdac356dbd75079594bbca81
    [BMED04_T1_L]=9e/179a99775d0ae7aeb37f30bae7369e
    [BMED05_T1_L]=f0/12b88b4eb67d6eb1455afbd930d874
    [BMED02_T1_R]=2f/123a00f8eda877b8d7972a60f059a0
    [BMED03_T1_R]=5b/9cb139afead91cbf0c945af0152185
    [BMED04_T1_R]=a8/863af13e046d2bd568aaa281fc3ef0
    [BMED05_T1_R]=b8/f125fc7be5efd798f64621726ab5be
    [BMED13_T2_L]=02/0f6eec0b0cd42dd910db69d33a548d
    [BMED15_T2_L]=87/4d6bf545133d6fcaaf83cf7e615b81
    [BMED17_T2_L]=7c/0fc5b550deb2a54a78b7c0e8714b19
    [BMED19_T2_L]=8f/923e68bc09375a1f1b4b763ecc7537
    [BMED13_T2_R]=f5/4185d3dee0a910795c5f8131c79b7a
    [BMED15_T2_R]=2a/335d786eb265a91473cf6bbfbbf89e
    [BMED17_T2_R]=9f/10fcaf0086f0bc4e82ec63ba104874
    [BMED19_T2_R]=e3/1bed1cabdeeb5426b91ed691f9113e
    [BMED35_T3_L]=6f/027efaae25ddb54accbd9308bc5389
    [BMED36_T3_L]=02/a493bc353f4ce29dfd38da8cf6e111
    [BMED37_T3_L]=c1/4356c45ea0ea52273ba46b5dd01065
    [BMED38_T3_L]=cd/6a715fd73831bed0106273cef295ad
    [BMED35_T3_R]=03/c4ab420f0592a1b63d632a6dd6f51c
    [BMED36_T3_R]=d1/5f6aa5344a992d51b9b390dcfc892f
    [BMED37_T3_R]=61/4fca4780dd513fd7d7cb6ea9cedde4
    [BMED38_T3_R]=80/ca9a9d8dd727cb9ebead9451344f4a
    [BMED56_T4_L]=8e/7061428a258529087cb91018c4efc4
    [BMED57_T4_L]=87/e7fe90183be896c6e18db92f867c45
    [BMED58_T4_L]=6e/0ef1ec23f0494f8d488c9c310827e1
    [BMED60_T4_L]=63/e02af6ea1db528e682c04356843d20
    [BMED56_T4_R]=e0/2247bcb11248b91a2e2a411505debf
    [BMED57_T4_R]=68/c5b0d10a4f0af2f3e5ffbe84afc1e9
    [BMED58_T4_R]=92/c7bf4a5c1aa55d07db2cbeedf21768
    [BMED60_T4_R]=a4/4b69c1a282ae213f2fb694d040f123
    [BMED75_T5_L]=fc/629027cb991787a31efbd60a1d5243
    [BMED77_T5_L]=bd/2284334a334954e9f1747d5b9f8a7d
    [BMED79_T5_L]=35/1e2d9b2624e66896f7307e5b0539a9
    [BMED80_T5_L]=fe/85a08a4ee1c978644172ec11fc88cf
    [BMED75_T5_R]=aa/22a3c907ea76275a1042cc5b0dbb70
    [BMED77_T5_R]=a2/da16bdf28cb9bb6124a7eed40567d9
    [BMED79_T5_R]=5b/ea3e84d371a2d842a389808ab2c61a
    [BMED80_T5_R]=7b/f2bb426c00637271cc35a80f85a928
)

for sample in "${!QUANT_DIRS[@]}"; do
    target="${WORKBASE}/${QUANT_DIRS[$sample]}/${sample}_quant"
    link="${TMPDIR}/${sample}_quant"
    ln -sfn "$target" "$link"
done

# ── Shared parameters ────────────────────────────────────────────────
SAMPLE_GROUPS="T1_L,T1_L,T1_L,T1_L,T1_R,T1_R,T1_R,T1_R,T2_L,T2_L,T2_L,T2_L,T2_R,T2_R,T2_R,T2_R,T3_L,T3_L,T3_L,T3_L,T3_R,T3_R,T3_R,T3_R,T4_L,T4_L,T4_L,T4_L,T4_R,T4_R,T4_R,T4_R,T5_L,T5_L,T5_L,T5_L,T5_R,T5_R,T5_R,T5_R"

SAMPLE_NAMES="BMED02_T1_L,BMED03_T1_L,BMED04_T1_L,BMED05_T1_L,BMED02_T1_R,BMED03_T1_R,BMED04_T1_R,BMED05_T1_R,BMED13_T2_L,BMED15_T2_L,BMED17_T2_L,BMED19_T2_L,BMED13_T2_R,BMED15_T2_R,BMED17_T2_R,BMED19_T2_R,BMED35_T3_L,BMED36_T3_L,BMED37_T3_L,BMED38_T3_L,BMED35_T3_R,BMED36_T3_R,BMED37_T3_R,BMED38_T3_R,BMED56_T4_L,BMED57_T4_L,BMED58_T4_L,BMED60_T4_L,BMED56_T4_R,BMED57_T4_R,BMED58_T4_R,BMED60_T4_R,BMED75_T5_L,BMED77_T5_L,BMED79_T5_L,BMED80_T5_L,BMED75_T5_R,BMED77_T5_R,BMED79_T5_R,BMED80_T5_R"

EQ_FILES=(
    BMED02_T1_L_quant/aux_info/eq_classes.txt
    BMED03_T1_L_quant/aux_info/eq_classes.txt
    BMED04_T1_L_quant/aux_info/eq_classes.txt
    BMED05_T1_L_quant/aux_info/eq_classes.txt
    BMED02_T1_R_quant/aux_info/eq_classes.txt
    BMED03_T1_R_quant/aux_info/eq_classes.txt
    BMED04_T1_R_quant/aux_info/eq_classes.txt
    BMED05_T1_R_quant/aux_info/eq_classes.txt
    BMED13_T2_L_quant/aux_info/eq_classes.txt
    BMED15_T2_L_quant/aux_info/eq_classes.txt
    BMED17_T2_L_quant/aux_info/eq_classes.txt
    BMED19_T2_L_quant/aux_info/eq_classes.txt
    BMED13_T2_R_quant/aux_info/eq_classes.txt
    BMED15_T2_R_quant/aux_info/eq_classes.txt
    BMED17_T2_R_quant/aux_info/eq_classes.txt
    BMED19_T2_R_quant/aux_info/eq_classes.txt
    BMED35_T3_L_quant/aux_info/eq_classes.txt
    BMED36_T3_L_quant/aux_info/eq_classes.txt
    BMED37_T3_L_quant/aux_info/eq_classes.txt
    BMED38_T3_L_quant/aux_info/eq_classes.txt
    BMED35_T3_R_quant/aux_info/eq_classes.txt
    BMED36_T3_R_quant/aux_info/eq_classes.txt
    BMED37_T3_R_quant/aux_info/eq_classes.txt
    BMED38_T3_R_quant/aux_info/eq_classes.txt
    BMED56_T4_L_quant/aux_info/eq_classes.txt
    BMED57_T4_L_quant/aux_info/eq_classes.txt
    BMED58_T4_L_quant/aux_info/eq_classes.txt
    BMED60_T4_L_quant/aux_info/eq_classes.txt
    BMED56_T4_R_quant/aux_info/eq_classes.txt
    BMED57_T4_R_quant/aux_info/eq_classes.txt
    BMED58_T4_R_quant/aux_info/eq_classes.txt
    BMED60_T4_R_quant/aux_info/eq_classes.txt
    BMED75_T5_L_quant/aux_info/eq_classes.txt
    BMED77_T5_L_quant/aux_info/eq_classes.txt
    BMED79_T5_L_quant/aux_info/eq_classes.txt
    BMED80_T5_L_quant/aux_info/eq_classes.txt
    BMED75_T5_R_quant/aux_info/eq_classes.txt
    BMED77_T5_R_quant/aux_info/eq_classes.txt
    BMED79_T5_R_quant/aux_info/eq_classes.txt
    BMED80_T5_R_quant/aux_info/eq_classes.txt
)

COMMON_ARGS="-i salmon_eq_classes -g $SAMPLE_GROUPS -n $SAMPLE_NAMES"

run_corset() {
    local label="$1"
    shift
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  Test: $label"
    echo "  $(date)"
    echo "════════════════════════════════════════════════════════════════"
    time "$CORSET" $COMMON_ARGS "$@" "${EQ_FILES[@]}"
    echo "  Finished: $(date)"
    echo ""
}

# ── Test 1: Leiden + kNN auto + -l 5 ──────────────────────────────────
run_corset "leiden + kNN auto + -l 5" \
    --algorithm leiden \
    --knn auto \
    -l 5 \
    -p test-leiden-knn-l5

# ── Test 2: Hierarchical + -l 5 ──────────────────────────────────────
run_corset "hierarchical + -l 5" \
    -l 5 \
    -p test-hier-l5

# ── Summary ──────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  All tests complete — $(date)"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "Output files:"
ls -lh test-*.txt 2>/dev/null || echo "  (none found)"
echo ""
echo "Line counts:"
wc -l test-*.txt 2>/dev/null || echo "  (none found)"

# ── Copy results from local scratch back to NFS ─────────────────────
echo ""
echo "Copying output files to $OUTDIR ..."
cp -v test-*.txt "$OUTDIR/" 2>/dev/null || echo "  (no output files to copy)"
echo "Done."
