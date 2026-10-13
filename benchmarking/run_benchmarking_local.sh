#!/bin/bash
# Runs one benchmarking workflow on a local computer, with a hard memory limit.
# Usage, from benchmarking/ inside tmux: bash run_benchmarking_local.sh semi
# Workflows: gen_refs, semi, pereira, synthetic, compare_all
# Options after the workflow name go to snakemake, for example:
#   bash run_benchmarking_local.sh semi --rerun-triggers mtime
#
# Set these variables before the command to fit a different computer.
# Memory values are in MB, except CGROUP_MAX.
#   CORES          CPU cores for snakemake and for each HoMi run (default 4)
#   HOMI_MEM_MB    memory for each HoMi run (default 24000)
#   KRAKEN_MEM_MB  memory for each kraken job inside HoMi (default 18000)
#   SKIP_HUMANN    true or false (default true). false also runs MetaPhlAn and HUMAnN.
#   HUMANN_MEM_MB  memory for each HUMAnN job inside HoMi (default 22000).
#                  Must not be more than HOMI_MEM_MB.
#   SMK_MEM_MB     memory for all jobs that run at the same time (default 26000).
#                  Must be at least HOMI_MEM_MB.
#   CGROUP_MAX     hard limit (default 28G). Above it, the kernel stops pipeline
#                  processes, not the OS or sshd. Keep it a few GB below total RAM.
# Example: CORES=8 HOMI_MEM_MB=40000 SMK_MEM_MB=48000 CGROUP_MAX=56G bash run_benchmarking_local.sh semi
set -euo pipefail
# --conda-frontend conda: Snakemake 7 cannot create environments with mamba 2,
# because it makes the environment directory before it calls mamba.

RUN=${1:?"usage: run_benchmarking_local.sh <gen_refs|semi|pereira|synthetic|compare_all>"}
CORES=${CORES:-4}
HOMI_MEM_MB=${HOMI_MEM_MB:-24000}
KRAKEN_MEM_MB=${KRAKEN_MEM_MB:-18000}
SKIP_HUMANN=${SKIP_HUMANN:-true}
HUMANN_MEM_MB=${HUMANN_MEM_MB:-22000}
SMK_MEM_MB=${SMK_MEM_MB:-26000}
CGROUP_MAX=${CGROUP_MAX:-28G}

if (( SMK_MEM_MB < HOMI_MEM_MB )); then
    echo "ERROR: SMK_MEM_MB ($SMK_MEM_MB) must be at least HOMI_MEM_MB ($HOMI_MEM_MB)." >&2
    exit 1
fi
if (( HUMANN_MEM_MB > HOMI_MEM_MB )); then
    echo "ERROR: HUMANN_MEM_MB ($HUMANN_MEM_MB) must not be more than HOMI_MEM_MB ($HOMI_MEM_MB)." >&2
    exit 1
fi

SMK_CONFIG=(run="$RUN" executor=local
            homi_cores="$CORES" homi_mem_mb="$HOMI_MEM_MB" kraken_mem_mb="$KRAKEN_MEM_MB"
            skip_humann="$SKIP_HUMANN" humann_mem_mb="$HUMANN_MEM_MB")

snakemake -s benchmarking.smk --config "${SMK_CONFIG[@]}" --unlock

systemd-run --user --scope --unit="homi-bench-$RUN-$$" \
    -p MemoryMax="$CGROUP_MAX" -p MemorySwapMax=0 \
    snakemake -s benchmarking.smk \
        --config "${SMK_CONFIG[@]}" \
        --cores "$CORES" \
        --resources mem_mb="$SMK_MEM_MB" homi_runs=1 \
        --use-conda --conda-frontend conda --rerun-incomplete \
        "${@:2}" \
    2>&1 | tee "local_$RUN.log"
