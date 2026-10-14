#!/bin/bash
# Run from this directory; never submits jobs or deletes remote/local results.
set -euo pipefail
remote="taravser@roihu-cpu.csc.fi"
remote_dir="/scratch/project_2004996/geosse-benchmark/publication-constant-100"
ssh_args=(-o BatchMode=yes -o ConnectTimeout=15 -i "$HOME/.ssh/roihu" -o "CertificateFile=$HOME/.ssh/roihu-cert.pub")
mkdir -p results
rsync -az --exclude='*.tmp' -e "ssh -o BatchMode=yes -o ConnectTimeout=15 -i $HOME/.ssh/roihu -o CertificateFile=$HOME/.ssh/roihu-cert.pub" "$remote:$remote_dir/results/" results/
ssh "${ssh_args[@]}" "$remote" 'sacct -j 2300270,2300275 --allocations -n -P --format=JobID,State,ExitCode,Elapsed' > accounting-inference.txt
ssh "${ssh_args[@]}" "$remote" 'sacct -j 2300286 --allocations -n -P --format=JobID,State,ExitCode,Elapsed' > accounting-summary.txt
Rscript status.R
Rscript inventory.R
