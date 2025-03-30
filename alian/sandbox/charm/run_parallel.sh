#!/bin/bash
# filepath: /Users/ploskon/devel/alian/alian/sandbox/charm/run_parallel.sh

n_jobs=$2
if [ -z "$n_jobs" ]; then
    n_jobs=5
fi

# Determine number of CPUs based on OS
if [[ "$(uname)" == "Darwin" ]]; then
    ncpu=$(sysctl -n hw.ncpu)
else
    ncpu=$(nproc)
fi

# Limit concurrent jobs to ncpus/2
concurrent_jobs=$(( ncpu - 1 ))
if [ $concurrent_jobs -lt 1 ]; then
    concurrent_jobs=1
fi
echo "Total requested job count: $n_jobs"
echo "Concurrent jobs limited to: $concurrent_jobs"

function make_job_string() {
    nev=$1
    seed=$2
    outfname=$3
    # replace .root with _seed.root
    outfname=${outfname/.root/_$seed.root}
    # job_string="./pythia_charm.py --nev $nev --py-seed $seed -o $outfname --py-ecm 5020 > log_$seed.txt 2>&1"
    # job_string="./pythia_charm_process.py --nev $nev --py-seed $seed -o $outfname --py-ecm 5020 > log_$seed.txt 2>&1"
    job_string="./pythia_charm_process.py --nev $nev --py-seed $seed -o $outfname --py-ecm 5020 --D0required > log_$seed.txt 2>&1"
    echo $job_string
}
export -f make_job_string

nev=$1
if [ -z "$nev" ]; then
    nev=10000
fi
jobs=()
for i in $(seq 1 $n_jobs); do
    seed=$((10000 + i))
    # outfname="pythia_charm.root"
    outfname="pythia_charm_process.root"
    job_string=$(make_job_string $nev $seed $outfname)
    jobs+=("$job_string")
done

njobs=${#jobs[@]}
echo "Running ${njobs} jobs in parallel"
for j in "${jobs[@]}"; do
    echo $j
done

parallel --progress --jobs $concurrent_jobs ::: "${jobs[@]}"
rm -rfv *.log
hadd -f pythia_charm_process.root pythia_charm_process_*.root
