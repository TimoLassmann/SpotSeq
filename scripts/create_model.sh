#!/usr/bin/env bash
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16

export PATH=~/code/SpotSeq/src:$PATH

INPUT=
NUMSTATES=100

function usage()
{
    cat <<EOF
usage: $0  -i <path to fasta files> -n <initial number of states
EOF
    exit 1;
}

while getopts i:n: opt
do
    case ${opt} in
        i) INPUT=${OPTARG};;
        n) NUMSTATES=${OPTARG};;
        *) usage;;
    esac
done

if [ "${INPUT}" = "" ]; then usage; fi

#
#   Sanity check 
#

programs=(spotseq_model spotseq_plot dot) 

printf "Running Sanity checks:\n";

for item in ${programs[*]}
do   
    if which $item >/dev/null; then
        printf "%15s found...\n"  $item;
    else
        printf "\nERROR: %s not found!\n\n" $item;
        exit 1;
    fi
done

echo "All dependencies found."

OUTMODEL=$INPUT".model"
OUTDOT=$INPUT".dot"
OUTPDF=$INPUT".pdf"

echo "Running: spotseq_model -in $INPUT --states 100 -out $OUTMODEL"
spotseq_model -in $INPUT --states $NUMSTATES -out $OUTMODEL

echo "Running: spotseq_plot -in $OUTMODEL  -out $OUTDOT"
spotseq_plot -in $OUTMODEL  -out $OUTDOT

echo "Running: dot -Tpdf $OUTDOT -o $OUTPDF"
dot -Tpdf $OUTDOT -o $OUTPDF




















