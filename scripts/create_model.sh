#!/usr/bin/env bash
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16

export PATH=~/code/SpotSeq/src:$PATH

INPUT=
NUMSTATES=100
NITER=
RUNNAME=

function usage()
{
    cat <<EOF
usage: $0  -i <path to fasta files> -s <initial number of states> -n <number of iterations> -r <runname>
EOF
    exit 1;
}

while getopts i:s:n:r: opt
do
    case ${opt} in
        r) RUNNAME=${OPTARG};;
        i) INPUT=${OPTARG};;
        s) NUMSTATES=${OPTARG};;
        n) NUMITER=${OPTARG};;
        *) usage;;
    esac
done

if [ "${INPUT}" = "" ]; then usage; fi

#
#   Sanity check 
#

programs=(spotseq_model spotseq_plot dot spotseq_score) 

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

OUTMODEL=$INPUT$RUNNAME".model"
OUTDOT=$INPUT$RUNNAME".dot"
OUTPDF=$INPUT$RUNNAME".pdf"
OUTSCORES=$INPUT$RUNNAME".scores.csv"


echo "Running: spotseq_model -in $INPUT --states 100 -out $OUTMODEL"
spotseq_model -in $INPUT --states $NUMSTATES -out $OUTMODEL -niter $NUMITER

echo "Running: spotseq_plot -in $OUTMODEL  -out $OUTDOT"
spotseq_plot -in $OUTMODEL  -out $OUTDOT

echo "Running: dot -Tpdf $OUTDOT -o $OUTPDF"
dot -Tpdf $OUTDOT -o $OUTPDF

echo "Running: spotseq_score  -m $OUTMODEL  -i $INPUT  -o $OUTSCORES"
spotseq_score  -m $OUTMODEL  -i $INPUT  -o $OUTSCORES 










