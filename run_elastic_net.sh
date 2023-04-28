set -e

n_jobs=16
query="cpu-b.q@cpu-b*"

enable="enable"
disable="disable"
src="./pmopt/pmopt"
outputs="a"

# data-dependent options
input=$1; maxpat=$2; loss=$3

# parameterized common options
# l1_step * n_l1 = l1_min => n_l1 = l1_min / l1_step = 2/l1_step
# l1_step_list="0.05 0.1 0.2 0.4" # => n_l1 = 40 20 10 5
l1_step_list="0.1 0.2 0.4" # => n_l1 = 40 20 10 5

# common options
maxiter="1000000"
maxtrial="1000"
verbose="1"
l1_min="1e-2"
maxfolds="10"
seed=0

cnt=1


run_one () {
    input=$1; maxpat=$2; loss=$3
    gap_select=$4; multi_screen=$5; multi_update=$6; multi=$7

    output="${outputs}/${input}-maxpat${maxpat}-${loss}"
    output="${output}-multiupdate-cv/l1step${l1_step}"
    output="${output}-multi${multi}"
    output="${output}-mupdate${multi_update}"
    output="${output}-mscreen${multi_screen}"
    output="${output}-gapselect${gap_select}"

    cmd="${src} mine predictive inputs/${input}.csv
        --loss ${loss}
        --${multi_update}-multi-update
        --${multi_screen}-multi-screen
        --${gap_select}-gap-select
        --multi-dynamic-screen ${multi}
        --l1-step ${l1_step}
        --l1-min ${l1_min}
        --maxiter ${maxiter}
        --maxpat ${maxpat}
        --maxtrial ${maxtrial}
        --verbose ${verbose}
        --output-dir ${output}
        --seed ${seed}
        --dry-run
    "

script="\
#!/bin/csh\n\
#$ -cwd\n\
#$ -V -S /bin/bash\n\
#$ -q ${query}\n\
#$ -N SPP${cnt}\n\
#$ -pe smp ${n_jobs}\n\
#$ -e ${output}.err\n\
#$ -o ${output}.out\n\
${cmd}
"
    mkdir -p ${output}

    echo -e ${script} > ${output}.csh
    qsub ${output}.csh
    cnt=`expr ${cnt} + 1`
    echo ${cnt}:${output} >> output.log
}

run_all () {
    echo

    for l1_step in ${l1_step_list}
    do
        # input loss maxpat gap-select multi-screen multi-update n-multi 
        run_one $1 $2 $3 disable disable disable 0 # existing
        run_one $1 $2 $3 enable enable disable 0 # proposed (no multi update)
        run_one $1 $2 $3 enable enable enable 1 # M=1
        run_one $1 $2 $3 enable enable enable 2 # M=2
        run_one $1 $2 $3 enable enable enable 4 # M=4
    done
}

run_all "rhodopsin" "50" "squared"
run_all "dna" "3" "squared"
run_all "a1a" "5" "squaredhinge"
run_all "w1a" "3" "squaredhinge"
run_all "splice" "3" "squaredhinge"
run_all "a9a" "5" "squaredhinge"