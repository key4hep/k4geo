set -e

# helper function to selectively print and run commands without a subshell
function run() {
    set -x
    "$@"
    # save exit code
    { rec=$?; } 2> /dev/null
    { set +x;   } 2> /dev/null
    # restore exit code
    (exit $rec)
}

export run

newpath=${1:-$GITHUB_WORKSPACE/new_hist.root}
refpath=${2:-$GITHUB_WORKSPACE/ref_hist.root}
outdir=${3:-$GITHUB_WORKSPACE}
mkdir -p $outdir
mkdir -p $outdir/html

# File to accumulate the histcmp results
histcmp_results=$outdir/histcmp_results.csv
echo -n "" > $histcmp_results

function run_histcmp() {
    a=$1
    b=$2
    title=$3
    html_path=$4
    plots_path=$5
    shift 5

    echo "::group::Comparing $a vs. $b"

    ec=0
    if [ ! -f "$a" ]; then
        echo "::error::histcmp failed: File $a does not exist"
        ec=1
    fi

    if [ ! -f "$b" ]; then
        echo "::error::histcmp failed: File $b does not exist"
        ec=1
    fi

    mkdir -p "$(dirname "$outdir/html/$html_path")"
    mkdir -p "$outdir/html/$plots_path"

    run histcmp $a $b \
        --label-reference=k4geo-master \
        --label-monitored=local-branch \
        --title="$title" \
        -o $outdir/html/$html_path \
        -p $outdir/$plots_path \
        "$@"

    this_ec=$?
    ec=$(($ec | $this_ec))

    if [ $this_ec -ne 0 ]; then
        echo "::error::histcmp failed: ec=$this_ec"
    fi

    echo "\"${title}\",html/${html_path},${this_ec}" >> $histcmp_results

    echo "::endgroup::"
}

run_histcmp \
    $newpath \
    $refpath \
    "Test comparison" \
    detector_compaisons.html \
    detector_comparison_plots \
    --config CI/config/test_config.yml