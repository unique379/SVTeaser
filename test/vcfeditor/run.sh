# Example run command for SURVIVOR to generate test data
# Please have SURVIVOR installed and in the $PATH var prior to running this script.


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
OUT_DIR=$SCRIPT_DIR/survivor_output

rm -rf $OUT_DIR
mkdir $OUT_DIR

SURVIVOR simSV $SCRIPT_DIR/ref.fa $SCRIPT_DIR/parameter_file 0.0 0 $OUT_DIR/simulated
python3 $SCRIPT_DIR/../../svteaser/vcfeditor.py -r $SCRIPT_DIR/ref.fa -i $OUT_DIR/simulated.insertions.fa -v $OUT_DIR/simulated.vcf -o $SCRIPT_DIR/out.vcf
