# Example run command for SURVIVOR to generate test data
# Please have SURVIVOR installed and in the $PATH var prior to running this script.


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

SURVIVOR simSV $SCRIPT_DIR/ref.fa $SCRIPT_DIR/parameter_file 0.0 0 $SCRIPT_DIR/simulated
python $SCRIPT_DIR/../../sv_simulator/vcfeditor.py -r $SCRIPT_DIR/ref.fa -i $SCRIPT_DIR/simulated.insertions.fa -v $SCRIPT_DIR/simulated.vcf -o $SCRIPT_DIR/out.vcf
