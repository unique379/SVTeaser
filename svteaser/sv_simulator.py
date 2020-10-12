
# Main interface to generate structural variant simulations

import pysam
import pandas as pd
import argparse
import os
import subprocess
from random import randint
from vcfeditor import update_vcf
import vcf

def parse_args():
    """Parse command line arguments.

    Return:
        args : parsed argument object.

    """
    parser = argparse.ArgumentParser(
        description='Generate structural variant simulations. Takes in a reference file, regions to \
                     simulate SVs on and outputs the following per alteration/simulation of SV \
                     1) the unaltered region of reference in fasta \
                     2) The altered region of reference in fasta \
                     3) The SVs in a VCF file \
                     NOTE: Please ensure that SURVIVOR is available in the PATH variable.')

    parser.add_argument('--ref_file', type=str,
                        help='Path to reference genome file',
                        required=True)
    parser.add_argument('--sv_regions', type=str,
                        help='Comma separated file containing (chr, region_start, region_end). \
                        For every row, an SV of length randint(50, mx_variation) is generated with the region \
                        specified by (chr, start, end).\
                        chr, start, end \
                        chr22, 1000, 20000 \
                        chr22, 50000, 80000', required=False)
    parser.add_argument('--num_sv_regions', type=int,
                        help='Alternatively to the csv file defined by --sv_regions, user can also \
                              provide number of regions to generate SVs for. The programme will randomly \
                              choose locations within the genome to introduce the SVs. --sv_regions will be given priority \
                              if both options are provided.',
                        required=False)
    parser.add_argument('--len_sv_region', type=int,
                        help='The length of regions to create.',
                        required=False)
    parser.add_argument('--survivor_param_file', type=str,
                        help='This programme uses survivor tool for simulations. We take in the same param file \
                              as needed by survivor.',
                        required=False,
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "parameter_file"))
    parser.add_argument('--out_dir', type=str,
                        help='Path to output directory.', required=True) 
    args = parser.parse_args()
    return args


def generate_regions_from_file(regions_file):
    regions = pd.read_csv(regions_file, sep=", ")
    region_list = []
    for index, row in regions.iterrows():
        # 1. Fetch the region from reference
        chrm = row['chr']
        start = row['start']
        end = row['end']
        region_list.append((chrm, start, end))
    return region_list

def generate_random_regions(ref_file, region_length, num_regions):
    def generate_region(ref, length):
        chroms = []
        for chrm in ref.references:
            chroms.append(chrm)

        chrmidx = randint(1, len(chroms))
        chrm = chroms[chrmidx -1]
        start = randint(1, ref.get_reference_length(chrm) - length)
        end = start + length

        # Offset by 1 since we started from 1 base for randint
        return chrm, start - 1, end - 1

    ref = pysam.Fastafile(ref_file)

    region_list = []
    for idx in range(num_regions):
        chrm, start, end = generate_region(ref, region_length)
        region_list.append((chrm, start, end))
    return region_list

def add_fasta_entry(name, seq, fasta_fh):
    fasta_fh.write(">{}\n".format(name))
    fasta_fh.write("{}\n".format(seq))

def update_altered_fa(ref_seq, altered_ref_seq, padding):
    begin_seq = ref_seq[0:padding]
    end_seq = ref_seq[len(ref_seq)-padding: ]
    return begin_seq + altered_ref_seq + end_seq


def process_regions(ref_file, regions, out_dir, param_file):
    out_vcf_path = os.path.join(out_dir, "svteaser.sim.vcf")
    out_ref_fa_path = os.path.join(out_dir, "svteaser.ref.fa")
    out_altered_fa_path = os.path.join(out_dir, "svteaser.altered.fa")

    out_vcf_fh = None
    out_ref_fa_fh = open(out_ref_fa_path, "w+")
    out_altered_fa_fh = open(out_altered_fa_path, "w+")

    ref = pysam.FastaFile(ref_file)

    # Define padding in reference region where SVs are not to be inserted.
    padding = 800

    for i, (chrom, start, end) in enumerate(regions):
        # Track status.
        if (i + 1) % 50 == 0:
            print("Processed {}/{} regions...".format(i + 1, len(regions)))

        # Temporari dir.
        temp_dir = os.path.join(out_dir, "temp")
        os.mkdir(temp_dir)

        # Extract ref seequebce.
        name = "{}_{}_{}".format(chrom, start, end)
        ref_seq = ref.fetch(chrom, start, end)

        # Remove first 800bp and last 800bp, so that the tails
        # do not contain SVs
        ref_seq_surv = ref_seq[padding:len(ref_seq)-padding]
        # Write ref sequence to temporary fa file.
        temp_ref_fa = os.path.join(temp_dir, "temp_ref.fa")
        with open(temp_ref_fa, "w") as fh:
            add_fasta_entry(name, ref_seq_surv, fh)

        # Run SURVIVOR.
        prefix = os.path.join(temp_dir, "simulated")
        survivor_cmd = ["SURVIVOR",
                        "simSV",
                        temp_ref_fa,
                        param_file,
                        "0.0",
                        "0",
                        prefix]
        subprocess.check_call(survivor_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Read output of SURVIVOR
        altered_fa_path = "{}.fasta".format(prefix)
        insertions_fa_path = "{}.insertions.fa".format(prefix)
        sim_vcf = "{}.vcf".format(prefix)

        # Update VCF
        temp_vcf = os.path.join(temp_dir, "temp.vcf")
        update_vcf(temp_ref_fa, insertions_fa_path, sim_vcf, temp_vcf, pos_padding=padding)

        # Merge seqs and variants entries into single FA/VCF files
        # Add the initial and last 800bp back to the altered fasta
        altered_seq = pysam.FastaFile(altered_fa_path).fetch(name)
        altered_seq = update_altered_fa(ref_seq, altered_seq, padding)
        add_fasta_entry(name, altered_seq, out_altered_fa_fh)

        add_fasta_entry(name, ref_seq, out_ref_fa_fh)

        vcf_reader = vcf.Reader(open(temp_vcf, 'r'))
        if not out_vcf_fh:
            out_vcf_fh = vcf.Writer(open(out_vcf_path, 'w'), vcf_reader)

        for record in vcf_reader:
            out_vcf_fh.write_record(record)

        # Remove temporary files.
        import shutil
        shutil.rmtree(temp_dir)

    out_altered_fa_fh.close()
    out_ref_fa_fh.close()
    out_vcf_fh.close()

def find_survivor():
    try:
        path = subprocess.check_output(["which", "SURVIVOR"])
    except:
        print("ERROR: Please add SURVIVOR to your PATH variable and rerun. Thanks!")
        exit(1)

def main():
    args = parse_args()

    # Check for SURVIVOR.
    find_survivor()
    
    # Read the reference file
    ref = pysam.Fastafile(args.ref_file)
    

    # Make output dir
    os.mkdir(args.out_dir)

    regions = None

    if args.sv_regions:
        # Read sv_regions file, if provided
        regions = generate_regions_from_file(args.sv_regions)
    elif args.num_sv_regions:
        #Choose a random chromosome, a random region of 10kb within the chromosome
        regions = generate_random_regions(args.ref_file,
                                          args.len_sv_region,
                                          args.num_sv_regions)

    assert(regions is not None), "No regions to process. Please provide at least 1 region."

    process_regions(args.ref_file, regions, args.out_dir, args.survivor_param_file)


if __name__ == '__main__':
    main()
