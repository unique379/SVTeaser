
# Main interface to generate structural variant simulations

import pysam
import pandas as pd
import argparse
import os
from random import randint

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
                     3) The SVs in a VCF file')

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
                              as needed by survivor.', required=True) 
    parser.add_argument('--out_dir', type=str,
                        help='Path to output directory.', required=True) 
    args = parser.parse_args()
    return args


def generate_region(ref, length):
    chroms = []
    for chrm in ref.references:
        chroms.append(chrm)
        
    chrmidx = randint(1, len(chroms))
    print (chrmidx)
    print (chroms)
    chrm = chroms[chrmidx -1]
    start = randint(1, ref.get_reference_length(chrm) - length)
    end = start + length

    # Offset by 1 since we started from 1 base for randint
    return chrm, start - 1, end - 1 


def main():
    args = parse_args()
    
    # Read the reference file
    ref = pysam.Fastafile(args.ref_file)
    

    # Make output dir
    os.mkdir(args.out_dir)

    # Read sv_regions file, if provided
    if args.sv_regions:
        regions = pd.read_csv(args.sv_regions, sep=", ")
        print (regions)
        for index, row in regions.iterrows():
            # 1. Fetch the region from reference
            chrm = row['chr']
            start = row['start']
            end = row['end']
            sv_region_ref = ref.fetch(chrm, start, end)
            # 2. Dump the unaltered region 
            region_name = chrm + "_" + str(start) + "_" + str(end)
            unaltered_regions_filename = os.path.join(args.out_dir, region_name + ".fa")

            with open(os.path.join(unaltered_regions_filename), "w") as rf:
                    rf.write(">" + region_name + "\n")
                    rf.write(sv_region_ref)

    elif args.num_sv_regions:
        #Choose a random chromosome, a random region of 10kb within the chromosome
        for idx in range(0, args.num_sv_regions):
            chrm, start, end = generate_region(ref, args.len_sv_region)
            sv_region_ref = ref.fetch(chrm, start, end)
            # 2. Dump the unaltered region
            region_name = chrm + "_" + str(start) + "_" + str(end)
            unaltered_regions_filename = os.path.join(args.out_dir, region_name + ".fa")

            with open(os.path.join(unaltered_regions_filename), "w") as rf:
                    rf.write(">" + region_name + "\n")
                    rf.write(sv_region_ref)


    # 3. Call Survivor
    #call_survivor()



if __name__ == '__main__':
    main()
