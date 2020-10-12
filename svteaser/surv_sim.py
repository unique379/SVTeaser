import os
import inspect
import logging
import argparse

from collections import OrderedDict

from acebinf import cmd_exe
from truvari import setup_logging
from svteaser.vcfeditor import update_vcf
from svteaser.utils import vcf_compress

def edit_surv_params(fn):
    """
    Edit the SURVIVOR simSV parameters file
    """
    params = OrderedDict()
    with open(fn, 'r') as fh:
        for line in fh:
            data = line.strip().split(":")
            name = data[0]
            try:
                val = int(data[1].strip())
            except ValueError:
                val = data[1].strip()
            params[name] = val
    
    # We can control things here. 
    # might make sense to expose all of these to command line?
    params["TRANSLOCATION_number"] = 0
    params["INVERSION_number"] = 0
    with open(fn, 'w') as fout:
        for key, val in params.items():
            fout.write(f"{key}: {val}\n")

def surv_sim_main(args):
    """
    Run survivor simSV commands. Output them into a directory
    """
    args = parseArgs(args)
    # check the SURVIVOR is in the environment
    ret = cmd_exe("SURVIVOR -h")
    if ret.ret_code != 0:
        logging.error("Cannot find SURVIVOR in environment")
        exit(ret.ret_code)
    logging.debug(f"Making outdir {args.output}")
    try:
        os.mkdir(args.output)
    except FileExistsError:
        logging.error(f"Output directory {args.output} already exists")
        exit(1)

    os.chdir(args.output)

    logging.debug(f"Running SURVIVOR")
    ret = cmd_exe("SURVIVOR simSV surv_params")
    logging.debug(ret.stderr)
    logging.debug(ret.stdout)
    if ret.ret_code != 0:
        logging.error("Problem running SURVIVOR")
        logging.error(ret.stderr)
        exit(ret.ret_code)

    # Should actually do some parsing/editing of the surv_params
    edit_surv_params("surv_params")

    ret = cmd_exe(f"SURVIVOR simSV {args.reference} surv_params 0 0 simulated.sv")
    logging.debug(ret.stderr)
    logging.debug(ret.stdout)
    if ret.ret_code != 0:
        logging.error("Problem running SURVIVOR")
        logging.error(ret.stderr)
        exit(ret.ret_code)
    
    /usr/local/lib/python3.7/site-packages/svteaser/surv_sim.py

    update_vcf(args.reference, "simulated.sv.insertions.fa", "simulated.sv.vcf", "simulated.sv.vcf")
    vcf_compress("simulated.sv.vcf")
    logging.info("Finished")

def parseArgs(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="surv_sim", description=inspect.getdoc(surv_sim_main),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("reference", metavar="REF", type=str,
                        help="Reference file overwhich to simulate SVs")
    parser.add_argument("output", metavar="OUT", type=str, default="output",
                        help="SVTeaser output basename (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    args.reference = os.path.abspath(args.reference)
    args.output = args.output + ".svt"
    setup_logging(args.debug)
    return args
