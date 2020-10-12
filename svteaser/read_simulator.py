import os
import inspect
import argparse

from truvari import setup_logging
from acebinf import cmd_exe


def sim_reads_main(args):
    """
    Run read simulators
    """
    args = parseArgs(args)
    # Run the commands
    ret = cmd_exe("read_simulator -h")
    if ret.ret_code != 0:
        logging.error("Cannot find read_simulator in environment")
        exit(ret.ret_code)

    os.chdir(args.workdir)

    ret = cmd_exe("read_simulator")
    if ret.ret_code != 0:
        logging.error("Problem running read_simulator")
        logging.error(ret.stderr)
        exit(ret.ret_code)

    logging.info("Finished")

def parseArgs(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="sim_reads", description=inspect.getdoc(sim_reads_main),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("workdir", metavar="REF", type=str,
                        help="SVTeaser working directory")
    parser.add_argument("output", metavar="OUT", type=str,
                        help="Output directory for all the files")
    args = parser.parse_args(args)
    setup_logging()
    return args

