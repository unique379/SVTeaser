"""
Misc tools
"""
from acebinf import cmd_exe

def vcf_compress(fn):
    """
    Run vcftools to sort/compress/index a vcf file
    """
    ret = cmd_exe(f"vcf-sort {fn} | bgzip > {fn}.gz && tabix {fn}.gz")
