"""Simulate structural variants at a given region of genome."""

import pysam

class sv_simulator():
    """ Structural Variantion simulator"""

    def __init__(self, ref_file, num_sv_regions=50, max_sv_region_len=100):
        """Initialize class.

        Args:
            ref_file: Reference fasta file to generate structural variants on.
            num_sv_regions: Number of regions within the genome to add structural variants to
            max_sv_region_len: Maximum number of bp of structural variants to simulate.

        """
        self.genome_file = pysam.Fastafile(ref_file)

    def __init__(self, ref_file, sv_regions=None, max_sv_region_len=100):
        """Initialize class.

        Args:
            ref_file: Reference fasta file to generate structural variants on.
            sv_regions: A comma separated file containing start and end of the structural variant region.
                        1 structural variant will be generated per row of variant length defined by
                        randint(50, max_sv_region_len)
            max_sv_region_len: Maximum number of bp of structural variants to simulate.

        """
        self.genome_file = pysam.Fastafile(ref_file)

    def get_genome_file(self):
        return self.genome_file


