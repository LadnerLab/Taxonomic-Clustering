#!/usr/bin/env python3

import optparse
import sys

import protein_oligo_library as oligo


def main():

    usage = "usage %prog [options]"
    option_parser = optparse.OptionParser( usage )

    add_program_options( option_parser )

    options, arguments = option_parser.parse_args()

    if options.query is None:
        print( "Fasta query file must be provided." )
        sys.exit() 
    elif options.tax is None:
        print( "Taxonomic lineage file must be provided." )
        sys.exit() 


def add_program_options( option_parser ):
    option_parser.add_option( '-q', '--query', help = "Fasta query file to read sequences from and do ordering of. [None, Required]" )
    option_parser.add_option( '-t', '--tax', help = "Taxonomic lineage file such as the one from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/" )

if __name__ == '__main__':
    main()
