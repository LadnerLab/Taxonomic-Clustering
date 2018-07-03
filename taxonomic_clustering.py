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

    current_rank = oligo.Rank[ options.start.upper() ]
    if options.start is None:
        current_rank = oligo.Rank[ 'FAMILY' ]


    names, sequences = oligo.read_fasta_lists( options.query )
    num_seqs = len( names )
    sequence_dict = { names[ index ]: sequences[ index ] for index in range( num_seqs ) }

    sequence_tax_id = set( [ oligo.get_taxid_from_name( item ) for item in names ] )

    tax_data = oligo.get_taxdata_from_file( options.tax )
    clusters = {}

    for index in range( current_rank.value, -1, -1 ):
        rank_data = oligo.group_seq_from_taxid( sequence_tax_id,
                                                tax_data,
                                                index
                                              )
        if len( sequence_dict ) > 0:
            for current_name in names:
                current_id = int( oligo.get_taxid_from_name( current_name ) )
                if current_id in rank_data:
                    current_rank_data = rank_data[ current_id ]
                    if not current_rank_data in clusters:
                        clusters[ current_rank_data ] = list()

                    clusters[ current_rank_data ].append( ( current_name, sequence_dict[ current_name ] ) )
                    del sequence_dict[ current_name ]
        

def add_program_options( option_parser ):
    option_parser.add_option( '-q', '--query', help = "Fasta query file to read sequences from and do ordering of. [None, Required]" )
    option_parser.add_option( '-t', '--tax', help = "Taxonomic lineage file such as the one from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/" )
    option_parser.add_option( '-n', '--number', type = int, default = 10000,
                              help = "Threshold value for determining cutoff of number of sequences that can be included in output. [10,000]"
                            )
    option_parser.add_option( '-s', '--start', default = 'family',
                              help = "Level of the taxonomic hierarchy at which to begin clustering. [family]"
                            )

if __name__ == '__main__':
    main()
