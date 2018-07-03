#!/usr/bin/env python3

import optparse
import sys
import os

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

    # Get the ranks descending order
    ranks = reversed( sorted(
                             [ oligo.Rank[ item.upper() ].value for item in options.start ]
                            )
                    )

    ranks = [ oligo.Rank( item ).name for item in ranks ]
    
    current_rank = oligo.Rank[ options.start[ 0 ].upper() ]
    if options.start is None:
        current_rank = oligo.Rank[ 'FAMILY' ]



    names, sequences = oligo.read_fasta_lists( options.query )
    num_seqs = len( names )
    sequence_dict = { names[ index ]: sequences[ index ] for index in range( num_seqs ) }

    sequence_tax_id = set( [ oligo.get_taxid_from_name( item ) for item in names ] )

    tax_data = oligo.get_taxdata_from_file( options.tax )
    clusters = {}

    merged_ids = { 10969 : 444185, 11619 : 216991, 11630 : 216993,
                   11806 : 353765, 45218 : 2169996, 45222 : 2169994,
                   45709 : 2169992
                 }

    for index in range( current_rank.value, -1, -1 ):
        rank_data = oligo.group_seq_from_taxid( sequence_tax_id,
                                                merged_ids,
                                                tax_data,
                                                index
                                              )
        if len( sequence_dict ) > 0:
            for current_name in list( sequence_dict.keys() ):
                current_id = int( oligo.get_taxid_from_name( current_name ) )

                if current_id in merged_ids:
                    current_id = merged_ids[ current_id ]

                if current_id in rank_data:
                    current_rank_data = rank_data[ current_id ]
                    if current_rank_data not in clusters:
                        clusters[ current_rank_data ] = list()

                    clusters[ current_rank_data ].append( ( current_name, sequence_dict[ current_name ] ) )
                    del sequence_dict[ current_name ]


    write_outputs( options.output, clusters, options.number )
                   
                
        

def write_outputs( out_directory, cluster_dict, threshold ):
    if not os.path.exists( out_directory ):
        os.mkdir( out_directory )
    os.chdir( out_directory )

    for cluster_key, cluster_value in cluster_dict.items():
        names_list = [ item [ 0 ] for item in cluster_value ] 
        sequence_list = [ item[ 1 ] for item in cluster_value ]

        overflow = len( cluster_value ) // threshold
        num_lists = overflow if overflow > 0 else 1
        overflow = len( cluster_value ) % threshold 
        seqs_per_file = len( sequence_list ) // num_lists

        start = 0
        end = seqs_per_file + overflow
        for index in range( num_lists ):
            oligo.write_fastas( names_list[ start:end ],
                                sequence_list[ start:end ],
                                cluster_key + "_" + str( index + 1 ) + "_.fasta"
                              )
            start += seqs_per_file + overflow
            end += seqs_per_file 
            overflow = 0

def add_program_options( option_parser ):
    option_parser.add_option( '-q', '--query', help = "Fasta query file to read sequences from and do ordering of. [None, Required]" )
    option_parser.add_option( '-t', '--tax', help = "Taxonomic lineage file such as the one from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/" )
    option_parser.add_option( '-n', '--number', type = int, default = 10000,
                              help = "Threshold value for determining cutoff of number of sequences that can be included in each output. [10,000]"
                            )
    option_parser.add_option( '-s', '--start', action = "append", 
                              help = ( "Level of the taxonomic hierarchy at which to begin "
                                       "clustering. If this option is given multiple times, "
                                       "e.g. -s family -s phylum, "
                                       "they will be processed in order of taxonomic rank, e.g., "
                                       "superkingdom, kingdom, phylum, class, order, family, genus, species [ family ]"

                                     )
                            )
    option_parser.add_option( '-o', '--output', default = 'tax_out',
                              help = "Directory to write grouped fasta files to, each file contains one rank-level grouping"
                            )
if __name__ == '__main__':
    main()
