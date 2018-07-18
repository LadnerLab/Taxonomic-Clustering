#!/usr/bin/env python3

import optparse
import sys
import os
import statistics

import protein_oligo_library as oligo

def main():

    usage = "usage %prog [options]"
    option_parser = optparse.OptionParser( usage )

    add_program_options( option_parser )

    options, arguments = option_parser.parse_args()

    if options.query is None:
        print( "Fasta query file must be provided." )
        sys.exit() 

    names, sequences = oligo.read_fasta_lists( options.query )
    num_seqs = len( names )
    sequence_dict = { names[ index ]: sequences[ index ] for index in range( num_seqs ) }

    if 'tax' in options.clustering:
        if options.lineage:
            clusters = cluster_taxonomically( options, sequence_dict )
        else:
            print( "Lineage file must be provided for taxonomic clustering, exiting" )
            sys.exit()
    else:
        clusters = {}
        clusters_ymers = {}
        ymer_dict = {}

        for current_seq in range( len( sequences ) ):
            current_ymers = oligo.subset_lists_iter( sequences[ current_seq ], 10, 1 )
            clusters_ymers[ names[ current_seq ] ] = current_ymers
            ymer_dict[ names[ current_seq ] ] = current_ymers

        clusters_with_names, total_ymers = cluster_by_kmers( options, sequence_dict, ymer_dict )

        for key, value in clusters_with_names.items():
            clusters[ key ] = [ ( current_name, sequence_dict[ current_name ] ) for current_name in value ]
            for current_name in value:
                if key not in clusters_ymers:
                    clusters_ymers[ key ] = set()
                clusters_ymers[ key ] = clusters_ymers[ key ] | ymer_dict[ current_name ]
                 


    min_cluster_size, median_cluster_size, avg_cluster_size, max_cluster_size = get_cluster_stats( clusters_ymers, total_ymers )
    print( "Number of unique ymers: %d." % len( total_ymers ) )
    print( "Number of clusters: %d." % len( clusters.keys() ) )
    print( "Minimum cluster size: %d." % min_cluster_size )
    print( "Median cluster size: %.2f." % median_cluster_size )
    print( "Average cluster size: %.2f." % avg_cluster_size )
    print( "Maximum cluster size: %d." % max_cluster_size )

    write_outputs( options.output, clusters, options.number )
        

def write_outputs( out_directory, cluster_dict, threshold ):
    """
        Writes program outputs to directory specified by the output option on the command line

        :param out_directory: directy to write output files to, created if it does not exist
        :param cluster_dict: dictionary of clusters, one cluster is written to each file
        :param threshold: numerical threshold, if a cluster contains more sequences than this number,
                          a new file is written in the format Rank_2_.fasta
        
    """
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
                                str( cluster_key ) + "_" + str( index + 1 ) + "_.fasta"
                              )
            start += seqs_per_file + overflow
            end += seqs_per_file 
            overflow = 0


def cluster_taxonomically( options, sequence_dict ):
    """
        Clusters sequences based on taxonomic rank. Ranks that have more than 
        the options.number amount of sequences will be split up evenly. 

        :param options: options object to get the program preferences from.
        :param sequence_dict: dictionary o fsequences with name: sequence key-value pairs
    
        :returns clusters: dictionary of cluster: sequence pairings
    """

    names = sequence_dict.keys()
    # Get the ranks descending order
    ranks = reversed( sorted(
                             [ oligo.Rank[ item.upper() ].value for item in options.start ]
                            )
                    )

    ranks = [ oligo.Rank( item ).name for item in ranks ]
    
    current_rank = oligo.Rank[ options.start[ 0 ].upper() ]
    if options.start is None:
        current_rank = oligo.Rank[ 'FAMILY' ]
    sequence_tax_id = set( [ oligo.get_taxid_from_name( item ) for item in names ] )

    tax_data = oligo.get_taxdata_from_file( options.lineage )
    clusters = {}

    merged_ids = { 10969 : 444185, 11619 : 216991, 11630 : 216993,
                   11806 : 353765, 45218 : 2169996, 45222 : 2169994,
                   45709 : 2169992
                 }

    current_rank = oligo.Rank[ ranks [ len( ranks ) - 1 ] ].value 
    index = 0
    for index in range( len( ranks ) ):
        current_rank = oligo.Rank[ ranks[ index ] ].value
        rank_data = oligo.group_seq_from_taxid( sequence_tax_id,
                                                merged_ids,
                                                tax_data,
                                                current_rank
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
                    if len( clusters[ current_rank_data ] ) > options.number and index < len( ranks ) - 1:
                        # Put the items back in the pool of choices if our cluster becomes too large
                        for item in clusters[ current_rank_data ]:
                            sequence_dict[ item[ 0 ] ] = item[ 1 ]

                        del clusters[ current_rank_data ]
                    else:
                            del sequence_dict[ current_name ]
                else:
                    print( "WARNING: An ID was not found in rank_data, this is likely to produce incorrect results" )

    return clusters

def cluster_by_kmers( options, sequence_dict, kmer_dict ):
    """
        Clusters sequences based on their number of shared kmers
    
        :param options: options object containing the necessary data
        :param sequence_dict: dictionary of sequences containing name: sequence mappings 

        :returns: dictionary of clusters created from sequences in sequence_dict
    """
    names_list = list( sequence_dict.keys() )
    sequence_list = list( sequence_dict.values() )
    kmer_clusters = {}
    out_clusters = {}
    total_kmers = set()


    names_list, sorted_seqs = oligo.sort_sequences_by_length( names_list, sequence_list, key = 'descending' )

    kmer_clusters[ 0 ] = kmer_dict[ names_list[ 0 ] ]
    out_clusters[ 0 ] = [ names_list[ 0 ] ]

    for index in range( 1, len( sorted_seqs ) ):
        current_seq_ymers = kmer_dict[ names_list[ index ] ]
        inserted = False
        total_kmers |= current_seq_ymers

        for current_cluster in list( kmer_clusters.keys() ):
            intersection = current_seq_ymers & kmer_clusters[ current_cluster ]
            percent_similar = ( len( intersection ) / len( current_seq_ymers ) )

            if percent_similar >= options.id:
                kmer_clusters[ current_cluster ] = \
                               ( current_seq_ymers | kmer_clusters[ current_cluster ] )
                if current_cluster not in out_clusters:
                    out_clusters[ current_cluster ] = list()
                out_clusters[ current_cluster ].append( names_list[ index ] )
                inserted = True
                break
                
        if not inserted:
            kmer_clusters[ index ] = current_seq_ymers
            out_clusters[ index ] = [ names_list[ index ] ] 

    return out_clusters, total_kmers


def get_cluster_stats( cluster_dict, kmer_dict ):
    """
        Gets minimum, average and maximum cluster sizes from a dictionary of clusters
    
        :param cluster_dict: dictionary of clusters, where the key is the cluster label
    
        :returns: integer minimum, and maximum cluster sizes, float median and average cluster size
    """
    min_cluster_size = len( list( cluster_dict.values() ) [ 0 ] )
    median_cluster_size = 0
    avg_cluster_size = 0
    max_cluster_size = len( list( cluster_dict.values() ) [ 0 ] )

    total_size = 0
    num_clusters = len( cluster_dict.keys() )

    for key, current_cluster in cluster_dict.items():
        cluster_len = len( current_cluster ) 
        total_size += cluster_len

        if cluster_len < min_cluster_size:
            min_cluster_size = cluster_len
        if cluster_len > max_cluster_size:
            max_cluster_size = cluster_len
            
    avg_cluster_size = total_size / num_clusters
    median_cluster_size = statistics.median( sorted( [ len( item ) for item in cluster_dict.values() ] ) )

    return min_cluster_size, median_cluster_size, avg_cluster_size, max_cluster_size
        
    
def add_program_options( option_parser ):
    option_parser.add_option( '-q', '--query', help = "Fasta query file to read sequences from and do ordering of. [None, Required]" )
    option_parser.add_option( '-l', '--lineage', help = "Taxonomic lineage file such as the one from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/" )
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
    option_parser.add_option( '-c', '--clustering', default = 'kmer',
                              help = ( "Method to use for clustering. Can be taxonomic or kmer-based. If taxonomic is selected, "
                                       "a taxonomic lineage file must also be provided. No lineage file is necessary for kmer "
                                       "clustering method. [kmer]"
                                     )
                            )
    option_parser.add_option( '--id', default = 0.8, type = float,
                              help = ( "Percentage of its kmers a sequence must share with a "
                                       "cluster in order for it to become a member of that cluster"
                                       " [0.8]"
                                     )
                            )

if __name__ == '__main__':
    main()
