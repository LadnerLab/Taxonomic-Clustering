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
            min_cluster_size, median_cluster_size, avg_cluster_size, max_cluster_size = get_cluster_stats( clusters, None )

            total_ymers = sequences 

            print( "Number of unique sequences: %d." % len( sequences ) )
            print( "Number of clusters: %d." % len( clusters.keys() ) )
            print( "Minimum cluster size: %d." % min_cluster_size )
            print( "Median cluster size: %.2f." % median_cluster_size )
            print( "Average cluster size: %.2f." % avg_cluster_size )
            print( "Maximum cluster size: %d." % max_cluster_size )

        else:
            print( "Lineage file must be provided for taxonomic clustering, exiting" )
            sys.exit()
    else:
        ymer_dict = {}
        sorted_ids = sorted( options.id.split( ',' ) )


        for current_seq in range( len( sequences ) ):
            current_ymers = frozenset( oligo.subset_lists_iter( sequences[ current_seq ], options.kmerSize, 1 ) )
            ymer_dict[ names[ current_seq ] ] = current_ymers

        clusters_with_names, clusters_with_kmers, total_ymers = cluster_by_kmers( float( sorted_ids[ 0 ] ), sequence_dict, ymer_dict )

        for current_id in sorted_ids[ 1:: ]:
            current_id = float( current_id )
            max_cluster_size = max( [ len( item ) for item in clusters_with_names.values() ] )

            if max_cluster_size > options.number:
                # Get the keys of the clusters that are too large
                re_cluster_kmers( sequence_dict, ymer_dict, clusters_with_names, clusters_with_kmers, current_id, options.number )

        min_cluster_size, median_cluster_size, avg_cluster_size, max_cluster_size = get_cluster_stats( clusters_with_kmers, total_ymers )

        print( "Id threshold: %s." % options.id )
        print( "Number of unique ymers: %d." % len( total_ymers ) )
        print( "Number of clusters: %d." % len( clusters_with_kmers.keys() ) )
        print( "Minimum cluster size: %d." % min_cluster_size )
        print( "Median cluster size: %.2f." % median_cluster_size )
        print( "Average cluster size: %.2f." % avg_cluster_size )
        print( "Maximum cluster size: %d." % max_cluster_size )

        output_clusters = {}
        for cluster, names_list in clusters_with_names.items():
            output_clusters[ cluster ] = [ ( name, sequence_dict[ name ] ) for name in names_list ]
        clusters = output_clusters

    write_outputs( options.output, clusters, clusters_with_kmers, sequence_dict, ymer_dict, options.number )

def write_outputs( out_directory, cluster_dict, kmer_cluster_dict, sequence_dict, kmer_name_dict, threshold ):
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

    overflow_clusters = 0

    for cluster_key, cluster_value in cluster_dict.items():
        names_list = [ item [ 0 ] for item in cluster_value ] 
        sequence_list = [ item[ 1 ] for item in cluster_value ]

        overflow = len( cluster_value ) // threshold
        num_lists = overflow if overflow > 0 else 1
        overflow = len( cluster_value ) % threshold 
        seqs_per_file = len( sequence_list ) // num_lists


        sub_clusters_from_kmers( { 3: kmer_cluster_dict[ 3 ] }, kmer_name_dict, names_list, sequence_list, threshold )
        if num_lists > 1:
            overflow_clusters += 1
            write_large_cluster( names_list, sequence_list, cluster_key )

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

    if overflow_clusters > 0:
        print( ( "WARNING: %d cluster(s) had more than %d sequences, and were split up "
                 "into sizes of %d. The original large clusters were written to "
                 "%s/large_clusters."
               )
               % ( overflow_clusters, threshold, threshold, out_directory )
             )


def sub_clusters_from_kmers( cluster_to_split, kmer_name_dict, names_list, sequence_list, int_thresh ):

    out_cluster_kmers = {}
    out_cluster_names = {}



    cluster_num = list( cluster_to_split.keys() )[ 0 ]
    cluster_size = len( cluster_to_split[ cluster_num ] )
    cluster_size_original = cluster_size

    index = 0
    sub_cluster = 2

    while cluster_size > int_thresh:
        cluster_size -= len( kmer_name_dict[ names_list[ index ] ] )

        cluster_name = str( cluster_num ) + "_" + str( sub_cluster )

        if cluster_name not in out_cluster_kmers:
            out_cluster_names[ cluster_name ] = list()
            out_cluster_kmers[ cluster_name ] = 0
        out_cluster_names[ cluster_name ].append( names_list[ index ] )
        out_cluster_kmers[ cluster_name ] += len( kmer_name_dict[ names_list[ index ] ] )

        if out_cluster_kmers[ cluster_name ] >= int_thresh:
            sub_cluster += 1
        
    return out_cluster_names
    
    
    
def write_large_cluster( names_list, sequence_list, file_name ):
    dir_for_clusters = "large_clusters"
    if not os.path.exists( dir_for_clusters ):
        os.mkdir( dir_for_clusters )
    os.chdir( dir_for_clusters )
    
    oligo.write_fastas( names_list, sequence_list, str( file_name ) + "_too_large.fasta" )

    os.chdir( ".." )

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

def cluster_by_kmers( id_threshold, sequence_dict, kmer_dict ):
    """
        Clusters sequences based on their number of shared kmers
    
        :param id: floating point identity threshold determining a sequence joins a cluster
        :param sequence_dict: dictionary of sequences containing name: sequence mappings 
        :param kmer_dict: dictionary of sequence name: kmer pairings

        :returns: dictionary of clusters created from sequences in sequence_dict
        :returns: set of all of the kmers from every sequence
    """
    names_list = list( sequence_dict.keys() )
    sequence_list = list( kmer_dict.values() )
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

        dict_items = kmer_clusters.items()
        for key, current_cluster in dict_items:
            intersection = current_seq_ymers & current_cluster
            percent_similar = ( len( intersection ) / len( current_seq_ymers ) )

            if percent_similar >= id_threshold:
                kmer_clusters[ key ] |= current_seq_ymers

                if key not in out_clusters:
                    out_clusters[ key ] = list()
                out_clusters[ key ].append( names_list[ index ] )
                inserted = True
                break
                
        if not inserted:
            cluster_number = len( kmer_clusters.keys() ) + 1
            
            kmer_clusters[ cluster_number ] = current_seq_ymers
            out_clusters[ cluster_number ] = [ names_list[ index ] ] 

    return out_clusters, kmer_clusters, total_kmers


def get_cluster_stats( cluster_dict, kmer_dict ):
    """
        Gets minimum, average and maximum cluster sizes from a dictionary of clusters
    
        :param cluster_dict: dictionary of clusters, where the key is the cluster label
        :param kmer_dict: dictionary of name: kmer pairings for sequences
    
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
        
    
def re_cluster_kmers( sequence_dict, ymer_dict, clusters_with_names, clusters_with_kmers, current_id, max_clust_size ):
    too_big_clusters = [ item for item in clusters_with_names.keys() \
                                     if len( clusters_with_names[ item ] ) > max_clust_size \
                       ]
    for current_cluster in too_big_clusters:
        current_seq_dict = {}
        current_ymer_dict = {}

        max_key = max( clusters_with_names.keys() ) + 1

        sub_names = {}
        sub_kmers = {}
        for sequence_name in clusters_with_names[ current_cluster ]:
            current_seq_dict[ sequence_name ] = sequence_dict[ sequence_name ]
            current_ymer_dict[ sequence_name ] = ymer_dict[ sequence_name ]


        sub_clusters_with_names, sub_clusters_with_kmers, total_ymers = cluster_by_kmers( current_id,
                                                                                  current_seq_dict,
                                                                                  current_ymer_dict
                                                                                ) 

        sub_names.update( sub_clusters_with_names )
        sub_kmers.update( sub_clusters_with_kmers )

        for current_key in sub_names.keys():
            clusters_with_names[ max_key ] = sub_names[ current_key ]
            clusters_with_kmers[ max_key ] = sub_kmers[ current_key ]

            max_key += 1


        del clusters_with_names[ current_cluster ]
        del clusters_with_kmers[ current_cluster ]


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
    option_parser.add_option( '--id', default = 0.8, type = str,
                              help = ( "Comma-separated list of identity thresholds to use for clustering. "
                                       "A sequence must share at least this proportion of its kmers with "
                                       "a cluster in order to join it. If a cluster is larger than the threshold "
                                       "specified by the --number flag, the next-biggest id will be used to break up "
                                       " the clusters larger than this number."
                                       " [0.8]"
                                     )
                            )
    option_parser.add_option( '-k', '--kmerSize', default = 10, type = int,
                              help = "Integer size of kmers to compare when doing clustering"
                            )

if __name__ == '__main__':
    main()
