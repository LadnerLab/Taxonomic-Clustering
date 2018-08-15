#!/usr/bin/env python3
import optparse
import sys
import os
import statistics
import copy

import protein_oligo_library as oligo
import clustering.cluster as cluster

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

    assert len( names ) == len( sequences )

    sequence_dict = { names[ index ]: sequences[ index ] for index in range( num_seqs ) }
    ymer_dict = {}
    total_ymers = set()

    for current_seq in range( len( sequences ) ):
        current_ymers = frozenset( oligo.subset_lists_iter( sequences[ current_seq ], options.kmerSize, 1 ) )
        total_ymers |= current_ymers
        ymer_dict[ names[ current_seq ] ] = current_ymers

    if 'tax' in options.clustering:
        if options.lineage:

            clusters_with_kmers = {}

            created_clusters = cluster_taxonomically( options, sequence_dict, ymer_dict )

        else:
            print( "Lineage file must be provided for taxonomic clustering, exiting" )
            sys.exit()
    else:
        sorted_ids = sorted( options.id.split( ',' ) )
        clusters_with_names, clusters_with_kmers, similar_clusters = cluster_by_kmers( float( sorted_ids[ 0 ] ), sequence_dict, ymer_dict )

        for current_id in sorted_ids[ 1:: ]:
            current_id = float( current_id )
            max_cluster_size = max( [ len( item ) for item in clusters_with_kmers.values() ] )

            if max_cluster_size > options.number:
                re_cluster_kmers( sequence_dict, ymer_dict, clusters_with_names, clusters_with_kmers, similar_clusters, current_id, options.number )

        print( "Id threshold: %s." % options.id )

        output_clusters = {}
        for cluster, names_list in clusters_with_names.items():
            output_clusters[ cluster ] = [ ( name, sequence_dict[ name ] ) for name in names_list ]
        clusters = output_clusters

    min_cluster_size, median_cluster_size, avg_cluster_size, max_cluster_size = get_cluster_stats( created_clusters, total_ymers )

    print( "Number of unique ymers: %d." % len( total_ymers ) )
    print( "Number of clusters: %d." % len( created_clusters.keys() ) )
    print( "Minimum cluster size: %d." % min_cluster_size )
    print( "Median cluster size: %.2f." % median_cluster_size )
    print( "Average cluster size: %.2f." % avg_cluster_size )
    print( "Maximum cluster size: %d." % max_cluster_size )

    write_outputs( options.output, created_clusters, options.number )

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

    overflow_clusters = 0
    cluster_size_file = open( 'cluster_sizes.txt', 'w' )

    for cluster_key, cluster_value in cluster_dict.items():

        names_list = [ item for item in cluster_value.names ]
        sequence_list = [ item for item in cluster_value.sequences ]
        sub_clusters = cluster.Cluster.split_clusters_bigger_than_threshold( cluster_value, threshold )

        num_lists = 0

        if len( sub_clusters ) > 0:
            num_lists = len( sub_clusters )
            for current_sub_cluster in sub_clusters:
                current_sub_cluster.write()

                sub_cluster_length = current_sub_cluster.get_num_kmers()
                cluster_size_file.write( current_sub_cluster.name + '.fasta|' + str( sub_cluster_length ) )
                cluster_size_file.write( '\n' )
        else:
            cluster_value.write()
            cluster_length = len( kmer_cluster_dict[ cluster_key ] )
            cluster_size_file.write( cluster_key + '.fasta|' + str( cluster_length ) )
            cluster_size_file.write( '\n' )

        if num_lists > 1:
            overflow_clusters += 1
            write_large_cluster( names_list, sequence_list, cluster_key )


    cluster_size_file.close()
    if overflow_clusters > 0:
        print( ( "WARNING: %d cluster(s) had more than %d kmers, and were split up. "
                 "The original large clusters were written to "
                 "%s/large_clusters."
               )
               % ( overflow_clusters, threshold, out_directory )
             )


def write_cluster( file_name, sequence_names, sequence_dict ):
    names_list = list()
    sequence_list = list()

    for current_seq in sequence_names:
        names_list.append( current_seq )
        sequence_list.append( sequence_dict[ current_seq ] )
    oligo.write_fastas( names_list, sequence_list, file_name + ".fasta" )
    

def sub_clusters_from_kmers( cluster_to_split, kmer_name_dict, names_list, sequence_list, int_thresh ):

    out_cluster_kmers = {}
    out_cluster_names = {}

    cluster_num = list( cluster_to_split.keys() )[ 0 ]
    cluster_size = cluster_to_split[ cluster_num ]
    cluster_size_original = cluster_size

    clustered_names = list()

    index = 0
    sub_cluster = 1

    while len( cluster_to_split[ cluster_num ] ) > int_thresh:
        cluster_to_split[ cluster_num ] -= kmer_name_dict[ names_list[ index ] ]
        cluster_name = str( cluster_num ) + "_" + str( sub_cluster )

        if cluster_name not in out_cluster_kmers:
            out_cluster_names[ cluster_name ] = list()
            out_cluster_kmers[ cluster_name ] = set()
        out_cluster_names[ cluster_name ].append( names_list[ index ] )
        out_cluster_kmers[ cluster_name ] |= kmer_name_dict[ names_list[ index ] ]

        clustered_names.append( names_list[ index ] )

        index += 1
        if len( out_cluster_kmers[ cluster_name ] ) > int_thresh:
            sub_cluster += 1

    if len( names_list ) - len( clustered_names ) != 0:
        cluster_name = str( cluster_num ) + "_" + str( sub_cluster + 1 )

        out_cluster_names[ cluster_name ] = list()
        out_cluster_kmers[ cluster_name ] = set()
        for index in range( len( names_list ) ):
            if names_list[ index ] not in clustered_names:
                out_cluster_names[ cluster_name ].append( names_list[ index ] )
                out_cluster_kmers[ cluster_name ] |= kmer_name_dict[ names_list[ index ] ]

    return out_cluster_names, out_cluster_kmers
    
    
    
def write_large_cluster( names_list, sequence_list, file_name ):
    dir_for_clusters = "large_clusters"
    if not os.path.exists( dir_for_clusters ):
        os.mkdir( dir_for_clusters )
    os.chdir( dir_for_clusters )
    
    oligo.write_fastas( names_list, sequence_list, str( file_name ) + "_too_large.fasta" )

    os.chdir( ".." )

def cluster_taxonomically( options, sequence_dict, kmer_dict ):
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
    sequence_dict = copy.deepcopy( sequence_dict )
    
    current_rank = oligo.Rank[ options.start[ 0 ].upper() ]
    if options.start is None:
        current_rank = oligo.Rank[ 'FAMILY' ]
    sequence_tax_id = set( [ oligo.get_taxid_from_name( item ) for item in names ] )

    tax_data = oligo.get_taxdata_from_file( options.lineage )
    created_clusters = {}
    clusters_created = list()

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
        deleted_clusters = list()

        if len( sequence_dict ) > 0:
            for current_name in list( sequence_dict.keys() ):

                current_id = int( oligo.get_taxid_from_name( current_name ) )
                current_id = check_for_id_in_merged_ids( merged_ids, current_id )

                current_rank_data = rank_data[ current_id ].lower()


                if current_id in rank_data and current_rank_data not in deleted_clusters:
                    if current_rank_data not in created_clusters:
                        new_cluster = cluster.Cluster( current_rank_data )
                        created_clusters[ current_rank_data ] = new_cluster

                    created_clusters[ current_rank_data ].add_sequence( current_name, sequence_dict[ current_name ] )
                    created_clusters[ current_rank_data ].add_sequence_kmers( current_name, kmer_dict[ current_name ] )

                    if created_clusters[ current_rank_data ].get_num_kmers() > options.number and index < len( ranks ) - 1:

                        # Put the items back in the pool of choices if our cluster becomes too large
                        put_large_cluster_back_in_pool( created_clusters, sequence_dict, current_rank_data )
                        deleted_clusters.append( current_rank_data )

                    else:
                            del sequence_dict[ current_name ]
                elif current_id not in rank_data:
                    print( "WARNING: An ID was not found in rank_data, this is likely to produce incorrect results" )

    return created_clusters

def cluster_by_kmers( id_threshold, sequence_dict, kmer_dict ):
    """
        Clusters sequences based on their number of shared kmers
    
        :param id: floating point identity threshold determining a sequence joins a cluster
        :param sequence_dict: dictionary of sequences containing name: sequence mappings 
        :param kmer_dict: dictionary of sequence name: kmer pairings

        :returns: dictionary of clusters created from sequences in sequence_dict
    """
    names_list = list( sequence_dict.keys() )
    sequence_list = list( kmer_dict.values() )

    out_clusters = {}

    names_list, sorted_seqs = oligo.sort_sequences_by_length( names_list,
                                                              sequence_list,
                                                              key = 'descending'
                                                            )

    cluster_name = "0_%f" % id_threshold

    out_clusters[ cluster_name ] = cluster.Cluster( cluster_name )
    out_clusters[ cluster_name ].add_sequence_and_its_kmers( names_list[ 0 ],
                                                             sorted_seqs[ 0 ],
                                                             kmer_dict[ names_list[ 0 ] ]
                                                           )


    for index in range( 1, len( sorted_seqs ) ):
        current_seq_name = names_list[ index ]
        current_seq_ymers = kmer_dict[ current_seq_name ]

        inserted = False

        if len( current_seq_ymers ) > 0:

            dict_items = out_clusters.items()
            for cluster_name, current_cluster in dict_items:

                intersection = current_seq_ymers & current_cluster.kmer_dict
                percent_similar = ( len( intersection ) / len( current_seq_ymers ) )

                if percent_similar >= id_threshold:
                    current_cluster.add_sequence_and_its_kmers( current_seq_name,
                                                                sequence_dict[ current_seq_name ],
                                                                current_seq_ymers
                                                              )


                    inserted = True
                    break
                
            if not inserted:
                cluster_number = len( kmer_clusters.keys() ) + 1
                cluster_name = "%d_%f" % ( cluster_number, id_threshold )

                out_clusters[ cluster_name ] = Cluster( cluster_name )
                out_cluster.add_sequence_and_its_kmers( current_seq_name,
                                                        sequence_dict[ current_seq_name ],
                                                        current_seq_ymers
                                                      )

    return out_clusters


def get_cluster_stats( cluster_dict, kmer_dict ):
    """
        Gets minimum, average and maximum cluster sizes from a dictionary of clusters
    
        :param cluster_dict: dictionary of clusters, where the key is the cluster label
        :param kmer_dict: dictionary of name: kmer pairings for sequences
    
        :returns: integer minimum, and maximum cluster sizes, float median and average cluster size
    """
    min_cluster_size = min( [ item.get_num_kmers() for item in cluster_dict.values() ] )
    median_cluster_size = 0
    avg_cluster_size = 0
    max_cluster_size = max( [ item.get_num_kmers() for item in cluster_dict.values() ] )

    total_size = 0
    num_clusters = len( cluster_dict.keys() )

    for key, current_cluster in cluster_dict.items():
        cluster_len = current_cluster.get_num_kmers()
        total_size += cluster_len
            
    avg_cluster_size = total_size / num_clusters
    median_cluster_size = statistics.median( sorted( [ item.get_num_kmers() for item in cluster_dict.values() ] ) )

    return min_cluster_size, median_cluster_size, avg_cluster_size, max_cluster_size
        
    
def re_cluster_kmers( sequence_dict, ymer_dict, clusters_with_names, clusters_with_kmers, similar_clusters, current_id, max_clust_size ):
    too_big_clusters = [ item for item in clusters_with_names.keys()
                                     if len( clusters_with_kmers[ item ] ) > max_clust_size
                       ]

    print( str( len( clusters_with_kmers[ '0_0.300000' ] ) ) )
    # least_similar_threshold = [  value for key, value in similar_clusters.items() 
    #                                  if len( clusters_with_kmers[ key ] ) > max_clust_size
    #                    ]
    least_similar_threshold = list()
    for key, value in similar_clusters.items():
        if len( clusters_with_kmers[ key ] ) > max_clust_size:
            least_similar_threshold.append( value )


    for current_cluster_index in range( len( too_big_clusters ) ):
        current_cluster = too_big_clusters[ current_cluster_index ]

        if least_similar_threshold[ current_cluster_index ] <= current_id:
            current_seq_dict = {}
            current_ymer_dict = {}

            max_key = str( max( [ int( item.split( '_' )[ 0 ] ) for item in list( clusters_with_names.keys() ) ] ) + 1 )
            max_key += str( "_%f" % current_id )

            sub_names = {}
            sub_kmers = {}
            for sequence_name in clusters_with_names[ current_cluster ]:
                current_seq_dict[ sequence_name ] = sequence_dict[ sequence_name ]
                current_ymer_dict[ sequence_name ] = ymer_dict[ sequence_name ]


            sub_clusters_with_names, sub_clusters_with_kmers, kmer_similarities = cluster_by_kmers( current_id,
                                                                                 current_seq_dict,
                                                                                 current_ymer_dict
                                                                               ) 

            sub_names.update( sub_clusters_with_names )
            sub_kmers.update( sub_clusters_with_kmers )

            for current_key in sub_names.keys():
                clusters_with_names[ max_key ] = sub_names[ current_key ]
                clusters_with_kmers[ max_key ] = sub_kmers[ current_key ]

                max_key = str( max( [ int( item.split( '_' )[ 0 ] ) for item in list( clusters_with_names.keys() ) ] ) + 1 )
                max_key += str( "_%f" % current_id )


            del clusters_with_names[ current_cluster ]
            del clusters_with_kmers[ current_cluster ]


def check_for_id_in_merged_ids( merged_ids, current_id ):
    if current_id in merged_ids:
        current_id = merged_ids[ current_id ]
    return current_id

def put_large_cluster_back_in_pool( clusters, sequence_dict, current_rank_data ):
    for item in clusters[ current_rank_data ].get_names_and_sequences():
        sequence_dict[ item[ 0 ] ] = item[ 1 ]

    del clusters[ current_rank_data ]
                        

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
    option_parser.add_option( '--id', default = '0.8', type = str,
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
