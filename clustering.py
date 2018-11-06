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

    sequence_dict = create_seq_dict( names, sequences )

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
        created_clusters = cluster_by_kmers( float( sorted_ids[ 0 ] ), sequence_dict, ymer_dict )

        for current_id in sorted_ids[ 1:: ]:
            current_id = float( current_id )
            max_cluster_size = max( [ item.get_num_kmers() for item in created_clusters.values() ] )

            if max_cluster_size > options.number:
                re_cluster_kmers( sequence_dict,
                                  ymer_dict,
                                  created_clusters,
                                  current_id,
                                  options.number
                                 )

        print( "Id threshold: %s." % options.id )

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

        num_lists = 0

        if cluster_value.get_num_kmers() > threshold:

            sub_clusters = cluster.Cluster.split_clusters_bigger_than_threshold( cluster_value, threshold )
            num_lists = len( sub_clusters )
            for current_sub_cluster in sub_clusters:
                current_sub_cluster.write()

                sub_cluster_length = current_sub_cluster.get_num_kmers()
                cluster_size_file.write( current_sub_cluster.name + '.fasta|' + str( sub_cluster_length ) )
                cluster_size_file.write( '\n' )
        else:
            cluster_value.write()
            cluster_length = cluster_value.get_num_kmers()
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
    file_name.replace( ' ', '_' )

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

    file_name = file_name.replace( ' ', '_' )
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

    names, sequences = seq_dict_to_names_and_seqs( sequence_dict )

    if options.start is None:
        start_ranks = [ 'FAMILY' ]
    else:
        start_ranks = options.start
    # Get the ranks descending order
    ranks = reversed( sorted(
                             [ oligo.Rank[ item.upper() ].value for item in start_ranks ]
                            )
                    )

    ranks = [ oligo.Rank( item ).name for item in ranks ]
    sequence_dict = copy.deepcopy( sequence_dict )

    current_rank = ranks[ 0 ]

    reference_names, reference_seqs = oligo.read_fasta_lists( options.unclustered )

    ref_dict          = create_seq_dict( reference_names, reference_seqs )

    # inverted dictionaries for missing id resolutions
    ref_dict_inverted = create_seq_dict( reference_names, reference_seqs, key = 'sequences' )
    seq_dict_inverted = create_seq_dict( names, sequences, key = 'sequences' )

    combined_dictionaries = combine_dicts( ref_dict_inverted, seq_dict_inverted )

    rank_map = oligo.parse_rank_map( options.rank_map )

    created_clusters = {}
    clusters_created = list()

    missing_seqs = list()

    sequence_tax_id = set()
    for current_name, current_seq in sequence_dict.items():
        if 'TaxID' not in current_name and 'OX' not in current_name:
            taxid  = resolve_missing_taxid( current_name, current_seq, combined_dictionaries )

            if taxid:
                sequence_tax_id.add( taxid )
            else:
                missing_seqs.append( ( current_name, current_seq ) ) 
        else:
            sequence_tax_id.add( oligo.get_taxid_from_name( current_name ) )


    tax_data = oligo.get_taxdata_from_file( options.lineage )
    tax_data = oligo.fill_tax_gaps( tax_data, rank_map )

    del reference_seqs


    merged_ids = {
                       10969   : 444185,
                       11619   : 2169991,
                       11630   : 2169993,
                       11806   : 353765,
                       45218   : 2169996,
                       45222   : 2169994,
                       45709   : 2169992,
                       489502  : 10407,
                       587201  : 10255,
                       587202  : 10255,
                       587203  : 10255,
                       1173522 : 11723,
                       1554474 : 1511807,
                       1554482 : 1330068,
                       1554483 : 1330491,
                       1554492 : 1330066,
                       1554494 : 1307800,
                       1554498 : 1511784,
                       1559366 : 1513237,
                       1560037 : 1131483,
                       2169701 : 11027
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

                current_id = oligo.get_taxid_from_name( current_name )
                if current_id:
                    current_id = int( current_id )
                else:
                    current_id = int( resolve_missing_taxid ) 

                current_id = check_for_id_in_merged_ids( merged_ids, current_id )

                if current_id:
                    current_rank_data = rank_data[ current_id ].lower()


                if current_id in rank_data and current_rank_data not in deleted_clusters:
                    if current_rank_data not in created_clusters:
                        new_cluster = cluster.Cluster( current_rank_data )
                        created_clusters[ current_rank_data ] = new_cluster

                    created_clusters[ current_rank_data ].add_sequence_and_its_kmers( current_name, sequence_dict[ current_name ], kmer_dict[ current_name ] )

                    if created_clusters[ current_rank_data ].get_num_kmers() > options.number and index < len( ranks ) - 1:

                        # Put the items back in the pool of choices if our cluster becomes too large
                        put_large_cluster_back_in_pool( created_clusters, sequence_dict, current_rank_data )
                        deleted_clusters.append( current_rank_data )

                    else:
                            del sequence_dict[ current_name ]
                elif current_id not in rank_data:
                    print( "WARNING: An ID was not found in rank_data, this is likely to produce incorrect results" )

    if missing_seqs:
        pass

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
    sequence_list = list( sequence_dict.values() )

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

                intersection = current_seq_ymers & current_cluster.kmers
                percent_similar = ( len( intersection ) / len( current_seq_ymers ) )

                if percent_similar >= id_threshold:
                    current_cluster.add_sequence_and_its_kmers( current_seq_name,
                                                                sequence_dict[ current_seq_name ],
                                                                current_seq_ymers
                                                              )


                    inserted = True
                    break
                
            if not inserted:
                cluster_number = len( out_clusters.keys() ) + 1
                cluster_name = "%d_%f" % ( cluster_number, id_threshold )

                out_clusters[ cluster_name ] = cluster.Cluster( cluster_name )
                out_clusters[ cluster_name ].add_sequence_and_its_kmers( current_seq_name,
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
        
def re_cluster_kmers( sequence_dict, ymer_dict, clusters, current_id, max_clust_size ):
    too_big_clusters = get_too_big_clusters( clusters, max_clust_size )

    for current_cluster_index in range( len( too_big_clusters ) ):
        current_cluster = too_big_clusters[ current_cluster_index ]
        least_similar_threshold = current_cluster.get_least_similar_sequence()

        if least_similar_threshold <= current_id:
            new_cluster = cluster.Cluster( current_cluster.name )

            max_key = str( max( [ int( item.name.split( '_' )[ 0 ] ) for item in list( clusters.values() ) ] ) + 1 )
            max_key += str( "_%f" % current_id )

            current_seq_dict = {}
            current_ymer_dict = {}
            sub_names = {}
            for sequence_name in clusters[ current_cluster.name ].names:
                current_seq_dict[ sequence_name ] = sequence_dict[ sequence_name ]
                current_ymer_dict[ sequence_name ] = ymer_dict[ sequence_name ]


            sub_clusters = cluster_by_kmers( current_id,
                                             current_seq_dict,
                                             current_ymer_dict
                                           ) 

            for current_sub_cluster in sub_clusters.values():
                current_name = current_sub_cluster.name
                max_key = str( max( [ int( item.name.split( '_' )[ 0 ] ) for item in list( clusters.values() ) ] ) + 1 )
                max_key += str( "_%f" % current_id )

                current_sub_cluster.name = max_key
                clusters[ max_key ] = current_sub_cluster

            del clusters[ current_cluster.name ]


def check_for_id_in_merged_ids( merged_ids, current_id ):
    if current_id in merged_ids:
        current_id = merged_ids[ current_id ]
    return current_id

def put_large_cluster_back_in_pool( clusters, sequence_dict, current_rank_data ):
    for item in clusters[ current_rank_data ].get_names_and_sequences():
        sequence_dict[ item[ 0 ] ] = item[ 1 ]

    del clusters[ current_rank_data ]
                        
def get_too_big_clusters( clusters, max_clust_size ):
    clusters_too_large = [
                           item for item in clusters.values() \
                           if item.get_num_kmers() > max_clust_size
                         ]
    return clusters_too_large
    
def get_least_similar_threshold_from_clusters( clusters_to_check ):
    return min( [ item.get_least_similar_sequence() for item in clusters_to_check.values() ] )
    

def get_repid_from_name( name ):

    try:
        repid = name.split( 'RepID=' )[ 1 ].strip()
    except:
        return None

    return repid

def create_seq_dict( names_list, seq_list, key = 'names' ):
    """
        Creates a dictionary from of list of names and a list of sequences,
        which can be either name: sequence mappings, or sequence: [ key ] mappings.
        :pre: names_list contains names, where the names[ i ] corresponds to sequences[ i ]
        :post: returns a dictionary containing the appropriate mapping, as specified by the 
               key parameter
    """
    return_dict = {}

    for index in range( len( names_list ) ):
        current_name = names_list[ index ]
        current_seq  = seq_list[ index ]

        if key == 'names':
            return_dict[ current_name ] = current_seq
        else:
            if current_seq not in return_dict:
                return_dict[ current_seq ] = list()
            return_dict[ current_seq ].append( current_name )
                
    return return_dict

def seq_dict_to_names_and_seqs( seq_dict, key = 'names' ):
    """
        Accepts a dictionary containing seq: name 
        or name: seq mappings, and returns a list containing
        the sequence names, and a list containing the sequences themselves
    
        :pre: seq_dict contains a mapping of key: value
              pairings, where key is a string and value is 
              either a string, list, or set.
        :post: returns a list of names and a list of sequences
               when key is 'names', names will be the first list 
               returned, and when key is 'sequences' 
               sequences will be the first list returned
    """
    names     = list()
    sequences = list()

    for current_name, current_seq in seq_dict.items():
        if type( current_seq ) == list or \
            type( current_seq ) == set:

            for item in current_seq:
                names.append( current_name )
                sequences.append( item )
        else:
            names.append( current_name )
            sequences.append( current_seq )
                
    if key == 'sequences':
        names, sequences = sequences, names

    return names, sequences

def combine_dicts( *dicts ):
    """
        Combines two or more dictionaries into one
        dictionary.

        :note: created dictionary will contain
               key: [ value ] mappings. If dictionaries
               share a key, then their values will be combined 
               into a list. 
    """
    return_dict = {}

    for current_dict in dicts:
        for key, value in current_dict.items():
            if key in return_dict:
                if type( return_dict[ key ] ) != list:
                    return_dict[ key ] = [ return_dict[ key ] ]
                return_dict[ key ].append( value )
            else:
                return_dict[ key ] = value if type( value ) == list else [ value ]
    return return_dict

def resolve_missing_taxid( name, sequence, combination_dict ):
    """
        Attempts to find a missing taxID where one cannot be found in name
    """
    return_id = ''
    names = combination_dict[ sequence ]

    for current in names:
        if 'RepID=' in current:
            rep_id = get_repid_from_name( current )

            if rep_id in current:
                return oligo.get_taxid_from_name( current )

        
    return return_id



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
    option_parser.add_option(
                              '-r', '--rank_map',
                              help = ( "Name of file containing mappings of "
                                       "taxids to level of taxonomic hierarchy."
                                     )
                            )
    option_parser.add_option( '-u', '--unclustered', 
                              help = "Name of unclusstered data file"
                            )

if __name__ == '__main__':
    main()
