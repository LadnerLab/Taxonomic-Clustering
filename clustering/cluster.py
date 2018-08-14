#!/usr/bin/env python3

class Cluster:
    def __init__( self, cluster_identifier):

        self.name = cluster_identifier
        self.names = list()
        self.sequences = list()
        self.kmers = set()
        self.kmer_dict = {}
        self.sequence_dict = {}
        self.kmer_size = 0
        self.sequence_size = 0

        # Minimum similarity between two sequences in a cluster
        self._least_similar_sequence = 1.0


    def _create_sequence_dict( self, names, sequence ):
        for index in range( len( names ) ):
            current_name = names[ index ]
            current_seq = sequences[ index ]

            self.sequence_dict[ current_name ] = current_seq

    def add_sequence( self, name, sequence ):
        self.names.append( name )
        self.sequences.append( sequence )
        self.sequence_dict[ name ] = sequence

        self.sequence_size += 1

    def add_sequence_kmers( self, name, new_kmers ):

        intersection = new_kmers & self.kmers
        percent_similar = ( len( intersection ) / len( new_kmers ) )

        self.kmers |= new_kmers
        self.kmer_dict[ name ] = new_kmers

        if self.sequence_size > 1 and percent_similar < self._least_similar_sequence:
            self._least_similar_sequence = percent_similar

    def add_sequence_and_its_kmers( self, name, sequence, kmer_set ):
        self.add_sequence( name, sequence )
        self.add_sequence_kmers( name, kmer_set )

    def set_cluster_size_threshold( self, new_thresh ):
        self.cluster_size_threshold = new_thresh

    def set_kmer_dict( self, kmer_dict ):
        self.kmer_dict = kmer_dict 

    def get_num_kmers( self ):
        self.kmer_size = len( self.kmers )
        return len( self.kmers )

    def get_num_sequences( self ):
        return self.sequence_size

    def get_names_and_sequences( self ):
        
        out_list = list()

        for current_name, current_seq in self.sequence_dict.items():
            out_list.append( ( current_name, current_seq ) )

        return out_list
           
    def remove_sequence( self, seq_name ):
        seq_to_remove = self.sequence_dict[ seq_name ]
        seq_to_remove_kmers = self.kmer_dict[ seq_name ]

        # Remove the sequence from our list of names and sequences
        self.names.remove( seq_name )
        self.sequences.remove( seq_to_remove )

        # Remove the sequence's kmers from our set of kmers
        self.kmers -= self.kmer_dict[ seq_name ]

        # Remove the sequence from our kmer dictionary and sequence dictionary
        del self.sequence_dict[ seq_name ]
        del self.kmer_dict[ seq_name ]

        self.sequence_size -= 1

        return seq_to_remove, seq_to_remove_kmers

    def write( self ):
        out_file = open( self.name + '.fasta', 'w' )

        out_file.write( str( self ) )
        out_file.close()

    @staticmethod
    def split_clusters_bigger_than_threshold( cluster_to_split, int_thresh ):
        
        out_clusters = list()
        current_cluster = None

        cluster_size = cluster_to_split.get_num_kmers()
        cluster_size_original = cluster_size

        original_cluster_name = cluster_to_split.name

        clustered_names = list()

        sub_cluster = 1

        while cluster_to_split.get_num_kmers() > int_thresh and cluster_to_split.sequence_size > 1:
            current_seq_to_remove_name = cluster_to_split.names[ 0 ]
            current_seq_to_remove, current_seq_kmers = cluster_to_split.remove_sequence(
                                                       current_seq_to_remove_name
                                                       )

            cluster_name = str( original_cluster_name ) + "_" + str( sub_cluster )

            if cluster_name not in clustered_names:
                new_cluster = Cluster( cluster_name ) 
                out_clusters.append( new_cluster )

                clustered_names.append( cluster_name )
                current_cluster = new_cluster

            current_cluster.add_sequence_and_its_kmers( current_seq_to_remove_name,
                                                        current_seq_to_remove,
                                                        current_seq_kmers
                                                      )

            if current_cluster.get_num_kmers() > int_thresh:
                sub_cluster += 1
        
        out_clusters.append( cluster_to_split )

        return out_clusters

    def __str__( self ):
        out_string = ""
        for name, sequence in self.sequence_dict.items():
            out_string += '> ' + name
            out_string += '\n'
            out_string += sequence
            out_string += '\n'
        return out_string.strip()
