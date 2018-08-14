#!/usr/bin/env python3
from enum import enum

class Rank( Enum ):
        SUPERKINGDOM = 8
        KINGDOM = 7
        PHYLUM = 6
        CLASS = 5
        ORDER = 4
        FAMILY = 3
        GENUS = 2
        SPECIES = 1
        TAX_NAME = 0

class Cluster:
    def __init__( self, cluster_identifier):

        self.name = cluster_identifier
        self.names = list()
        self.sequence = list()
        self.kmers = set()
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
        self.sequence.append( sequence )
        self.sequence_dict[ name ] = sequence

        self.sequence_size += 1

    def add_sequence_kmers( self, new_kmers ):
        self.kmers |= new_kmers

    def set_cluster_size_threshold( self, new_thresh ):
        self.cluster_size_threshold = new_thresh

    def set_kmer_dict( self, kmer_dict ):
        self.kmer_dict = kmer_dict 

    def get_num_kmers( self ):
        return len( self.kmers )

    def get_num_sequences( self ):
        return self.sequence_size

    def get_names_and_sequences( self ):
        
        out_list = list()

        for current_name, current_seq in self.sequence_dict.items():
            out_list.append( ( current_name, current_seq ) )

            

