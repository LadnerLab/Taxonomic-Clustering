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
    def __init__( self, names, sequences ):
        self.names = names
        self.sequence = sequences
        self.kmer_dict = {}
        self.sequence_dict = {}
        self.cluster_size_threshold = 10000


        self._create_sequence_dict( names, sequences )

   



    def _create_sequence_dict( self, names, sequence ):
        for index in range( len( names ) ):
            current_name = names[ index ]
            current_seq = sequences[ index ]

            self.sequence_dict[ current_name ] = current_seq

    def set_cluster_size_threshold( self, new_thresh ):
        self.cluster_size_threshold = new_thresh


class TaxonomicCluster( Cluster ):
    def __init__( self, names, sequences, lineage_file = None ):
        super().__init( self, names, sequences )__

        self.lineage_file = lineage_file


    


class KmerCluster( Cluster ):
    def __init__( self ):
        super().__init( self, names, sequences )__
