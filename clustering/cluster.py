

class Cluster:
    def __init__( self, names, sequences ):
        self.names = names
        self.sequence = sequences
        self.kmer_dict = {}
        self.sequence_dict = {}


        self._create_sequence_dict( names, sequences )


    def _create_sequence_dict( self, names, sequence ):
        for index in range( len( names ) ):
            current_name = names[ index ]
            current_seq = sequences[ index ]

            self.sequence_dict[ current_name ] = current_seq






        

class TaxonomicCluster( Cluster ):
    def __init__( self, names, sequences ):
        super().__init( self, names, sequences )__


class KmerCluster( Cluster ):
    def __init__( self ):
        super().__init( self )__
