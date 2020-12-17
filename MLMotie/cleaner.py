import pandas as pd



class cleaner:
    
    def get_moties(self, dfDocuments, kamer = 2.0):

        motie_values = ['Motie (gewijzigd/nader)', 'Motie']

        motie_boolean = dfDocuments["Soort"].isin(motie_values)

        dfMoties = dfDocuments[motie_boolean]
        dfMoties = dfMoties[dfMoties['Kamer'] == kamer]

        return dfMoties
    
    def get_stem_results(self, dfBesluit):

        stem_values = ['Stemmen - aangenomen', 
                                'Stemmen - verworpen', 
                                'Stemmen - uitstellen', 
                                'Stemmen - ingetrokken', 
                                'Stemmen - aangehouden', 
                                'Stemmen - zonder stemming aannemen']

        stem_boolean = dfBesluit["BesluitSoort"].isin(stem_values)

        dfStemmingen = dfBesluit[stem_boolean]

        return dfStemmingen