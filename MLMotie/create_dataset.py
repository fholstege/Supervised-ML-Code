import pandas as pd


from gatherer import gatherer
from cleaner import cleaner

gatherer = gatherer()
cleaner = cleaner()

# set to max width
pd.set_option('display.max_columns', None)

# get df with moties
#dfDocuments = gatherer.get_all_data("Activiteit")
#dfMoties = cleaner.get_moties(dfDocuments)

data_types = ["Besluit", "Agendapunt", "Document", "DocumentActor"]

gatherer.get_all_data(data_types)