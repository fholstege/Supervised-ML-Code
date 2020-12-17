import pandas as pd
import requests
from cleaner import cleaner
import json 

df_BesluitAgendapuntDocument = pd.read_csv("BesluitAgendapuntDocumentDocumentActor.csv")

df_BesluitAgendapuntDocument = df_BesluitAgendapuntDocument.drop( 'Unnamed: 0', 1)
df_AgendapuntDocument = pd.json_normalize(df_BesluitAgendapuntDocument['Agendapunt'].apply(json.loads))

print(df_AgendapuntDocument)
print(df_AgendapuntDocument.columns)