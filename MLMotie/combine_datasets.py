import pandas as pd


dfAgendapunt = pd.read_csv("Agendapunt.csv")
dfBesluit = pd.read_csv("Besluit.csv")

dfDocument = pd.read_csv("Document.csv")

dfAgendapunt = dfAgendapunt.rename(columns = {'Id':'Agendapunt_Id'})


dfAgendapuntBesluit = pd.merge(dfAgendapunt, dfBesluit, on="Agendapunt_Id")



dfAgendapuntBesluitDocument = pd.merge(dfAgendapuntBesluit, dfDocument, on = "Onderwerp")

print(len(dfAgendapuntBesluitDocument))
print(dfAgendapuntBesluitDocument.head(5))
print(dfAgendapuntBesluitDocument.columns)

dfAgendapuntBesluitDocument.to_csv("AgendapuntBesluitDocument.csv")