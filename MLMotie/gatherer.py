import requests
import pandas as pd


class gatherer:

    def __init__(self):

        self.api_root = "http://37.97.129.67/OData/v4/2.0/"
    
    def get_all_data(self, type_data):

        name_csv = ''.join(type_data)+ ".csv"

        index = 0
        query_size = 250

        keep_querying = True

        complete_data = self.request_data(type_data, index, query_size)
        index = index + query_size

        while keep_querying:

            partial_data = self.request_data(type_data,index,query_size)
            index = index + query_size

            complete_data = complete_data.append(partial_data)

            print(len(complete_data))

            if index % 50000 == 0:
                complete_data.to_csv(name_csv)

            if len(partial_data) < 250:
                keep_querying = False

        complete_data.to_csv(type_data + ".csv")
        return complete_data

    def request_data(self, type_data, index, query_size):

        type_data[0] + "?$expand="

        main_type = type_data[0]+ "?$expand="
        secondary_type = type_data[1]
        all_data_types = main_type + secondary_type

        if len(type_data)>2:
            for type_expansion in type_data[2:]:
                 all_data_types = all_data_types + "($expand=" + type_expansion
        
        all_data_types = all_data_types + ((len(type_data)-2) * ")")
        
        index_link = "&?$top=" + str(query_size) + "&$skip=" + str(index)

        self.link = self.api_root + all_data_types + index_link

        response = requests.get(self.link)
        response_dict = response.json()

        data = response_dict.get("value")

        return pd.DataFrame.from_dict(data)




