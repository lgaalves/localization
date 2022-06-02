#!/usr/bin/python
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd 

def correcting_index(df_raw):
    df=df_raw.copy()
    df.ix[3,"Intercountry Input-Output Table"]="Country"
    df.ix[1,"Intercountry Input-Output Table"]="Code"
    df.ix[2,"Intercountry Input-Output Table"]="xxx"
    df.ix[4,"Intercountry Input-Output Table"]="yyy"
    df.rename(columns={'Intercountry Input-Output Table':"Code", 'Unnamed: 2':"Country"},inplace=True)
    df.set_index(['Country', 'Code'],inplace=True)
    df=df.T
    indextuples=[tuple(x) for x in np.array(df[[(np.nan,"Code"),(np.nan,"Country")]])]
    df.set_index(pd.MultiIndex.from_tuples(indextuples, names=('Code','Country')),inplace=True)
    df=df.T
    df.dropna(inplace=True)
    df=df.T
    df=df.ix[df.index.dropna()]
    df=df.T
    return df

def save_edges(year):
    
    # Read label maps
    countries=pd.read_csv("../data/countries.csv",delimiter=",",index_col=None)
    sectors=pd.read_csv("../data/sectors.csv",delimiter="\t",index_col=None)
    
    # Remove $ and \
    listofcountries=list(countries.Abbreviation)
    listofcodes=list(sectors.Code)

    listofcodes[4]='C10-C12'
    listofcodes[5]='C13-C15'
    listofcodes[21]='C31_C32'
    listofcodes[25]='E37-E39'
    listofcodes[37]='J59_J60'
    listofcodes[39]='J62_J63'
    listofcodes[44]='M69_M70'
    listofcodes[48]='M74_M75'
    listofcodes[53]='R_S'

    countries.Abbreviation=listofcountries
    sectors.Code=listofcodes
    
    # Make multiple index to be used to select data
    indexlist=[]
    for country in listofcountries:
        for code in listofcodes:
            indexlist.append(tuple([country,code]))

    # Reads raw data
    df_raw=pd.read_excel("../data/{}/WIOT{}_Nov16_ROW.xlsx".format(year,year))
    
    # Correct index from raw data
    df=correcting_index(df_raw)
    
    # Select data based on list of countries and activities
    df=df.ix[indexlist].T.swaplevel().ix[indexlist].T
    
    # Make edge list
    data=[]
    for i in range(0,len(indexlist)):
        for j in range(0,len(indexlist)):
            print(i,j)
            data.append([indexlist[i][0],
                   indexlist[i][1],
                   indexlist[j][0],
                   indexlist[j][1],
                   df.ix[indexlist[i]][indexlist[j]]])
            
    edge_df=pd.DataFrame(data,columns=['source_node','source_layer','destination_node','destination_layer','weight'])  
    
    # Map names to numbers
    countriesmap=pd.Series([i for i in range(1,len(listofcountries)+1)],index=listofcountries)
    productsmap=pd.Series([i for i in range(1,len(listofcodes)+1)],index=listofcodes)

    edge_df.source_node=edge_df.source_node.map(countriesmap)
    edge_df.destination_node=edge_df.destination_node.map(countriesmap)

    edge_df.source_layer=edge_df.source_layer.map(productsmap)
    edge_df.destination_layer=edge_df.destination_layer.map(productsmap)

    # Save edge list
    data_path="../data/{}/multilayer/edges-{}.csv".format(year,year)
    edge_df.to_csv(data_path)
    
    return edge_df

def main():
    for year in range(2000,2015):
        save_edges(year)

if __name__ == "__main__":
    main()