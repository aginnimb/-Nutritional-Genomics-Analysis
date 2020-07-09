import time
import csv
import pandas as pd
import collections

#read in the data files
chebi = pd.read_csv("/Users/anu/Desktop/chebi_pathways.csv")#, names = ["chebi_id","chebi_pathways"])
ncbi = pd.read_csv("/Users/anu/Desktop/ncbi_pathways.csv") #names = ["ncbi_id","ncbi_pathways"])

pathway_weights = collections.Counter()

def cleanse(df,col_name,id_col):
    """to clean the data and zip the columns as dictionary"""
    pathway_list= []
    for pathway in df[col_name]:
        clean = pathway.strip('}{').replace('"','').replace(' ','').split(',')
#         print("clean:",clean)
        pathway_weights.update(clean)
        #print("2nd pw:",pathway_weights)
        pathway_list.append(clean)
    dict_data = dict(zip(id_col,pathway_list))
    return dict_data

def jaccard_weight(x):
    """Function to calculate the weight of each pathway with in the list """
    return 1.0

def assoc_weight(x):
    """ Function to calculate the weighted scores for the pathway association """
    #print("rturn:", 1.0/pathway_weights[x])
    return 1.0/pathway_weights[x]

def similarity(l1, l2, w):
    """function to calculate the similarity scoring using Jaccard Index algorithm """
    s1 = set(l1)
    s2 = set(l2)
    numer = 0
    denom = 0
    for a in s1.intersection(s2):
        numer += w(a)
    for a in s1.union(s2):
        denom += w(a)
    if denom==0:
        return 0.0
    return numer/denom

def loop_and_score(self1,self2):
    """traverse through the list of lists(pathways) and jaccard scoring and writing to a csv file """
    header = ['GeneID','chebi_id','Jaccard_Similarity_Score','Weighted_Jaccard_score']
    with open("pathway_score.csv",'w') as infile:
        writer = csv.writer(infile,delimiter = ',')
        writer.writerow(header)
        for (keys1,vals1) in self1.items():
                for (keys2,vals2) in self2.items():
                        jscore = similarity(vals1, vals2, jaccard_weight)
                        #print("vals1:",vals1, "vals2:", vals2,"assoc_weight:",assoc_weight)
                        wscore = similarity(vals1, vals2, assoc_weight)
                        #print("wscore:",wscore)
                        if wscore >0.0:
#                             print(keys1,keys2,jscore,wscore)
                            writer.writerow([keys1,"CHEBI:" + str(keys2),jscore,wscore]) #adding prefix to second column using str()
    infile.close()

ncbi_dict = cleanse(ncbi,"Pathways",ncbi["GeneID"])
# print(ncbi_dict)
chebi_dict = cleanse(chebi,"Pathways",chebi["Chebi_id"])
loop_and_score(ncbi_dict,chebi_dict)
