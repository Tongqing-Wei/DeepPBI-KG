#k_mean algorithm
import pandas as pd
import numpy as np
import os
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import StandardScaler 
from joblib import parallel_backend
import argparse

#The parameters are defined and encapsulated
parser = argparse.ArgumentParser()
parser.add_argument('--sample_file', type=str, help = 'input all negtive sample csv file')
parser.add_argument('--phage_feature', type=str, help = 'input phage feature vector file')
parser.add_argument('--host_feature', type=str, help = 'input host feature vector file')
parser.add_argument('--k', type=int, help = 'set kmeans with k cluster')
parser.add_argument('--output', type=str, help = 'output kmeans result')
opt = parser.parse_args()  


#Parameter initialization and sample feature generation
phage = pd.read_csv(opt.phage_feature, header = 0, index_col = 0)
host = pd.read_csv(opt.host_feature, header = 0, index_col = 0)
phi = pd.read_csv(opt.sample_file, sep = '\t', header = 0)
fea = []
for i in range(len(phi)):
    fea.append(np.concatenate((phage.loc[phi.iloc[i, 0],:].values, host.loc[phi.iloc[i, 1],:].values), axis = 0))
fea = np.asarray(fea)

inputfile = pd.DataFrame(fea) 
outputfile = opt.output + os.sep +  'data_type_train.csv' #save file result name
k = opt.k#The category of clustering
iteration = 3 #Maximum number of iteration of clustering
# data = pd.read_excel(inputfile, index_col = 'Id') #read data
data = inputfile
stdScale = StandardScaler().fit(data)
data_zs = pd.DataFrame(stdScale.transform(data)) #Data standardization, std() means to find the population sample variance (divided by n-1), in numpy std() is divided by n
data_zs_nan = np.isnan(data_zs)
data_zs[data_zs_nan] = 0
print(data_zs)


#batch_size = 100
#model = MiniBatchKMeans(n_clusters=k, max_iter = iteration, batch_size=batch_size, random_state=28) 
with parallel_backend('multiprocessing', n_jobs = 8):
    model = KMeans(n_clusters = k, max_iter = iteration, random_state = 123) #Divided into k class, the number of concurrent 8
    model.fit(data_zs) #Start clustering

#Simple print result
r1 = pd.Series(model.labels_).value_counts() #Count the number of categories
r2 = pd.DataFrame(model.cluster_centers_) #The center obtained by finding the clustering center r2 is not necessarily the original data point, but generally the center of mass of each cluster
r = pd.concat([r2, r1], axis = 1) #Horizontal concatenation (0 is vertical) to get the number under the category corresponding to the cluster center
print(r)
r.columns = list(data.columns) + [u'cluster categories'] #rename table head
r.to_csv(outputfile, index = False)

#Output raw data and its categories in detail
r = pd.concat([data, pd.Series(model.labels_, index = data.index)], axis = 1)  #Output the categories corresponding to each sample in detail
r.columns = list(data.columns) + [u'cluster categories'] #rename table head
# r.to_csv(outputfile, index = False) #save result

# Pick negative samples from each cluster
def distance(vec1,vec2) :
    return np.sqrt(sum(np.power(vec1-vec2,2)))
point = []
for i in range(k):
    d = data_zs[r[u'cluster categories'] == i]
    center = r2.iloc[i,:]
    ls = {}
    for j in range(len(d)):
        ls[d.iloc[j,:].name] = distance(np.asarray(d.iloc[j,:]), np.asarray(center))
    min_key = min(ls, key = ls.get)
    dic1SortList = sorted(ls.items(),key = lambda x:x[1],reverse = False)
    point.extend(pd.DataFrame(dic1SortList).iloc[:1,0].values.tolist())

data.loc[point].to_csv(opt.output + os.sep + 'negsample_train.csv', index = True)
