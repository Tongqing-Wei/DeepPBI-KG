from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from itertools import cycle
import numpy as np
import pandas as pd
import re
import os
import torch
#torch.cuda.current_device() 
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, classification_report, confusion_matrix
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc, mean_squared_error, precision_recall_curve
import argparse
import pickle
import time
start_time = time.time()

# set random seed
np.random.seed(0)
torch.manual_seed(0)
#torch.cuda.manual_seed_all(0)
    
    
#define model
class DNN(nn.Module):
    def __init__(self):
        super(DNN,self).__init__()
        self.encoder  =  nn.Sequential(
            nn.Linear(2676,1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        x = self.encoder(x)
        return x


if __name__ == '__main__':

    #The parameters are defined and encapsulated
    parser = argparse.ArgumentParser()
    parser.add_argument('--phage', type=str, help = 'input phage feature data file path')
    parser.add_argument('--host', type=str, help = 'input host feature data file path')
    parser.add_argument('--output', type=str, help = 'output predict result file path')
    parser.add_argument('--scaler', type=str, help = 'input standScaler file path')
    parser.add_argument('--parameter', type=str, help = 'input model parameter file path')
    opt = parser.parse_args()         
    
    # read data
    phage = pd.read_csv(opt.phage, index_col = 0, header = 0)
    host= pd.read_csv(opt.host, index_col = 0, header = 0)


    print("Dataloader")
    phage_name = phage.index.tolist()
    host_name = host.index.tolist()
    
    interact = pd.DataFrame(index = list(range(len(phage) * len(host))), columns = ['phage', 'host'])
    interact['score'] = 0
    
    fea = []
    n = 0
    for i in range(len(phage)):
        for j in range(len(host)):
            fea.append(np.concatenate((phage.iloc[i,:].values, host.iloc[j, :].values), axis = 0))
            interact.iloc[n, 0] = phage_name[i]
            interact.iloc[n, 1] = host_name[j]
            n += 1
    fea = np.asarray(fea)

    scaler = pickle.load(open(opt.scaler, 'rb'))
    val_data = scaler.transform(fea)


    # loading model
    #device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = DNN()
    print(model)

    model.load_state_dict(torch.load((opt.parameter), map_location = torch.device("cpu")))
    
    result = model(torch.tensor(val_data).float())
    interact['score'] = result.detach().cpu().numpy().flatten()
    out = pd.DataFrame(index = phage_name, columns = host_name)
    for i in range(len(interact)):
        out.loc[interact.iloc[i,0], interact.iloc[i,1]] = interact.iloc[i,-1]
    out.to_csv(opt.output + os.sep + 'result.csv', index = True)
