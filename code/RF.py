# RF
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split, cross_val_score
import pandas as pd 
import numpy as np
import re
import os
import joblib
import argparse


if __name__ == '__main__':
    #The parameters are defined and encapsulated
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_phage', type=str, help = 'input phage key gene vector hash table csv file')
    parser.add_argument('--input_host', type=str, help = 'input host key gene vector hash table csv file')
    parser.add_argument('--input_PHI', type=str, help = 'input PHI excel file')
    parser.add_argument('--type', type=str, help = 'train or predict')
    parser.add_argument('--output', type=str, help = 'output file path')
    opt = parser.parse_args()     
        
    phage = pd.read_csv(opt.input_phage, index_col = 0, header = 0)
    host = pd.read_csv(opt.input_host, index_col = 0, header = 0)
    phi = pd.read_excel(opt.input_PHI, sheet_name = 'ID', header = 0)
    gene_name = phage.columns.tolist() + host.columns.tolist()
    fea = []
    label = []
    for i in range(len(phi)):
        if phi.iloc[i,0] in phage.index.tolist() and phi.iloc[i,1] in host.index.tolist():
            fea.append(np.concatenate((phage.loc[phi.iloc[i,0]].values, host.loc[phi.iloc[i,1]].values), axis = 0))
            label.append(phi.iloc[i,-1])
    fea = np.asarray(fea)
    
    train_data, test_data, train_label, test_label = train_test_split(fea, label, test_size = 0.3, random_state = 123)
    # define the model
    model = RandomForestClassifier(random_state = 123)
    if opt.type == 'train':
        # fit the model
        model.fit(train_data, train_label)
        joblib.dump(model, opt.output + os.sep + 'rf.pkl')
    else:
        model = joblib.load(opt.output + os.sep + 'rf.pkl')

    tra_label=model.predict(train_data) #train set predict label
    tes_label=model.predict(test_data) #test set predict label
    tes_prob = model.predict_proba(test_data)[:,1]
    print("train set：", accuracy_score(train_label,tra_label))
    print("test set：", accuracy_score(test_label,tes_label))
    print(cross_val_score(model,np.asarray(fea),np.asarray(label),scoring = 'roc_auc', cv=5).mean())
    # get importance
    importance = model.feature_importances_
    result = pd.DataFrame(index = list(range(len(gene_name))), columns = ['gene', 'importance'])
    result = result.replace(np.nan,0)
    result['gene'] = gene_name
    result['importance'] = importance
    result.to_csv(opt.output + os.sep + 'key_gene_importance.csv', index = True)
