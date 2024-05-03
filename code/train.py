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
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, classification_report, confusion_matrix
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc, mean_squared_error, precision_recall_curve
import argparse
import pickle
import time
start_time = time.time()

#set random seed
np.random.seed(0)
torch.manual_seed(0)
#torch.cuda.manual_seed_all(0)

# construct dataset
print("construct dataset")
class Dataset(Dataset):
    def __init__(self, input_data, input_label):
        data = []
        for i,j in zip(input_data, input_label):
            data.append((i,j))
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        d, l = self.data[index]
        return d, l
    
    
    
#define model
class DNN(nn.Module):
    def __init__(self):
        super(DNN,self).__init__()
        self.encoder  =  nn.Sequential(
            nn.Linear(input_dim,1024),
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


# define and plot confusion matrix
def get_confusion_matrix(trues, preds):
    labels = [0,1,2,3,4,5,6,7,8,9]
    conf_matrix = confusion_matrix(trues, preds)
    return conf_matrix
  
def plot_confusion_matrix(conf_matrix, path):
    plt.imshow(conf_matrix, cmap=plt.cm.Greens)
    indices = range(conf_matrix.shape[0])
    labels = [0,1]
    plt.xticks(indices, labels)
    plt.yticks(indices, labels)
    plt.colorbar()
    plt.xlabel('y_pred')
    plt.ylabel('y_true')
    # show data
    for first_index in range(conf_matrix.shape[0]):
        for second_index in range(conf_matrix.shape[1]):
            plt.text(first_index, second_index, conf_matrix[first_index, second_index])
    plt.savefig(path + os.sep + 'heatmap_confusion_matrix.png', dpi = 400)
    plt.show()    

if __name__ == '__main__':

    #The parameters are defined and encapsulated
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str, help = 'input train data file path')
    parser.add_argument('--test', type=str, help = 'input test data file path')
    parser.add_argument('--output', type=str, help = 'output relevant result file path')
    parser.add_argument('--epoch', type=int, help = 'train epoch')
    parser.add_argument('--lr', type=float, help = 'train learning rate')
    opt = parser.parse_args()         
    
    # read data
    #data = pd.read_excel('.../key_gene_train_random_set.xlsx', sheet_name = 'data', header = 0)
    data = pd.read_excel(opt.data, sheet_name = 'data', header = 0)
    #data = pd.read_excel('.../RBP_train_set.xlsx', sheet_name = 'data', header = 0)
    #data = pd.read_csv('.../key_CDD_kmeans.csv', header = 0)
    #test = pd.read_excel('.../key_gene_union_test_set.xlsx', sheet_name = 'data', header = 0)
    #test = pd.read_csv('.../key_CDD_kmeans.csv', header = 0)
    test = pd.read_excel(opt.test, sheet_name = 'data', header = 0)
    #test = pd.read_excel('...\\external\\case.xlsx', sheet_name = 'key_gene_PHI', header = 0)
    #test = pd.read_excel('.../RBP_union_test_set.xlsx', sheet_name = 'data', header = 0)
    #test = pd.read_excel('...\\external\\external_rbp_PHI.xlsx', sheet_name = 'data', header = 0)
    data


    print("loading data")
    #data, label = data.iloc[:,:-1], data.iloc[:,-1]
    train_data, train_label = data.iloc[:,:-1], data['label']
    test_data, test_label = test.iloc[:,:-1], test['label']
    #scaler = StandardScaler()
    #data = scaler.fit_transform(data)
    #train_data, test_data, train_label, test_label = train_test_split(data, label, test_size=.3, random_state=123)
    stdScale = StandardScaler().fit(train_data)
    train_data = stdScale.transform(train_data)
    test_data = stdScale.transform(test_data)
    pickle.dump(stdScale, open(opt.output + os.sep + 'scaler.pkl', 'wb'))
    train_data, test_data, train_label, test_label = np.asarray(train_data), np.asarray(test_data), np.asarray(train_label), np.asarray(test_label)
    print('train data：', train_data.shape)
    print('test data：', test_data.shape)


    # define relate parameter
    print("define relate parameter")
    epochs = opt.epoch
    batch_size = train_data.shape[0]
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu") 
    input_dim = train_data.shape[1]


    # DataLoader
    trainDataset = Dataset(train_data, train_label)
    testDataset = Dataset(test_data, test_label)
    # print(trainDataset[0])
    # print(trainDataset[0])
    trainDataLoader = DataLoader(trainDataset, batch_size=batch_size, shuffle=True, num_workers=0)
    testDataLoader = DataLoader(testDataset, batch_size=batch_size, shuffle=False, num_workers=0)


    # define loss function、optimizer and initialize relate parameter 
    model = DNN()
    print(model)
    model.to(device)

    print("define loss function、optimizer")
    criterion = nn.BCELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=opt.lr, weight_decay=1e-4)

    print("initialize relate parameter ")
    for param in model.parameters():
        nn.init.normal_(param, mean=0, std=0.01)



    # train and test
    print("start main iteration")
    model.train()
    for epoch in range(epochs):
        tot_loss = 0.0
        train_preds = []
        train_trues = []
        # model.train()
        for i,(train_data_batch, train_label_batch) in enumerate(trainDataLoader):
            train_data_batch = train_data_batch.float().to(device) # double data transform float
            train_label_batch = train_label_batch.to(device)
            outputs = model(train_data_batch)
            outputs = outputs.squeeze(dim=-1)
            # _, preds = torch.max(outputs.data, 1)
            loss = criterion(outputs, train_label_batch.float())
            # print(loss)
            #Backpropagation optimizes network parameters
            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
            #Add up the losses in each step
            tot_loss += loss.data
            train_outputs = outputs.ge(0.5).float()

            train_preds.extend(train_outputs.detach().cpu().numpy())
            train_trues.extend(train_label_batch.detach().cpu().numpy())


        sklearn_accuracy = accuracy_score(train_trues, train_preds) 
        sklearn_precision = precision_score(train_trues, train_preds, average='macro')
        sklearn_recall = recall_score(train_trues, train_preds, average='macro')
        sklearn_f1 = f1_score(train_trues, train_preds, average='macro')
        print("[sklearn_metrics] Epoch:{} loss:{:.4f} accuracy:{:.4f} precision:{:.4f} recall:{:.4f} f1:{:.4f}".format(epoch, tot_loss, sklearn_accuracy, sklearn_precision, sklearn_recall, sklearn_f1))
    
    
    test_preds = []
    test_trues = []
    test_probs = []
    torch.save(model.state_dict(), opt.output + os.sep + 'model.pth')
    model.eval()
    with torch.no_grad():
        for i,(test_data_batch, test_data_label) in enumerate(testDataLoader):
            test_data_batch = test_data_batch.float().to(device) # double data transform float
            test_data_label = test_data_label.to(device)
            test_outputs = model(test_data_batch)
            test_output = test_outputs.ge(0.5).float()
            test_preds.extend(test_output.detach().cpu().numpy())
            test_trues.extend(test_data_label.detach().cpu().numpy())
            test_probs.extend(test_outputs.detach().cpu().numpy())

        sklearn_accuracy = accuracy_score(test_trues, test_preds)
        sklearn_precision = precision_score(test_trues, test_preds, average='macro')
        sklearn_recall = recall_score(test_trues, test_preds, average='macro')
        sklearn_f1 = f1_score(test_trues, test_preds, average='macro')

        print(classification_report(test_trues, test_preds))
        conf_matrix = get_confusion_matrix(test_trues, test_preds)
        print(conf_matrix)
        plot_confusion_matrix(conf_matrix, opt.output)
        print("[sklearn_metrics] accuracy:{:.4f} precision:{:.4f} recall:{:.4f} f1:{:.4f}".format(sklearn_accuracy, sklearn_precision, sklearn_recall, sklearn_f1))




    precision, recall, _thresholds = precision_recall_curve(test_label,test_probs)
    roc_auc = auc(recall, precision)
    print("AUPR : ", roc_auc)

    # Binary ROC curve
    # roc_curve:（True Positive Rate , TPR）or（sensitivity）
    # x-axis：（False Positive Rate , FPR）
    fpr, tpr, thresholds_keras = roc_curve(test_trues, test_probs)
    roc_auc = auc(fpr, tpr)
    print("AUC : ", roc_auc)
    plt.subplots(figsize=(7,5.5));
    plt.plot(fpr, tpr, color='darkorange',
             lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc="lower right")
    plt.savefig(opt.output + os.sep + "ROC.png", dpi = 400)
    plt.show()

    print("--- %s seconds ---" % (time.time() - start_time))
