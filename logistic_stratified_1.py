#from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
#from sklearn.datasets import load_iris
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import classification_report
from sklearn.model_selection import StratifiedKFold

#X, y = load_iris(return_X_y=True)
#print(y)
#print(X)
#print(y)
#data = pd.read_csv('TumorOnly_Training_Overlap_compared_RNA_WES_RNA_pipeline_data_10_26.txt', sep=" ", header=None)
###############Training#######################
with open("TumorOnly_Training_Overlap_compared_RNA_WES_RNA_pipeline_data_10_26.txt", "r") as r1:
    list_sample=[]
    list_chrom=[]
    list_POS=[]
    list_ref=[]
    list_alt=[]
    list_ref_reads=[]
    list_alt_reads=[]
    list_VAF=[]
    list_status=[]
    dict_DNA={}
    dict_DNA["A"] = [1, 0, 0, 0]
    dict_DNA["T"] = [0, 1, 0, 0]
    dict_DNA["G"] = [0, 0, 1, 0]
    dict_DNA["C"] = [0, 0, 0, 1]
    count=0
    for lines in r1:
        count+=1
        if count>2:
            if "chrX" in lines:
                list_chrom+=[99]
            else:
                list_chrom+=[int(lines.split("\t")[1].strip("chr"))]
            list_POS+=[int(lines.split("\t")[2])]
            temp_list=[]
            #temp_list2=[]
            for item in lines.split("\t")[3]:
                #print(lines)
                temp_list+=dict_DNA[item]
            temp_list2=0
            for x in range(len(temp_list)):
                if temp_list[x]==1:
                    temp_list2+=2**(len(temp_list)-x)
            #print(temp_list2)
            list_ref+=[temp_list2]
            temp_list=[]
            temp_list2=[]
            for item in lines.split("\t")[4]:
                temp_list+=dict_DNA[item]
            temp_list2=0
            for x in range(len(temp_list)):
                if temp_list[x]==1:
                    temp_list2+=2**(len(temp_list)-x)
            #if temp_list2>200:
                #print(temp_list2)
            list_alt+=[temp_list2]
            #print(list_alt)
            #list_ref+=[lines.split("\t")[3]]
            #list_alt+=[lines.split("\t")[4]]
            list_ref_reads+=[float(lines.split("\t")[5])]
            list_alt_reads+=[float(lines.split("\t")[6])]
            list_VAF+=[float(lines.split("\t")[7])]
            #print(lines)
            list_status+=[lines.split("\t")[8]]#+" "+lines.split("\t")[9]]
            #print(lines.split("\t")[6]+" "+lines.split("\t")[7])
#print(len(list_ref))
#print(len(list_alt))

data = pd.DataFrame(
    {'CHROM': list_chrom,
     'POS': list_POS,
     'REF': list_ref,
     'ALT': list_alt,
     'REF_READS': list_ref_reads,
     'ALT_READS': list_alt_reads,
     'VAF': list_VAF,
    })

dict_status={}
x=0
for item in list_status:
    if item not in dict_status:
        dict_status[item]=x
        x+=1

list_status2=[]
for item in list_status:
    list_status2+=[dict_status[item]]

data_truth = pd.DataFrame(
    {'STATUS': list_status2})

new_data=[]
list_ref2=[]
list_alt2=[]
#print(max([list_ref]))
for item in list_ref:
    list_ref2+=[item/max(list_ref)]
for item in list_alt:
    list_alt2+=[item/max(list_alt)]    
for x in range(len(list_chrom)):
    new_data+=[[list_chrom[x],list_POS[x],list_ref2[x],list_alt2[x],list_ref_reads[x],list_alt_reads[x],list_VAF[x]]]
#print(list_status2)
#for item in new_data:
    #for item2 in item:
        #if type(item2)==str:
        #    print(item2)
        #print(type(item2))
#print(data_truth.values.ravel())
#print(data)
#print(data_truth.values.ravel())
#print(data_truth)
skf = StratifiedKFold(n_splits=5)
skf.get_n_splits(new_data, list_status2)

target_names=[]
for item in list_status2:
    target_names+=[item]

for train_index, test_index in skf.split(new_data, list_status2):
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    clf3 = make_pipeline(StandardScaler(), SVC(gamma='auto'))
    print("SVC")
    clf3.fit(new_data[train_index], list_status2[train_index])
    #clf.predict(new_data_train2)
    print(clf3.score(new_data[test_index], list_status2[test_index]))
    #target_names = ['class 0', 'class 1', 'class 2']
    print(classification_report(list_status2[test_index], clf3.predict(new_data[test_index]), target_names=target_names))
    print("GaussianNB")
    NB = GaussianNB().fit(new_data[train_index], list_status2[train_index])
    print(NB.score(new_data[test_index], list_status2[test_index]))
    print(classification_report(list_status2[test_index], NB.predict(new_data[test_index]), target_names=target_names))
    print("LogisticRegression")
    clf = LogisticRegression(random_state=0).fit(new_data[train_index], list_status2[train_index])
    print(clf.score(new_data[test_index], list_status2[test_index]))
    print(classification_report(list_status2[test_index], clf.predict(new_data[test_index]), target_names=target_names))
    print("KNeighborsClassifier")
    neigh = KNeighborsClassifier(n_neighbors=9)
    neigh.fit(new_data[train_index], list_status2[train_index])
    print(neigh.score(new_data[test_index], list_status2[test_index]))
    print(classification_report(list_status2[test_index], neigh.predict(new_data[test_index]), target_names=target_names))
    clf2 = RandomForestClassifier(max_depth=2, random_state=0)
    print("RandomForestClassifier")
    clf2.fit(new_data[train_index], list_status2[train_index])
    print(clf2.score(new_data[test_index], list_status2[test_index]))
    print(classification_report(list_status2[test_index], clf2.predict(new_data[test_index]), target_names=target_names))
