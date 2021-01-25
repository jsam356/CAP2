#!/usr/bin/python
# coding: utf-8

# In[1]:


#So far the data is divided in an initial set and a series of additionla batches of data (so far one).

#ids will save the ids of the sample (SRAs), projects will save the bioproject identifiers (PJRNA). SRRs will save
#the SRR (run) identifier.
ids=[]
projects=[]
SRRs=[]

#The ids and corresponding project of the first batch is loaded first from Proj_UID. Proj_UID was generated
#in the ProjectMatch notebook. in the Match folder.

FirstBatch=open("Proj_UID.csv","r")

for line in FirstBatch:
    line=line.strip("\n").split(";")
    ids.append(line[0])
    projects.append(line[1])
    SRRs.append(line[2])


# In[2]:


#The second batch is loaded in the same way. Proj_UIDBatch2 was generated
#in the CheckforNew notebook in the Distances folder.
SecondBatch=open("Proj_UIDBatch2.csv","r")

for line in SecondBatch:
    line=line.strip("\n").split(";")
    
    if line[0] in ids:
        continue
    else:
        ids.append(line[0])
    projects.append(line[1])
    SRRs.append(line[2])


# In[3]:


#Now we load the Mash distances. Here they are loaded in the tabular output format of Mash Triangle 
#generated with the -E flag (if I'm not wrong).
import numpy as np

phy=open("NewDistances.tab","r")

#The distances are saved in a matrix. The dimensions of the matrix can be found by determining the length
#of ids.
distances=np.zeros((1617,1617))

#At each line the tab file is parsed and the distances added. The file includes, amongst others, the 
#id of all pairs of samples and their distance. We parse the ids and distance and place the distance in
#the corresponding spot in the distance matrix, so that the order of samples in the list ids corresponds 
#to the distances in the matrix.
for line in phy:
    line=line.strip("\n").split("\t")
    id1=line[0].split("_")[0].split("/")[1]
    id2=line[1].split("_")[0].split("/")[1]
    dist=float(line[2])
    
    distances[ids.index(id1),ids.index(id2)]=dist
    distances[ids.index(id2),ids.index(id1)]=dist
phy.close()


# In[4]:


#Hereon we're interested in a dereplicated list of all projects, so we dereplicate projects ans save as
#Ordprojects (for now).
Ordprojects=list(set(projects))
print(len(Ordprojects))


# In[5]:


#We generate col, a list which will map every sample to a number, and that number to a given project, based
#on the order of appareance of the projects in Proj_UID and Proj_UIDBatch2. More on col below.
col=[]

for i in projects:
    col.append(Ordprojects.index(i))


# In[6]:


#The original ordered project list, with the ith positions having the project of the ith id, is saved in
#UnOrdprojects, while the dereplicated list is now saved simply as projects.
UnOrdprojects=projects
projects=Ordprojects


# In[7]:


#Now we move on to load the keywords and MeSH terms for the bioprojects which have them. In general almost
#all bioprojects have keywords, while only a handful have MeSH Terms.

#For keywords the file Papers-Grid_Min.csv (custom made based on the Airtable data) contains the title of the
#project (or paper if there's a paper), the bioproject identifier and the keywords from airtable.
fpapers=open("Papers-Grid_Min.csv","r",encoding="utf-8")

#We'll save the keywords per project in keywords and the projects corresponding to each set of keywords in
#ProjforKeywords. Note that we're mapping projects to keywords, not samples to keywords. Also, since 
#Papers-Grid_Min.csv was downloaded from Airtable, the order of the projects is not necessarily the same as 
#in the list projects generated above (that's why we need another list).
ProjforKeywords=[]
keywords=[]

for line in fpapers:
#We use the PRJNA in the project id as a separator and get the list of keywords from every project.
    if "PRJNA" in line:
        line=line.strip("\n").split('PRJNA')
        projectID="PRJNA"+line[1].split(",")[0].split(";")[0].split(".pdf")[0].split(")")[0]

#Since we're only interested in the keywords of the projects whose samples have Mash distances, we filter
#the projects and get the keywords of those which are present in Proj_UID.csv or Proj_UIDBatch2.csv.
        if projectID in projects:
            ProjforKeywords.append(projectID)
            Prelkeywords=line[1].split(",")[1::]
            Realkeywords=[]
            for i in Prelkeywords:
                if i=="":
                    continue
                else:
#keywords are added in lowercase.
                    Realkeywords.append(i.lower())
            keywords.append(Realkeywords)
fpapers.close()


# In[8]:


#Now we define a dictionary which will store the projects (value) associated to each keyword (key) present
#in the dataset. 
kwordtoProj={}

for i in range(len(keywords)):
    for j in keywords[i]:
        if j in kwordtoProj.keys():
            kwordtoProj[j].append(ProjforKeywords[i])
        else:
             kwordtoProj[j]=[ProjforKeywords[i]]


# In[9]:


#For MeSH terms we employ the file MeSHTerms.csv which contains the major and minor MeSH terms per 
#article and project. This csv is partly made with notebook Title_Description_Extraction.csv in 
#Match/Keywords_MeSH, and partly done manually.

MeShTerms=open("MeSHTerms.csv","r",encoding="utf-8")

#We generate lists to hold the major and minor MeSH terms per bioproject (MeSHMaj and MeSHMin), and the 
#bioprojects which have them in the order they are found in MeSHTerms.csv
ProjforMeSH=[]
MesHMaj=[]
MesHMin=[]

#These two are temporary lists for holding the terms associated to single bioprojects.
thisMaj=[]
tempMin=[]

#Due to the way the csv with the MeSH terms is structured, a variable with the first project id in the file
#needs to be defined.
Currentproject="PRJNA450123"

#Now the meSH terms are loaded.
while True:
    line = MeShTerms.readline()
    if line is None or line=='':
        break
    else:
#The file has one line per Major-Minor MeSh term combination in each specific project, so there are many
#lines associated to each project.
        line=line.strip("\n").split(";")
        if line[1]==Currentproject:
            thisMaj.append(line[2])
            thisMin=[]
#While the same bioproject appears in the following line, we add Major and Minor terms into temporary lists.
            for i in range(3,len(line)):
                if line[i]!="":
                    thisMin.append(line[i])
            tempMin.append(thisMin)
#When the bioproject changes we save the terms of the prior bioproject in their corresponding lists.
        else:
            MesHMaj.append(thisMaj)
            MesHMin.append(tempMin)
            Currentproject=line[1]
            ProjforMeSH.append(line[1])
            thisMaj=[line[2]]
            tempMin=[]
            thisMin=[]
            for i in range(3,len(line)):
                if line[i]!="":
                    thisMin.append(line[i])
            tempMin.append(thisMin)
    
MeShTerms.close()


# In[10]:


#We generate dictionaries with the Major and Minor MeSH terms (key) and the bioprojects the're included
#in (value).

#Note that for minor MeSH terms we're considering the combination of the major and the minor term as a key,
#to avoid ambiguities.
MajMeSHtoProj={}
MinMeSHtoProj={}

for i in range(len(MesHMaj)):
    for j in range(len(MesHMaj[i])):
        if MesHMaj[i][j] in MajMeSHtoProj.keys():
            MajMeSHtoProj[MesHMaj[i][j]].append(ProjforMeSH[i])
        else:
             MajMeSHtoProj[MesHMaj[i][j]]=[ProjforMeSH[i]]
#Major and minor MeSH combinations are separated by _ to make them easily distinguishable.
        for z in MesHMin[i][j]:
            if MesHMaj[i][j]+"_"+z in MinMeSHtoProj.keys():
                MinMeSHtoProj[MesHMaj[i][j]+"_"+z].append(ProjforMeSH[i])
            else:
                MinMeSHtoProj[MesHMaj[i][j]+"_"+z]=[ProjforMeSH[i]]


# In[11]:


#Some of the projects we're working with have abstracts in Airtable which might contain the keywords which
#are associated to the projects. Therefore, it makes sense to look for keywords in the abstracts of the 
#projects, to make sure that all bioprojects which have a given keyword in their abstracts (but not in their
#original set of keywords) now have that keyword.

#Note that abstracts in Airtable might be shorts descriptions or paper abstracts per se, depending
#on the case.

#The abstracts are saved in Abstracts.csv, which contains the title of the project, the project ID, and the
#abstract itself.
abst=open("Abstracts.csv","r",encoding="utf-8")

while True:
    line = abst.readline()
    if line is None or line=='':
        break
#Due to the way the Airtable information must be downloaded, abstracts can span multiple lines in the file.
#We look for the PJRNA in the bioproject identifier and use that to identify the beggining of an abstract
#and save all its lines.
    else:
        if "PRJNA" in line:
            line=line.strip("\n").split("PRJNA")
            projectID="PRJNA"+line[1].split(",")[0].split(";")[0].split(".pdf")[0].split(")")[0]
            #print(projectID)
#Naturally we only consider projects whose samples have Mash distances.
            if projectID in projects:
                try:
#Some of the abstracts begin with double quotes, which are used to mark the beggining and end of them, so 
#they're used to parse them and mark the end of the abstract.
                    text=line[1].split('"')[-2]
                    while '"' not in text:
#When an abstract is found all its line are searched for the presence of any of the keywords. if any its find
#and its not already associated to the bioproject, the bioproject is added to the list of projects in the 
#keyword dictionary, and the keyword added to the list of keywords for that bioproject.
                        for j in kwordtoProj.keys():
                            if j in text.lower() and projectID not in kwordtoProj[j]:
                                kwordtoProj[j].append(projectID)
                                keywords[ProjforKeywords.index(projectID)].append(j)
                        text=abst.readline().strip("\n")
#For abstracts without double quotes we use commas.
                except IndexError:
                    text=line[1].split(',')[-1]
                    for j in kwordtoProj.keys():
                            if j in text.lower() and projectID not in kwordtoProj[j]:
                                kwordtoProj[j].append(projectID)
                                keywords[ProjforKeywords.index(projectID)].append(j)
abst.close()


# In[12]:


#Now we move on to consolidate keywords if necessary. Consolidate means group together keywords which directly 
#mean the same thing (Permanent or Primary Consolidation) and those which are related (Secondary Consolidation).

#Both consolidations are given by two csv files which are custom made: PermanentConsolidation.csv and 
#SecondaryConsolidation.csv. They contain in the first column the word in which the keywords will be consolidated 
#(consolidated term), and in the rest of the columns the keywords already present which will be consolidated 
#into the word in the first columns.

#We do the same procedure for both files. First we open them.
PC=open("PermanentConsolidation.csv","r")

#Then a dictionary is defined with the keys being the word to consolidate and the values lists of the present
#keywords to be consolidated into that word.
FirstConsolidate={}

for line in PC:
    line=line.strip("\n").split(";")
    
#The keywords present are appended in lower case to coincide with the ones in kwordtoProj.
    thisList=[]
    for i in line[1::]:
        if i!="":
            thisList.append(i.lower())
    FirstConsolidate[line[0]]=thisList

PC.close()

#The secondary consolidate file is parsed in the same way.
SC=open("SecondaryConsolidation.csv","r")

SecondaryConsolidate={}

for line in SC:
    line=line.strip("\n").split(";")
    
    thisList=[]
    for i in line[1::]:
        if i!="":
            thisList.append(i.lower())
#There are two terms which need to be handled uniquely due to they being consolidated terms for the first consolidation
#but being then consolidated into other term in the second consolidation.
        if "cabbage looper" in thisList:
            thisList[thisList.index("cabbage looper")]="Cabbage looper"
        if "16s" in thisList:
            thisList[thisList.index("16s")]="16S"
    SecondaryConsolidate[line[0]]=thisList

SC.close()


# In[13]:


#The function generate consolidation takes the dictionaries loaded from either PermanentConsolidation.csv or 
#SecondaryConsolidation.csv and updates kwordtoProj, keywords and ProjforKeywords so that the presence of the
#consolidated terms in the projects now shown, and the terms consolidated are erased from those objects.

#For the sake of generality, kwordtoProj, keywords and ProjforKeywords are not directly used in the function 
#definition, but they should be pased as the parameters wordDict,ProjWordList, and ProjtoWord.
def generateConsolidation(ConsolidateDict,wordDict,ProjWordList,ProjtoWord):
#The function first adds to wordDict the new consolidated terms as keys and the list of the bioprojects they're in
#as values. The consolidated terms are also added to the ordered list with keywords per bioproject. The terms to 
#be consolidated are added to the list toRemove.
    toRemove=[]
    for NewTerm in ConsolidateDict.keys():
        thisList=[]
        for OldTerm in ConsolidateDict[NewTerm]:
            thisList.extend(wordDict[OldTerm])
            toRemove.append(OldTerm)
        wordDict[NewTerm]=list(set(thisList))
#Here the consolidated terms are added to the ordered list of keywords per bioproject, and the terms to be consolidated
#are removed from it.
    for proj in ProjtoWord:
        for NewTerm in ConsolidateDict.keys():
            if proj in wordDict[NewTerm]:
                ProjWordList[ProjtoWord.index(proj)].append(NewTerm)
        for OldTerm in toRemove:
            if proj in wordDict[OldTerm]:
                while ProjWordList[ProjtoWord.index(proj)].count(OldTerm)!=0:
                    ProjWordList[ProjtoWord.index(proj)].remove(OldTerm)
#Finally, the terms to be consolidated are removed from the dictionary of keywords.
    for key in toRemove:
        try:
            del wordDict[key]
        except KeyError:
            pass


# In[14]:


#We run the first consolidation.
generateConsolidation(FirstConsolidate,kwordtoProj,keywords,ProjforKeywords)


# In[15]:


#Now we move to the process of generating coloring schemes for graphing. Earlier we defined col, an ordered list
#in which each sample is assigned a number which corresponds to the position of its corresponding bioproject in
#the bioprojects list. That will be useful for coloring per project in the graphs. However, we might want to
#color by the presence of specific terms.

#change_col generates coloring schemes based on the presence of keywords or MeSH terms. dictOne refers to a MeSh or
#keyword dictionary, words the terms to be used for the coloring shceme, and all_col whether we should consider
#simultaneous presence of different terms in single samples (all_col=1) or not (all_col=0).
def change_col(dictOne,words,all_col=0):
#First a number is assigned to every term. The range of the number can be changed by changing the value of
#z and the line where z number is added to z. The idea of first using 20 and adding 10 is to increase the spread
#of the values for having a more divergent coloring in matplotlib.

#colmap maps terms to numbers, and numMap numbers to terms.
    colmap={}
    numMap={}
    z=20
    for j in words:
        colmap[j]=z
        numMap[z]=j
        z+=10
#After the colors are assigned and saved in colmap and numMap then each sample is is checked (using cols to locate
#its project) to see if it harbors the term in question. 
    new_cols=[]
    for i in range(len(col)):
        proj=projects[col[i]]
#If all_col=1 then all terms are checked and the number associated to the sample is the average of the terms
#which its bioproject contains. The term itself will be a combination of all the terms present. 
        if all_col==1:
            this_col=[]
            thisTerm=""
            for j in words:
                if proj in dictOne[j]:
                    if thisTerm=="":
                        thisTerm=j
                    else:
                        thisTerm=thisTerm+" & "+j
                    this_col.append(colmap[j])
            if thisTerm=="":
#The samples with no terms are assigned 2 as an arbitrary value (this can be changed).
                new_cols.append(2)
                continue
            elif thisTerm not in colmap.keys():
                colmap[thisTerm]=[mean(this_col)]
                numMap[mean(this_col)]=thisTerm
            new_cols.append(mean(this_col))
#If all_col=0 then the sample is assigned the first term it has from the ones considered.
        else:
            found=0
            for j in words:
                if proj in dictOne[j]:
                    new_cols.append(colmap[j])
                    found=1
                    break
            if found==0:
                new_cols.append(2)
#The function returns the dictionaries mapping colors and terms and the ordered list of new numbers (colors)                
#for each sample.
    return new_cols,colmap,numMap


# #### Test for Significant Differences in Clustering for Keywords

# In[16]:


#Now we can do a more individualized and comprehensive test for each of the keywords. We'll consider each keyword 
#individually and perform three tests. 
import statistics 
import random
import copy

#One of the tests is going to be an aleatorization analysis. For that we'll define a function to perform an 
#aleatorizations. We'll need the distance matrix (dist), the number of times the aleatorization should take place
#(times) and the labels of the groups (labels).
def aleatorizationForTerms(dist,times,labels):
    differences=[]
    newlabels=copy.deepcopy(labels)
#What we do is shuffle the list of labels the number of times indicated by times.
    for i in range(0,times):
        random.shuffle(newlabels)
        group_distances=[]
        diff_distances=[]
#With the shuffling done we get the distances of samples with the term and without the term.
        for i in range(1,dist.shape[0]):
            for j in range(0,i):
                if newlabels[i]==20 and newlabels[j]==20:
                    group_distances.append(dist[i,j])
                elif newlabels[i]==2 and newlabels[j]==2:
                    diff_distances.append(dist[i,j])
#We compute each tieme the difference of the medians of the distance of both groups (with the term and without it) 
#based on the labels.
        differences.append(statistics.median(diff_distances)-statistics.median(group_distances))
#At the end the function returns a list with all the results of this median difference.
    differences.sort()
    return differences

# In[ ]:

output=open("SigKeywords.txt","w")
#Now we do the tests themselves. We'll do the same two tests done in the previous analyses (Wilcoxon and Kolmogorov-
#Smirnov) apart from the aleatorization. 
from scipy.stats import ranksums,kstest
significant=[]
abSignificant=[]
abundance=[]

#Running aleatorizations for all the keywords would be particularly costly computationally, so we establish an 
#abundance cut. This cut can be set based on the overall abundance of keywords per projects or samples. 
#A histogram of the abundance in bioprojects or samples of the keywords is shown below.
byProj=0

for term in kwordtoProj.keys():
#If we'll work with the abundance by bioproject then we set byProj=1 and estimate abundances based on this.
#If we set byProj=0 then we count abundances by sample and filter using samples. Note that regardless of this
#the distance analysis is always done using samples.

#The default cut is 5, but this can be changed here. The list abundances will save the abundance of all keywords
#in terms of samples or bioprojects.
    if byProj==1:
        abundance.append(len(kwordtoProj[term]))
        if len(kwordtoProj[term])<5:
            continue
        ncol,thisColMap,Wol=change_col(kwordtoProj,[term],0)
    else:
        ncol,thisColMap,Wol=change_col(kwordtoProj,[term],0)
        abundance.append(ncol.count(20))
        if abundance[-1]>25 or abundance[-1]<20:
            continue
    print(term)

#For each term we get the distances of samples which have the keyword or don't as we've done before.
    group_distances=[]
    diff_distances=[]
    for i in range(1,distances.shape[0]):
        for j in range(0,i):
            if ncol[i]==20  and ncol[j]==20:
                group_distances.append(distances[i,j])
            elif ncol[i]==2  and ncol[j]==2:
                diff_distances.append(distances[i,j])

#We then move on to carry out the three tests. First the Wilcoxon and Kolmogorov-Smirnov.
    a,b=ranksums(group_distances,diff_distances)
    c,d=kstest(group_distances,diff_distances)
#Then we get the real difference of medians between samples with the same term and different terms.
    actual_diff=statistics.median(diff_distances)-statistics.median(group_distances)
    #After that we get the result of the aleatorization and compare with the median distribution.
    Randdiff=aleatorizationForTerms(distances,1000,ncol)
    SignPosLow=int((len(Randdiff)*0.01))
    SignPosHigh=int((len(Randdiff)*0.99))
#We compare the p-values of the first two tests (with bonferroni correction) and check the aleatorization. For the
#aleatorization we check if the distances is smaller than the 1% percentile of the random median distance distribution.
    if b<0.01/(2*len(kwordtoProj)) and d<0.01/(2*len(kwordtoProj)) and (Randdiff[SignPosLow]>actual_diff or Randdiff[SignPosHigh]<actual_diff):
#We print the terms which pass the criteria with their abundance
        if Randdiff[SignPosLow]>actual_diff:
            where="low"
        else:
            where="high"
        if byProj==1:
            abSignificant.append(len(kwordtoProj[term]))
            print(term+"\t"+str(len(kwordtoProj[term]))+"\t"+where)
        else:
            abSignificant.append(abundance[-1])  
            print(term+"\t"+str(abundance[-1])+"\t"+where)
        significant.append(term)
        output.write(term+"\n")

output.close()