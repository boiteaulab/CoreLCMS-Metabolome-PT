# Python script for attributing molecular formula to a 'metabolome' file using a CoreMS featurelist library
# RMB Last updated  03/08/2024
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau

# Import modules
import os
from tempfile import tempdir
import warnings
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import time
#import scipy.stats
import scipy
from scipy import stats
import tracemalloc
from sklearn.cluster import KMeans
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


def metabolome_stats(metabolome,filter_a,filter_b):
    metabolome_a = metabolome.filter(like=filter_a)
    metabolome_b = metabolome.filter(like=filter_b)
    fchange=metabolome_a.mean(axis=1)/metabolome_b.mean(axis=1)
    fchange[fchange>2**8]=2**8
    fchange[fchange<(1/2**8)]=(1/2**8)
    metabolome['lf_change'] = np.log2(fchange)
    _, p = stats.ttest_ind(metabolome_a.transpose(), metabolome_b.transpose())
    metabolome['pvalue']=p
    metabolome=metabolome.sort_values(by='pvalue',ascending=True)

    metabolome['adj p-value']=metabolome['pvalue']*len(metabolome)/range(1,len(metabolome)+1)
    metabolome['Sig_Flag']=False
    #metabolome['Sig_flag'][(p<0.01) & (abs(metabolome['lf_change'])>1)]=True
    metabolome.loc[(metabolome['adj p-value']<0.05) & (abs(metabolome['lf_change'])>1),'Sig_Flag']=True
    metabolome=metabolome.sort_values(by='Alignment ID')

    #metabolome['adj pvalue'] = stats.false_discovery_control(p)  #Need Scipy v. 1.11 or later...

    return(metabolome)


def stoichiometric_classification(featurelist):
    
    assignedresults=featurelist
    assignedresults['C']=assignedresults['C'].fillna(0)
    assignedresults['H']=assignedresults['H'].fillna(0)
    assignedresults['O']=assignedresults['O'].fillna(0)
    assignedresults['N']=assignedresults['N'].fillna(0)
    assignedresults['P']=assignedresults['P'].fillna(0)

    # Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
    assignedresults['O/C']=assignedresults['O']/assignedresults['C']
    assignedresults['H/C']=assignedresults['H']/assignedresults['C']
    assignedresults['N/C']=assignedresults['N']/assignedresults['C']
    assignedresults['P/C']=assignedresults['P']/assignedresults['C']
    assignedresults['N/P']=assignedresults['N']/assignedresults['P']

    assignedresults['NOSC'] =  4 -(4*assignedresults['C'] + assignedresults['H'] - 3*assignedresults['N'] - 2*assignedresults['O'])/assignedresults['C']

    assignedresults['Stoichiometric Classification']='Unassigned'

    assignedresults.loc[(assignedresults['C']>0)
                        ,'Stoichiometric Classification'] = 'Other'
    

    assignedresults.loc[(assignedresults['O/C']<=0.6) & 
                        (assignedresults['H/C']>=1.32) & 
                        (assignedresults['N/C']<=0.126) &
                        (assignedresults['P/C']<0.35)
                        ,'Stoichiometric Classification'] = 'Lipid'

    assignedresults.loc[(assignedresults['O/C']<=0.6) & 
                        (assignedresults['H/C']>=1.32) & 
                        (assignedresults['N/C']<=0.126) &
                        (assignedresults['P/C']<0.35) &
                        (assignedresults['P']>0)
                        ,'Stoichiometric Classification'] = 'Phospholipid'

    '''
    assignedresults.loc[(assignedresults['O/C']>=0.61) & 
                        (assignedresults['H/C']>=1.45) & 
                        (assignedresults['N/C']>0.07) & 
                        (assignedresults['N/C']<=0.2) & 
                        (assignedresults['P/C']<0.3) & 
                        (assignedresults['O']>=3) &
                        (assignedresults['N']>=1)
                        ,'Stoichiometric Classification'] = 'Amino-Sugar'
    '''
    assignedresults.loc[(assignedresults['O/C']>=0.8) & 
                        (assignedresults['H/C']>=1.65) & 
                        (assignedresults['H/C']<2.7) &
                        (assignedresults['O']>=3) &
                        (assignedresults['N']==0)
                        ,'Stoichiometric Classification'] = 'Carbohydrate'

    assignedresults.loc[(assignedresults['O/C']>=0.5) & 
                        (assignedresults['O/C']<1.7) & 
                        (assignedresults['H/C']>1) & 
                        (assignedresults['H/C']<1.8) &
                        (assignedresults['N/C']>=0.2) & 
                        (assignedresults['N/C']<=0.5) & 
                        (assignedresults['N']>=2) &
                        (assignedresults['P']>=1) &
                        (assignedresults['S']==0) &
                        (assignedresults['Calculated m/z']>305) &
                        (assignedresults['Calculated m/z']<523)
                        ,'Stoichiometric Classification'] = 'Nucleotide'
    
    assignedresults.loc[(assignedresults['S']>0)
                        ,'Stoichiometric Classification'] = 'Organosulfur'
    
    assignedresults.loc[(assignedresults['O/C']<=1.15) & 
                        (assignedresults['H/C']<1.32) & 
                        (assignedresults['N/C']<0.126) &
                        (assignedresults['P/C']<=0.2) 
                        ,'Stoichiometric Classification'] = 'Phytochemical'



    assignedresults.loc[(assignedresults['O/C']>0.12) & 
                        (assignedresults['O/C']<=0.6) & 
                        (assignedresults['H/C']>0.9) & 
                        (assignedresults['H/C']<2.5) & 
                        (assignedresults['N/C']>=0.126) & 
                        (assignedresults['N/C']<=0.7) & 
                        (assignedresults['P/C']<0.17) & 
                        (assignedresults['N']>=1)
                        ,'Stoichiometric Classification'] = 'Peptide'

    assignedresults.loc[(assignedresults['O/C']>0.6) & 
                        (assignedresults['O/C']<=1) & 
                        (assignedresults['H/C']>1.2) & 
                        (assignedresults['H/C']<2.5) & 
                        (assignedresults['N/C']>=0.2) & 
                        (assignedresults['N/C']<=0.7) & 
                        (assignedresults['P/C']<0.17) & 
                        (assignedresults['N']>=1)
                        ,'Stoichiometric Classification'] = 'Peptide'

    assignedresults.loc[(assignedresults['Is Isotopologue']==1)
                        ,'Stoichiometric Classification'] = 'Isotopologue'

    return(assignedresults)

# Define CoreMS LCMS functions
def metabolome_annotation(featurelist,metabolome,rt,mz,offset):
    # Settings for particular file type
    rt='Average Rt(min)' #metabolome_file retention time column header
    mz='Average Mz' #metabolome_file m/z column header
    threshold=3 #Mass accuracy of metabolomic data (ppm). 

    elements=['C','H','O','N','P','S','Na']
    
    #metabolome[rt]=metabolome[rt]+offset

    timebins=featurelist.Time.unique()
    metabolite_annotations=[]
    for i in metabolome.iterrows():
        current=i[1].to_dict()
        ctime=current[rt]+offset
        cmass=current[mz]
        match=-1
        if(ctime>timebins.min()):
            match=timebins[timebins<ctime].max()
        annotations=featurelist[(featurelist['Time']==match) & ((abs(featurelist['m/z']-cmass)/cmass*1e6)<threshold)]
                
        current['all library matches']=len(annotations)
        current['annotated library matches']=len(annotations[annotations['C']>0])
        current['Annotation']='Undetected'
        current['Stoichiometric Classification']='Undetected'

        if(cmass<200):
            current['Annotation']='Below m/z range'

        if(cmass>900):
            current['Annotation']='Above m/z range'


        if len(annotations)>0:
            if len(annotations)>1:
                annotations=annotations[annotations['S/N']==max(annotations['S/N'])]

            current['Annotation']='Unassigned'
            if (annotations['C']>1).any():
                current['Annotation']='Assigned'    
            current['Stoichiometric Classification']=annotations['Stoichiometric Classification'].to_numpy()[0]
            current['Molecular Formula']=annotations['Molecular Formula'].to_numpy()[0]
            current['Theor m/z']=annotations['Calculated m/z'].to_numpy()[0]
            current['Library Time']=annotations['Time'].to_numpy()[0]
            current['Library m/z Error']=annotations['m/z Error (ppm)'].to_numpy()[0]
            current['Molecular Class']=annotations['Molecular Class'].to_numpy()[0]
            current['Library S/N']=annotations['S/N'].to_numpy()[0]
            current['Library Ion Charge']=annotations['Ion Charge'].to_numpy()[0]
            current['Library Confidence Score']=annotations['Confidence Score'].to_numpy()[0]
            current['Library Is Isotopologue']=annotations['Is Isotopologue'].to_numpy()[0]
            current['Metabolome m/z Error']=(cmass-annotations['Calculated m/z'].to_numpy()[0])/cmass*1e6
            current['O/C']=annotations['O/C'].to_numpy()[0]
            current['H/C']=annotations['H/C'].to_numpy()[0]
            current['N/C']=annotations['N/C'].to_numpy()[0]
            current['P/C']=annotations['P/C'].to_numpy()[0]
            current['DBE']=annotations['DBE'].to_numpy()[0]
            current['NOSC']=annotations['NOSC'].to_numpy()[0]

            for element in elements:
                current[element]=annotations[element].to_numpy()[0]

        metabolite_annotations.append(current)

    annotated_metabolome=pd.DataFrame(metabolite_annotations)
    return(annotated_metabolome)

if __name__ == '__main__':

    #### Change file settings here
    global data_dir
    data_dir='/CoreMS/usrdata/PT_metabolome_EDTA/'

    global featurelist_file
    #featurelist_file='Pt_featurelist.csv' 
    featurelist_file='20221101_LBA_Boiteau_Zorbax3p5_EDTApooled_T16_35.csv'

    global metabolome_file
    metabolome_file='PT_EDTA_metabolome.csv'
    
    rt='Average Rt(min)' #featurelist retention time column header
    mz='Average Mz' #featurelist m/z column header

    ##### End user input

    # starting the monitoring
    starttime = time.time()


    #Load metabolome and generate stats comparing groups
    metabolome=pd.read_csv(data_dir+metabolome_file)
    metabolome=metabolome_stats(metabolome,'EDTA1_t16_10nM','EDTA1_t16_1uM')

    #Load featurelist and generate stoichiometric classifications
    featurelist=pd.read_csv(data_dir+featurelist_file)
    featurelist=stoichiometric_classification(featurelist)
    
    #Annotate metabolome
    annotated_metabolome=metabolome_annotation(featurelist,metabolome,rt,mz,-1.43)
    annotated_metabolome.to_csv(data_dir+metabolome_file.replace('.csv','_annotated.csv'))


    print('Total execution time: %.2f min' %((time.time()-starttime)/60))
