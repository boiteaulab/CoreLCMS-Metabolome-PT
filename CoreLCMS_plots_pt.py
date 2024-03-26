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


def rt_assign_plot(annotated_metabolome,filename,rt):

	#### Plot library assignments over time
    param='Stoichiometric Classification'

    assign_summary=[]

    for time in annotated_metabolome[rt].round().unique():
        current={}
        current['Time']=time
        for mol_class in class_order:
            current[mol_class]=len(annotated_metabolome[(annotated_metabolome[param]==mol_class) & (annotated_metabolome[rt].round()==time)])        
        assign_summary.append(current)

    #my_cmap = sns.color_palette("colorblind", as_cmap=True)
    #plt.style.use(my_cmap)

    df=pd.DataFrame(assign_summary)
    df=df.sort_values(by='Time')
    df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks',xlabel='Retention Time (min)')
    #sns.barplot(x='Retention Time (min)',y='Peaks',hue='Stoichiometric Classification',hue_order=class_order,data=df)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    plt.savefig(filename, bbox_inches='tight',format='pdf')

def metabolome_assign_plot(annotated_metabolome,filename,rt,mz):
    #### Plot and save annotated metabolome stochiometric classifications (rt vs 'O/C' and 'm/z')
    param='Molecular Class'
    #Set up a two panel plot
    fig, ((ax1, ax2)) = plt.subplots(1,2)
    fig.set_size_inches(12, 6)
    '''
    annotated_metabolome[param][annotated_metabolome['Na']>0]='Na Adduct'
    annotated_metabolome[param][annotated_metabolome[param]=='Unassigned']='Other'

    annotated_metabolome[param][annotated_metabolome['Library Is Isotopologue']==1]='Other'
    #annotated_metabolome[param][annotated_metabolome['P']>0]='CHOP'
    horder=['CHO','CHON','CHOS','CHOP','CHONS','CHONP','CHOSP','CHONSP','Na Adduct','Other']
    '''

    #Panel A
    #sns.scatterplot(x=rt,y=mz,hue=param,hue_order=horder,data=annotated_metabolome,ax=ax1, edgecolor='none')
    sns.scatterplot(x=rt,y=mz,hue='Stoichiometric Classification',hue_order=class_order,data=annotated_metabolome,ax=ax1, edgecolor='none')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax1.set_title('(a)', fontweight='bold', loc='left')
    ax1.set(xlabel='Retention time (min)',ylabel='m/z')

    #Panel B
    #sns.scatterplot(x='Metabolome m/z Error',y='Theor m/z',data=annotated_metabolome,hue=param,hue_order=horder,ax=ax2,legend=False, edgecolor='none')
    sns.scatterplot(x='O/C',y='H/C',hue='Stoichiometric Classification',hue_order=class_order,data=annotated_metabolome,ax=ax2, edgecolor='none',legend=False)

    ax2.set_title('(b)', fontweight='bold', loc='left')
    ax2.set(xlabel='O/C',ylabel='H/C')

    fig.tight_layout()

    fig.savefig(filename, dpi=300,format='pdf')
    '''
    sns.jointplot(x='Metabolome m/z Error',y='Theor m/z',data=annotated_metabolome,hue=param,hue_order=horder,ax=ax2,legend=False, edgecolor='none')

    plt.savefig(filename.replace('.pdf','_joint.pdf'), dpi=300,format='pdf')
    '''

def metabolome_error_plot(annotated_metabolome,filename,rt,mz):
    #### Plot and save annotated metabolome stochiometric classifications (rt vs 'O/C' and 'm/z')
  
    me_compare=annotated_metabolome[['Theor m/z','Metabolome m/z Error']].dropna()
    me_compare['Instrument']='Orbitrap MS'
    me_compare2=annotated_metabolome[['Theor m/z','Library m/z Error']].rename(columns={'Library m/z Error': 'Metabolome m/z Error'}).dropna()
    me_compare2['Instrument']='21T FT-ICR MS'
    
    g=sns.jointplot(x='Metabolome m/z Error',y='Theor m/z',data=pd.concat([me_compare,me_compare2]),hue='Instrument')
    g.set_axis_labels('m/z Error (ppm)', 'Theoretical m/z')
    #g.ax_marg_y.remove()

    g.savefig(filename, dpi=300,format='pdf')

def mz_bar_plot(annotated_metabolome,filename,rt,mz):
    #### Plot and save annotated metabolome stochiometric classifications (rt vs 'O/C' and 'm/z')
    
    #param='Annotation'
    assign_summary=[]
    ms=[100,200,300,400,500,600,700,800,900]

    for m in ms:
        current=annotated_metabolome[(annotated_metabolome[mz]>m) & (annotated_metabolome[mz]<m+100)]
        total=len(current)
        msms_collected=len(current[current['MS/MS assigned']==True])
        msms_assigned=len(current[current['MS/MS matched']==True])
        library_detected=len(current[current['all library matches']>0]) 
        library_assigned=len(current[current['Theor m/z']>0]) 
        assign_summary.append({'m/z':str(m)+'-'+str(m+100),'Total':total,'MS/MS collected':msms_collected,'MS/MS assigned':msms_assigned,'Library detected':library_detected,'Library annotated':library_assigned})

    df=pd.DataFrame(assign_summary)
    df=df.sort_values(by='m/z')
    df.plot.bar(x='m/z',y=df.columns[1:],stacked=False,ylabel='# of Features',xlabel='m/z range')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    plt.savefig(filename, bbox_inches='tight',format='pdf')

def volcano_plot(metabolome,filename):

    fig, axs = plt.subplot_mosaic([['a','b'],['a','c']], figsize=(10,6), constrained_layout=True)
    fig.set_size_inches(12, 6)
    
    metabolome['-log p']=-np.log10(metabolome['adj p-value'])

    sig_metabolome=metabolome[(metabolome['adj p-value']<0.05) & (abs(metabolome['lf_change'])>1)]
    
    #Panel A
    sns.scatterplot(x='lf_change',y='-log p',color='lightgray',data=metabolome,ax=axs['a'], edgecolor='none')
    sns.scatterplot(x='lf_change',y='-log p',hue='Stoichiometric Classification',hue_order=class_order,data=sig_metabolome,ax=axs['a'], edgecolor='none')
    axs['a'].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    axs['a'].set_title('(a)', fontweight='bold', loc='left')
    axs['a'].set(xlabel='log2 fold change (Fe deficient / Fe replete)',ylabel='-log10 adjusted p-value')

    #Panel B
    sub_metabolome=sig_metabolome[sig_metabolome['lf_change']<-1]
    assign_summary=[]
    for c in class_order:
        current=sub_metabolome[sub_metabolome['Stoichiometric Classification']==c]
        assign_summary.append({'Stoichiometric Classification':c,'Features':len(current)})
    df=pd.DataFrame(assign_summary)
    df.plot.bar(x='Stoichiometric Classification',y='Features',ax=axs['b'],legend=None)
    #axs['b'].set_xticklabels(axs['b'].get_xticklabels(),rotation=0)
    axs['b'].set_title('(b) More abundant in Fe replete treatment', fontweight='bold', loc='left')
    #axs['a'].set_ylim(0,30000)
    axs['b'].set(xlabel='Stoichiometric Classification',ylabel='# of Features')

    #Panel C
    sub_metabolome=sig_metabolome[sig_metabolome['lf_change']>1]
    assign_summary=[]
    for c in class_order:
        current=sub_metabolome[sub_metabolome['Stoichiometric Classification']==c]
        assign_summary.append({'Stoichiometric Classification':c,'Features':len(current)})
    df=pd.DataFrame(assign_summary)

    df.plot.bar(x='Stoichiometric Classification',y='Features',ax=axs['c'],legend=None)
    #axs['c'].set_xticklabels(axs['c'].get_xticklabels(),rotation=0)
    axs['c'].set_title('(c) More abundant in Fe deficient treatment', fontweight='bold', loc='left')
    #axs['a'].set_ylim(0,30000)
    axs['c'].set(xlabel='Stoichiometric Classification',ylabel='# of Features')

    #sns.scatterplot(x='O/C',y='H/C',data=metabolome[(metabolome['pvalue']<0.01) & (metabolome['lf_change']<-1)],color='red',ax=ax2, edgecolor='none')
    #sns.scatterplot(x='O/C',y='H/C',data=metabolome[(metabolome['pvalue']<0.01) & (metabolome['lf_change']>1)],color='green',ax=ax2, edgecolor='none')
    #ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    #fig.tight_layout()

    fig.savefig(filename, dpi=300,format='pdf')


    #metabolome['adj pvalue'] = stats.false_discovery_control(p)  #Need Scipy v. 1.11 or later...


def produced_metabolite_barchart(metabolome,filename):

    #fig, ((ax1)) = plt.subplots(1,2)
    #fig.set_size_inches(10, 6)

    assign_summary=[]
    for c in class_order:
        current=metabolome[metabolome['Stoichiometric Classification']==c]
        current_both=len(current[(current['1uM T16']/current['1uM T1']>2) & (current['10nM T16']/current['10nM T1']>2)])
        current_replete=len(current[(current['1uM T16']/current['1uM T1']>2)])-current_both
        current_deplete=len(current[(current['10nM T16']/current['10nM T1']>2)])-current_both

        assign_summary.append({'Stoichiometric Classification':c,'Both':current_both,'Fe Replete Only':current_replete,'Fe Deficient Only':current_deplete})
    df=pd.DataFrame(assign_summary)

    df.plot.bar(x='Stoichiometric Classification',y=df.columns[1:],stacked=True,xlabel='Stoichiometric Classification',ylabel='# of Features')
    #df.plot.bar(x='Stoichiometric Classification',y=df.columns[1:],stacked=True,xlabel='Stoichiometric Classification',ylabel='# of Features')
    #axs['b'].set_xticklabels(axs['b'].get_xticklabels(),rotation=0)
    #ax1.set_title('(Fe replete metabolites', fontweight='bold', loc='left')
    #axs['a'].set_ylim(0,30000)
    #ax1.set(xlabel='Stoichiometric Classification',ylabel='# of Features')
    plt.tight_layout()

    plt.savefig(filename, dpi=300,format='pdf')



def Fe_cond_comparison_scatter2(metabolome,filename):

    fig, ((ax1,ax2)) = plt.subplots(1,2)
    fig.set_size_inches(10, 4)

    current=metabolome[(metabolome['1uM T16']/metabolome['1uM T1']>2) | (metabolome['10nM T16']/metabolome['10nM T1']>2)]

    current['condition']='Both'
    current.loc[(current['1uM T16']/current['1uM T1']<2),'condition']='Fe deficient only'
    current.loc[(current['10nM T16']/current['10nM T1']<2),'condition']='Fe replete only'
    current=current[current['condition']!='Both']
    #Panel A
    sns.scatterplot(x='O/C',y='H/C',hue='condition',data=current,ax=ax1, edgecolor='none',legend=False)
    #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax1.set_title('(a)', fontweight='bold', loc='left')
    #ax1.set(xlabel='Retention time (min)',ylabel='m/z')
    
    #Panel B
    sns.scatterplot(x='N/C',y='H/C',hue='condition',data=current,ax=ax2, edgecolor='none')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax2.set_title('(b)', fontweight='bold', loc='left')

    plt.tight_layout()

    plt.savefig(filename, dpi=300,format='pdf')


def Fe_cond_comparison_scatter(metabolome,filename):


    current=metabolome[(metabolome['1uM T16']/metabolome['1uM T1']>2) | (metabolome['10nM T16']/metabolome['10nM T1']>2)]
    #current=current[current['Stoichiometric Classification']=='Lipid']
    current.loc[(current['1uM T16']/current['1uM T1']>2),'condition']='Fe replete only'
    current.loc[(current['10nM T16']/current['10nM T1']>2),'condition']='Fe deficient only'
    current.loc[(current['1uM T16']/current['1uM T1']>2)&(current['10nM T16']/current['10nM T1']>2),'condition']='Both'

    #current=current[current['condition']!='Both']

    sns.jointplot(x='O/C',y='H/C',hue='condition',data=current, edgecolor='none', marginal_kws={'common_norm':False})

    plt.tight_layout()

    plt.savefig(filename, dpi=300,format='pdf')


def Fe_cond_comparison_scatter3(metabolome,filename):

    fig, ((ax1,ax2)) = plt.subplots(1,2)
    fig.set_size_inches(10, 4)

    current=metabolome[(metabolome['1uM T16']/metabolome['1uM T1']>2) | (metabolome['10nM T16']/metabolome['10nM T1']>2)]

    current['condition']='Both'
    current.loc[(current['1uM T16']/current['1uM T1']<2),'condition']='Fe deficient only'
    current.loc[(current['10nM T16']/current['10nM T1']<2),'condition']='Fe replete only'
    #Panel A
    sns.violinplot(x='H/C',y='condition', data=current,ax=ax1, edgecolor='none',legend=False)
    #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax1.set_title('(a)', fontweight='bold', loc='left')
    #ax1.set(xlabel='Retention time (min)',ylabel='m/z')
    
    #Panel B
    sns.violinplot(x='N/C',y='condition', data=current,ax=ax2, edgecolor='none',legend=False)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax2.set_title('(b)', fontweight='bold', loc='left')

    plt.tight_layout()

    plt.savefig(filename, dpi=300,format='pdf')

def class_bar_plot(annotated_metabolome,filename,rt,mz):
    #### Plot and save annotated metabolome stochiometric classifications (rt vs 'O/C' and 'm/z')
    
    param='Stoichiometric Classification'
    assign_summary=[]
    LCMS_annotations=annotated_metabolome

    LCMS_annotations=LCMS_annotations[LCMS_annotations['cluster']!='None']
    LCMS_annotations=LCMS_annotations[LCMS_annotations['Stoichiometric Classification']!='Isotopologue']

    for c in LCMS_annotations['cluster'].unique():
        current=LCMS_annotations[LCMS_annotations['cluster']==c]
        d={}
        d['cluster']=c
        for m in LCMS_annotations[param].unique():
            d[m]=len(current[current[param]==m])
        assign_summary.append(d)
    df=pd.DataFrame(assign_summary)
    
    plt.figure(figsize=(3,2))
    df.plot.bar(x='cluster',y=df.columns[1:],stacked=False,xlabel='Treatment',ylabel='Number of metabolites most abundant in treatment',legend=True,fontsize=14)
    plt.savefig(filename+'cluster_bars.pdf', bbox_inches='tight',format='pdf')

### Make a heat map w/ heirarchical clustering
def metabolome_cluster_map(metabolome,clustermethod,filterstring,savefile):
    abundances=metabolome.filter(regex=filterstring)
    abundances=abundances[abundances.max(axis=1)>1]
    norm_abundances=abundances.div(abundances.max(axis=1),axis=0)
    norm_abundances=norm_abundances.transpose()

    h=sns.clustermap(norm_abundances,row_cluster=True,cmap='mako',method=clustermethod)
    h.savefig(savefile,dpi=300,format='pdf')

def feature_summary(metabolome):
    summary={}
    summary['Total']=len(metabolome)
    summary['MS/MS assigned']=len(metabolome[metabolome['MS/MS matched']==True])
    summary['Library assigned']=len(metabolome[metabolome['C']>0])
    summary['Outside library mass range']=len(metabolome[metabolome['Annotation']=='Below m/z range'])+len(metabolome[metabolome['Annotation']=='Above m/z range'])
    summary['Inside range but not detected in library']=len(metabolome[metabolome['Annotation']=='Undetected'])
    summary['Detected in library but not assigned']=len(metabolome[metabolome['Annotation']=='Unassigned'])
    return(summary)
    #Details of global assignments


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


    #Annotate metabolome
    annotated_metabolome=pd.read_csv(data_dir+metabolome_file.replace('.csv','_annotated.csv'))

    #class_order=['Undetected','Unassigned','Peptide','Lipid','Organosulfur', 'A-Sugars', 'Phospholipid', 'Phytochemical', 'Carbohydrate', 'Other', 'Isotopologue']
    class_order=['Undetected','Unassigned','Peptide','Lipid','Organosulfur', 'Phospholipid', 'Phytochemical', 'Carbohydrate', 'Other', 'Isotopologue']

    #Get subset of 'blank subtracted' and 'significant' features between Fe deficient vs Fe replete treatments
    annotated_metabolome_blanksubtract=annotated_metabolome[annotated_metabolome['Blank Flag']<0.5]
    annotated_metabolome_sig=annotated_metabolome_blanksubtract[annotated_metabolome_blanksubtract['Sig_Flag']==True]
    
    Table1={}
    Table1['All metabolome features']=feature_summary(annotated_metabolome)
    Table1['Blank subtracted features']=feature_summary( annotated_metabolome_blanksubtract)
    Table1['Fe replete features']=feature_summary(annotated_metabolome_blanksubtract[annotated_metabolome_blanksubtract['1uM T16']/annotated_metabolome_blanksubtract['1uM T1']>2])
    Table1['Fe deficient features']=feature_summary(annotated_metabolome_blanksubtract[annotated_metabolome_blanksubtract['10nM T16']/annotated_metabolome_blanksubtract['10nM T1']>2])
    Table1['Overlap Fe features']=feature_summary(annotated_metabolome_blanksubtract[(annotated_metabolome_blanksubtract['1uM T16']/
                                                  annotated_metabolome_blanksubtract['1uM T1']>2) & 
                                                  (annotated_metabolome_blanksubtract['10nM T16']/
                                                   annotated_metabolome_blanksubtract['10nM T1']>2)])
    Table1['Significant features']=feature_summary( annotated_metabolome_sig)
    pd.DataFrame(Table1).to_csv(data_dir+'Table1_Summary.csv')

    Table2=annotated_metabolome_sig[(annotated_metabolome_sig['lf_change']>1) & (annotated_metabolome_sig['C']>1)]
    Table2=Table2.sort_values(by=mz)
    Table2.to_csv(data_dir+'Table2_FeSig.csv')

    annotated_metabolome=annotated_metabolome.sort_values(by='C')

    ### Generate Plots

    sns.set_palette(None)

    metabolome_error_plot(annotated_metabolome,data_dir+'Figure3_assignment_error_plot.pdf',rt,mz)

    rt_assign_plot(annotated_metabolome,data_dir+'Figure4_assignment_bargraph.pdf',rt)

    metabolome_assign_plot(annotated_metabolome,data_dir+'FigureS2_assignment_scatterplots.pdf',rt,mz)

    mz_bar_plot(annotated_metabolome,data_dir+'Figure5_mzbar_plot.pdf',rt,mz)

    metabolome_cluster_map(annotated_metabolome_sig,'ward','EDTA1_',data_dir+'Figure6_heatmap.pdf')

    volcano_plot(annotated_metabolome[annotated_metabolome['Blank Flag']<0.5],data_dir+'Figure7_volcano_plot.pdf')

    produced_metabolite_barchart(annotated_metabolome_sig,data_dir+'FigureS4_produced_metabolite_barplot.pdf')

    Fe_cond_comparison_scatter(annotated_metabolome_blanksubtract,data_dir+'FigureS5_Fe_cond_scatterplot.pdf')
    Fe_cond_comparison_scatter3(annotated_metabolome_blanksubtract,data_dir+'FigureS5_Fe_cond_scatterplot3.pdf')

    print('Total metabolome features: ')
    print(len(annotated_metabolome))
    print('Annotated metabolome features: ')
    print(len(annotated_metabolome[annotated_metabolome['C']>0]))

    print('Total execution time: %.2f min' %((time.time()-starttime)/60))
