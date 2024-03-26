# Python script for comparing internal standard peak area across samples. Outliers are flagged.
# RMB Last updated  2/07/2024
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau, Maria Christina Alvarez Rodriguez

##########
# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import warnings
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.rcParams['pdf.fonttype'] = 42


warnings.filterwarnings("ignore")
import sys
sys.path.append("./")

# Import required CoreMS modules
os.chdir('/CoreMS/')
from corems.mass_spectra.input import rawFileReader


def eic_plot(samplelist,filename):

    """
    Plots the total ion chromatogram (TIC) for each sample in the sample list.

    Args:
        samplelist (pandas.DataFrame): A DataFrame containing a 'File' column with file paths.
        filename (str): The filename to save the plot as.

    Returns:
        None
    """
        
    tics=[]
    for file in samplelist['File'].unique():
        parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
        tic=parser.get_tic(ms_type='MS')[0]
        tic_df=pd.DataFrame({'Time': tic.time,'Intensity': tic.tic,'Sample':file.replace('.raw','')})
        tics.append(tic_df)

    tics=pd.concat(tics)
    fig, (ax) = plt.subplots(1)
    plt.figure(figsize=(10,6))
    
    sns.lineplot(x='Time',y='Intensity',data=tics,ax=ax, hue='Sample')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Total Ion Current Intensity')
    ax.set_xlim(0,36)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    fig.savefig(filename,dpi=300,format='pdf')

def std_qc(samplelist,stdmass,std_timerange,filename):
        
    """
    Plots the extracted ion chromatogram (EIC) of the internal standard for each sample,
    calculates the peak area and retention time, flags outliers based on standard deviation,
    and saves the results to a CSV file and plots.

    Args:
        samplelist (pandas.DataFrame): A DataFrame containing a 'File' column with file paths.
        stdmass (float): The m/z value of the internal standard.
        std_timerange (list): A list containing the start and end time (in minutes) of the retention time range for the internal standard peak.
        filename (str): The filename to save the plot and results as.

    Returns:
        pandas.DataFrame: The sample list with additional columns for QC area, retention time, and QC pass/fail flag.
    """

    area={}
    rt={}

    fig, axs = plt.subplot_mosaic([['a','b']], figsize=(11,5), constrained_layout=True)
    axs['a'].set(xlabel='Time (min)',ylabel='Intensity',title='Internal Standard EIC = '+str(stdmass) + ' m/z')

    for file in samplelist['File'].unique():
        try:
            parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
            EIC=parser.get_eics(target_mzs=[stdmass],tic_data={},peak_detection=False,smooth=False)
            
            df=pd.DataFrame({'EIC':EIC[0][stdmass].eic,'time':EIC[0][stdmass].time})
            df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]
            area[file]=(sum(df_sub['EIC']))
            rt[file]=(df_sub.time[df_sub.EIC==df_sub.EIC.max()].max())
            axs['a'].plot(df_sub['time'],df_sub['EIC']/1e7,label=file[11:])
            print(file)
        except:
            print('No file found: ' + file)

    axs['a'].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs['a'].set_title('a', fontweight='bold', loc='left')
    axs['a'].set_ylabel('Intensity (x 1e7)')

    samplelist=samplelist.set_index('File')

    samplelist['qc_area']=pd.Series(area)
    samplelist['QC Retention time']=pd.Series(rt)

    # Flag outliers with peak area greater than 2x standard deviation of the mean 

    peak_stdv=samplelist.qc_area.std()
    peak_mean=samplelist.qc_area.mean()

    samplelist['qc_pass']=0
    for i in samplelist.index:
        if (abs(samplelist.qc_area[i]-peak_mean)<2*peak_stdv):
            samplelist.qc_pass[i]=1

    print(str(samplelist.qc_pass.sum()) + ' pass of ' + str(len(samplelist)))

    peak_stdv=samplelist[samplelist.qc_pass==1].qc_area.std()

    print(str(round(peak_stdv/peak_mean*100,1))+' % std dev')

    #Create plot of overlaid standard EICs
    sns.histplot(x='qc_area',data=samplelist,ax=axs['b'])
    axs['b'].set_xlabel('Internal Standard Peak Area')
    #axs['b'].set_xlim(0,20e7)
    axs['b'].set_title('b', fontweight='bold', loc='left')

    plt.savefig(filename,dpi=300,format='pdf')

    return(samplelist)


if __name__ == '__main__':

    # Set data directory and file names here
    global data_dir
    data_dir='/CoreMS/usrdata/PT_metabolome_EDTA/'

    global sample_list_name
    sample_list_name='Pt_metabolome_samples.csv' #Sample list must contain column with header 'File'


    # Read in sample list and load MS data (Either as existing sample list, or create sample list from .raw files in data_dir)
    samplelist=pd.read_csv(data_dir+sample_list_name)
    #samplelist=pd.DataFrame({'File':[f for f in os.listdir(data_dir) if '.raw' in f]})

    eic_plot(samplelist,data_dir+'pooled_tic.pdf')

    stdmass=678.2918 # m/z of cycanocobalamin
    std_timerange=[8,11] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_B12QCextended_EICplot.pdf'))
    samplelist2.to_csv(data_dir+sample_list_name.replace('.csv','_QC.csv'))
    
    sample_list_name='Pt_metabolome_samples_1uMFe.csv' #Sample list must contain column with header 'File'

    # Read in sample list and load MS data (Either as existing sample list, or create sample list from .raw files in data_dir)
    samplelist=pd.read_csv(data_dir+sample_list_name)

    stdmass=231.1015708 # m/z of feature
    rt=12.75
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_231_EICplot.pdf'))
    samplelist2.to_csv(data_dir+sample_list_name.replace('.csv','_QC.csv'))

    stdmass=219.10155 # m/z of feature
    rt=15.301
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_219_EICplot.pdf'))

    stdmass=249.11244 # m/z of feature
    rt=12.883
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_249_EICplot.pdf'))

    stdmass=253.12975 # m/z of feature
    rt=7.59
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_253_EICplot.pdf'))

    stdmass=267.0687 # m/z of feature
    rt=16.68
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_267_EICplot.pdf'))

    stdmass=267.0687 # m/z of feature
    rt=16.68
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_267_EICplot.pdf'))

    stdmass=269.2225 # m/z of feature
    rt=19.351
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_269_EICplot.pdf'))

    
    stdmass=277.07605 # m/z of feature
    rt=10.99
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_277_EICplot.pdf'))

    stdmass=281.16104# m/z of feature
    rt=9.063
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_281_EICplot.pdf'))

    stdmass=308.11795 # m/z of feature
    rt=8.287
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_308_EICplot.pdf'))


    stdmass=309.06552 # m/z of feature
    rt=9.794
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_309_EICplot.pdf'))

    stdmass=324.11267 # m/z of feature
    rt=7.837
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_324_EICplot.pdf'))

    
    stdmass=456.22778 # m/z of feature
    rt=12.912
    std_timerange=[rt-1.5,rt+1.5] # retention time range of peak (min)
    samplelist2=std_qc(samplelist,stdmass,std_timerange,data_dir+sample_list_name.replace('.csv','_456_EICplot.pdf'))