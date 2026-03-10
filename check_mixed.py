import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from datetime import datetime
import argparse


# check extrema in allele freq distributions per timepoint for a given patient
def check_extrema_pat(pat,threshold=300,width=0.05,min=0.1): 
    def check_extrema_date(data,date,threshold=300,width=0.05,min=0.1):
        current_df=data[data['date']==date]
        x = np.linspace(0, 1, 500)
        kde = gaussian_kde(current_df['allele_freq'])
        y = kde(x)
    
        # Find local minima and maxima
        maxima = argrelextrema(y, np.greater)[0]
        max_x=[]
        peak_strengths=[]
        for m in x[maxima]:
            count=((current_df['allele_freq']>(m-width)) & (current_df['allele_freq']<(m+width))).sum()
            if count>threshold and m>min:
                max_x.append(m)
                peak_strengths.append(count)
        return max_x, peak_strengths

    df_pat=load_pat(pat)

    num_peaks=[]
    for date in df_pat['date'].unique():
        mx, ps=check_extrema_date(df_pat,date,threshold=threshold,width=width,min=min)
        num_peaks.append(len(mx))
    
    #plt.hist(num_peaks, bins=range(0, max(num_peaks)+2), align='left')
    #plt.xlabel('Number of peaks per timepoint')
    #plt.ylabel('Frequency')

    return num_peaks

# load patient data, filter fixed and 0 freq snp
def load_pat(pat):
    def is_valid_date(date_str): # check if date is valid
        try:
            datetime.strptime(date_str, "%d-%m-%Y")
            return True
        except ValueError:
            return False
    # load df
    df_pat=pd.read_csv(f'{pat}/mut_dyn/mut_evolution_dates.txt',sep='\t')
    df_pat=df_pat[df_pat['date'].map(is_valid_date)] # filter invalid dates
    df_pat['date']=pd.to_datetime(df_pat['date'],format='%d-%m-%Y')
    df_pat.sort_values('date',inplace=True)

    # get series of fixed snp. Index pos and value True if fixed
    df_gcalls=df_pat[df_pat['passed_filter']==True]
    fixed_series=df_gcalls.groupby('position').min()['allele_freq']>=0.9
    # filter dataframe to remove fixed snp
    fixed_positions=fixed_series[fixed_series].index
    df_pat=df_pat[~df_pat['position'].isin(fixed_positions)]

    # remove snp were at 0 freq. If snp only at few time points might be error or contamination. 
    df_pat=df_pat[df_pat['allele_freq']>0]

    return df_pat

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Check mixed infection by looking at allele freq distributions per timepoint and counting extrema above a threshold")
    parser.add_argument('-p', '--patient', type=str, required=True,
                        help='patient name')
    parser.add_argument('-k', '--peak', type=int, required=True,
                        help='peak strength threshold')
    parser.add_argument('-w', '--width', type=float, required=True,
                        help='width around peak to count snp')
    parser.add_argument('-n','--peak_num', type=int, required=False, default=2,
                        help='number of peaks  to call mixed sample')
    parser.add_argument('-s', '--samples', type=int, required=True,
                        help='number of mixed samples to call patient mixed')
    parser.add_argument('--min', type=float, required=False, default=0.1,
                        help='minimum allele freq for peak to be considered')
    
    args = parser.parse_args()
    patient=args.patient
    peak_strength=args.peak
    width=args.width
    peak_num=args.peak_num
    samples=args.samples
    min_af=args.min

    # get number of peaks per timepoint passing threshold
    num_peaks=check_extrema_pat(patient,threshold=peak_strength,width=width,min=min_af)
    mixed_timepoints=sum(n>=peak_num for n in num_peaks)
    if mixed_timepoints>=samples:
        mixed=True
    else:
        mixed=False
    print(f'{patient}\t{mixed_timepoints}\t{mixed}')
