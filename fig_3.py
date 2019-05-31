import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy import signal
import random

import pickle
import os

from mfm import MFM

def fig_3():
    def get_data():
        """
        Loads data from file, or generates data if no file found
        """
        #First, look for saved data
        path = 'data/fig_3.npz'
        if os.path.isfile(path):
            print('Loading data from cache')
            data_dict = pickle.load(open(path,'rb'))

            data = data_dict['data']
            t    = data_dict['t']
            dt   = data_dict['dt']

        #If no data, generate data
        else:
            print('Generating data for Figure 3')
            conditions = [{'DD':False},
                          {'DD':True},
                          {'DD':True,'cDBS':True,'cDBS_amp':4.13}]

            for c in conditions:
                c['tstop'] = 100

            mfms = []
            data = []
            for i in range(3):
                mfms.append(MFM(**conditions[i]))
                mfms[i].run()

                dt = mfms[i].params['dt']
                time_series = mfms[i].S[:,mfms[i].struct['p2']]
                time_series = np.split(time_series,5)[-1] #get last 5th
                time_series -= np.mean(time_series)
                t = np.arange(len(time_series))
                t = t*dt

                data.append(time_series)

            data_dict = {
                'data' : data,
                't'    : t,
                'dt'   : dt
            }

            if not os.path.isdir(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path))
            pickle.dump(data_dict,open(path,'wb'))

        return data_dict, data, t, dt

    def plot(data, t, dt):
        # Figures
        # Timeseries Figure
        #-----------------------------------------------------------------------
        # Create figure and subplots
        fig = plt.figure(figsize=(10,4))

        gs = gridspec.GridSpec(1,2)

        gs0 = gridspec.GridSpecFromSubplotSpec(3,1, subplot_spec=gs[0])
        ax0 = [plt.subplot(gs0[0,0])]
        ax0.append(plt.subplot(gs0[1,0],sharey=ax0[0]))
        ax0.append(plt.subplot(gs0[2,0],sharey=ax0[0]))

        gs1 = gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec=gs[1],hspace=0.5)
        ax1 = [plt.subplot(gs1[0,0])]
        ax1.append(plt.subplot(gs1[1,0],sharex=ax1[0],sharey=ax1[0]))

        c = 'k'
        ax0[2].set_xlabel('time (s)')
        ax0[1].set_ylabel('LFP (mV)')

        ts_legends=['naive','DD','cDBS']
        for i in range(len(data)):
            ax0[i].plot(t,data[i],label=ts_legends[i],c=c)
            leg = ax0[i].legend(bbox_to_anchor=(1.02,1.3),loc='upper right', handlelength=0, handletextpad=0, borderpad=0, frameon=False)
            for item in leg.legendHandles: item.set_visible(False)

            ax0[i].set_xlim((10,10.5))
            ax0[i].set_ylim((-2.5,2.5))
            ax0[i].set_yticks([-2,0,2])

        # Remove spines
        for axes in ax0[:-1]:
            axes.spines['bottom'].set_visible(False)
            axes.xaxis.set_ticks([])

        # PSD Figure
        #------------------------------------------------------------------------
        ax1[1].set_xlabel('frequency (Hz)')

        for i in range(len(data)):
            fmax = 100    
            f,Pxx_den = signal.welch(data[i],1/dt,nperseg=4096)
            Pxx_den = 10*np.log10(Pxx_den**2)
            if i == 0:
                ax1[0].plot(f[f<fmax],Pxx_den[f<fmax],color='C0')
            elif i == 1:
                ax1[0].plot(f[f<fmax],Pxx_den[f<fmax],color='C3')
                ax1[1].plot(f[f<fmax],Pxx_den[f<fmax],color='C3')
            else:
                ax1[1].plot(f[f<fmax],Pxx_den[f<fmax],color='C2')


        ax1[0].set_xlim((0,100))

        x = ax1[0].figbox.bounds[0]
        y_upper = ax1[0].figbox.bounds[1] + ax1[0].figbox.bounds[3]
        y_lower = ax1[1].figbox.bounds[1]
        y_center = (y_upper - y_lower) / 2

        ax1[0].set_ylabel('power (dB/Hz)     ',ha='right')

        gs.tight_layout(fig,w_pad=2)

        # Add legends
        ax1[0].legend(['naive','DD'], frameon=False, borderpad=0, bbox_to_anchor=(1.05,1.05))
        ax1[1].legend(['DD','cDBS'],  frameon=False, borderpad=0, bbox_to_anchor=(1.05,1.05))

        # Add subfigure labels (a, b)
        fontdict={'size': 'large',
                  'weight' : 'bold'}
        fig.text(0.0, 1, 'a', fontdict=fontdict, verticalalignment='top')
        fig.text(0.5, 1, 'b', fontdict=fontdict, verticalalignment='top')
                
            
    data_dict, data, t, dt = get_data()
    
    plot(data, t, dt)

def main():
    random.seed(0)
    fig_3()
    plt.show()

if __name__ == '__main__':
    main()
