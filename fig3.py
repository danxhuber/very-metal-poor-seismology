import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from echelle import interact_echelle 
from echelle import plot_echelle
from echelle import echelle as calc_echelle
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel

plt.rcParams['xtick.major.size'] = 9
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 6
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 9
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 6
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size']=18
plt.rcParams['mathtext.default']='regular'
plt.rcParams['lines.markersize']=8
plt.rcParams['xtick.major.pad']='3'
plt.rcParams['ytick.major.pad']='3'
plt.rcParams['ytick.minor.visible'] = 'True'
plt.rcParams['xtick.minor.visible'] = 'True'
plt.rcParams['xtick.direction'] = 'inout'
plt.rcParams['ytick.direction'] = 'inout'
plt.rcParams['ytick.right'] = 'True'
plt.rcParams['xtick.top'] = 'True'

colors=['#FF9408','#DC4D01','#00A9E0','#016795']

freqs=ascii.read('data/freqs.csv')
freqs_mod=ascii.read('data/Freqs_best_fit.dat')

dat=ascii.read('data/KIC8144907_ps.txt')
freq=dat['col1']
amp=dat['col2']
dnu=17.47

gauss_kernel = Gaussian1DKernel(2)
pssm = convolve(amp, gauss_kernel)
echx, echy, echz=calc_echelle(freq,pssm, dnu,fmin=125, fmax=283)

echz2=np.zeros((echz.shape[0],echz.shape[1]*2))
echz2[:,0:echz.shape[1]]=echz
echz2[:,echz.shape[1]:echz.shape[1]*2]=echz

echx2=np.concatenate((echx,echx+np.max(echx))) 

plt.ion()
plt.clf()

plt.subplot(1,3,1)
plt.plot(freq,pssm,color='black')
plt.xlim([160,260])
plt.ylim([0,7000])
plt.xlabel('Frequency ($\mu$Hz)')
plt.ylabel('Power Density (ppm$^{2}$ $\mu$Hz$^{-1}$)')
plt.subplot(1,3,2)
alf=0.5

plt.imshow(np.log10(echz),aspect="auto",extent=(echx.min(), echx.max(), echy.min(), echy.max()),origin="lower",cmap='Greys',interpolation='None')

um=np.where(freqs['l'] == 0)[0]
plt.plot(freqs['freq'][um] % dnu,freqs['freq'][um],'o',color=colors[3])
um=np.where(freqs['l'] == 1)[0]
plt.plot(freqs['freq'][um] % dnu,freqs['freq'][um],'s',color=colors[1])
um=np.where(freqs['l'] == 2)[0]
plt.plot(freqs['freq'][um] % dnu,freqs['freq'][um],'D',color=colors[0])

#um=np.where(freqs['l'] == 0)[0]
#plt.plot(dnu+(freqs['freq'][um] % dnu),freqs['freq'][um],'o',color=colors[3])
#um=np.where(freqs['l'] == 1)[0]
#plt.plot(dnu+(freqs['freq'][um] % dnu),freqs['freq'][um],'s',color=colors[1])
#um=np.where(freqs['l'] == 2)[0]
#plt.plot(dnu+(freqs['freq'][um] % dnu),freqs['freq'][um],'D',color=colors[0])

plt.xlabel('Frequency mod '+str(dnu)[0:6]+' $\mu$Hz')
plt.ylim([125,280])
plt.ylabel('Frequency ($\mu$Hz)')


plt.subplot(1,3,3)
alf=0.3
alf2=0.1

um=np.where(freqs['l'] == 0)[0]
plt.scatter(freqs['freq'][um] % dnu,freqs['freq'][um],marker='o',color='grey',alpha=alf)
um=np.where(freqs['l'] == 1)[0]
plt.scatter(freqs['freq'][um] % dnu,freqs['freq'][um],marker='s',color='grey',alpha=alf)
um=np.where(freqs['l'] == 2)[0]
plt.scatter(freqs['freq'][um] % dnu,freqs['freq'][um],marker='D',color='grey',alpha=alf)

um=np.where(freqs_mod['l'] == 0)[0]
plt.scatter(freqs_mod['Sur.Corr.'][um] % dnu,freqs_mod['Sur.Corr.'][um],marker='o',color=colors[3],facecolors='none')
um=np.where(freqs_mod['l'] == 1)[0]
plt.scatter(freqs_mod['Sur.Corr.'][um] % dnu,freqs_mod['Sur.Corr.'][um],marker='s',color=colors[1],facecolors='none')
um=np.where(freqs_mod['l'] == 2)[0]
plt.scatter(freqs_mod['Sur.Corr.'][um] % dnu,freqs_mod['Sur.Corr.'][um],marker='D',color=colors[0],facecolors='none')

um=np.where(freqs_mod['l'] == 0)[0]
plt.scatter(freqs_mod['Uncorr'][um] % dnu,freqs_mod['Uncorr'][um],marker='o',color=colors[3],facecolors='none',alpha=alf)
um=np.where(freqs_mod['l'] == 1)[0]
plt.scatter(freqs_mod['Uncorr'][um] % dnu,freqs_mod['Uncorr'][um],marker='s',color=colors[1],facecolors='none',alpha=alf)
um=np.where(freqs_mod['l'] == 2)[0]
plt.scatter(freqs_mod['Uncorr'][um] % dnu,freqs_mod['Uncorr'][um],marker='D',color=colors[0],facecolors='none',alpha=alf)

plt.ylabel('Frequency ($\mu$Hz)')

plt.xlabel('Frequency mod '+str(dnu)[0:6]+' $\mu$Hz')
plt.ylim([125,280])

plt.tight_layout()

plt.savefig('ps.png',dpi=200)

res=(freqs_mod['Observ.']-freqs_mod['Sur.Corr.'])/freqs_mod['Observ.']
res_abs=(freqs_mod['Observ.']-freqs_mod['Sur.Corr.'])

print(np.std(res_abs),np.std(res))

np.median((freqs_mod['Observ.']-freqs_mod['Sur.Corr.'])/freqs_mod['Observ.'])
