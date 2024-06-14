import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

plt.rcParams['xtick.major.size'] = 9
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 6
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 9
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 6
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size']=20
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

# use a color-blind friendly palette
# orange, red, light blue, dark blue
colors=['#FF9408','#DC4D01','#00A9E0','#016795']


#### APOKASC-2 (Pinsonneault et al. 2014)
apo=ascii.read('data/apokasc2.tsv')
apo['FeH'].fill_value=-99
apo['Numax'].fill_value=-99

### APO-K2 (Schonhut-Stasik et al. 2024)
apok2=ascii.read('data/apo_k2_public_catalog.csv')
um=np.where((apok2['feh_err'] < 0.1) & (apok2['alpham_err'] < 0.1))[0]
apok2=apok2[um]

## Kepler dwarfs & subgiants from Serenelli et al. (2017)
dw=ascii.read('data/dwarfs.csv')

# nu Indi (Chaplin et al. 2020)  
ids=['$\\nu$ Indi']
numaxv=[350.]
fehv=[-1.45]

# K2 sample from Valentini et al. (2019)
numaxv_val=[20.2,51.2,41.8,190.0]
fehv_val=[-2.01,-1.5,-1.56,-2.23]
ids_val=['EPIC201359581','EPIC205997746','EPIC206034668','EPIC206443679']
numaxv_val2=[190.0]
fehv_val2=[-2.23]
ids_val2=['EPIC206443679']

# Matsuno et al. (2021)
matsuno=ascii.read('data/matsuno.txt')

# Sample from this work
sample=ascii.read('data/sample.csv')
um=np.where((sample['feh_hds'] > -99) & (sample['feh_hires'] > -99))[0]
keep=np.zeros(len(sample))
sample_feh=np.zeros(len(sample))
for i in range(0,len(sample)):
	if ((sample['feh_hires'][i] < -1) & (sample['feh_hires'][i] > -99) & (sample['feh_hds'][i] == -99)):
		sample_feh[i]=sample['feh_hires'][i]
	if ((sample['feh_hds'][i] < -1) & (sample['feh_hds'][i] > -99) & (sample['feh_hires'][i] == -99)):
		sample_feh[i]=sample['feh_hds'][i]
	if ((sample['feh_hds'][i] < -1) & (sample['feh_hds'][i] > -99) & (sample['feh_hires'][i] > -99) & (sample['feh_hires'][i] < -1)):
		sample_feh[i]=(sample['feh_hds'][i]+sample['feh_hires'][i])/2

um=np.where(sample_feh != 0.)[0]
sample=sample[um]
sample_feh=sample_feh[um]

keep=np.zeros(len(sample))
for i in range(0,len(sample)):
	um=np.where(matsuno['col1'] == sample['kic'][i])[0]
	um2=np.where(apo['KIC'] == sample['kic'][i])[0]
	if ((len(um) ==0) & (len(um2) == 0)):
		keep[i]=1
		
um=np.where(keep == 1)[0]
sample=sample[um]
sample_feh=sample_feh[um]

um=np.where(sample['kic'] != 8144907)[0]
sample=sample[um]
sample_feh=sample_feh[um]

for i in range(0,len(sample)):
	if (sample['numax'][i] > 0):
		print(sample['kic'][i],'&',sample['numax'][i],'&',sample_feh[i])



ss=3

plt.ion()
plt.clf()

plt.subplot(1,2,1)
plt.semilogx(apo['Numax'],apo['FeH'],'o',color='black',ms=ss)
plt.semilogx(apok2['numax'],apok2['feh'],'o',color='black',ms=ss)

plt.semilogx(dw['col0'],dw['col1'],'o',color='black',ms=ss)
plt.semilogx(numaxv,fehv,'o',color=colors[2])
plt.semilogx(numaxv_val,fehv_val,'o',color='black',ms=ss)
plt.semilogx(numaxv_val2,fehv_val2,'o',color=colors[3])
plt.semilogx(matsuno['col2'],matsuno['col3'],'o',color='black',ms=ss)
plt.semilogx(sample['numax'],sample_feh,'o',color=colors[1])
plt.semilogx([328],[-1.62],'o',color=colors[0])

plt.ylabel('[Fe/H]')
plt.xlabel('Frequency of Maximum Power ($\mu$Hz)')
plt.ylim([-3,0.6])

plt.annotate(ids[0],(numaxv[0]+0.1*numaxv[0],fehv[0]),xycoords='data',fontsize=11,color=colors[2])
#plt.annotate(ids[1],(numaxv[1]+0.002*numaxv[0],fehv[1]),xycoords='data',fontsize=11,color='orange')
plt.annotate(ids_val[3],(numaxv_val[3]+0.1*numaxv_val[3],fehv_val[3]),xycoords='data',fontsize=11,color=colors[3])
plt.annotate('KIC7341231',(350,-1.6),xycoords='data',fontsize=11,color=colors[0])

plt.annotate('KIC8144907',(230,-2.6),xycoords='data',fontsize=11,color=colors[1])
plt.plot([212],[-2.66],'*',ms=20,color=colors[1])

plt.annotate('KIC5938297',(25,-2.85),xycoords='data',fontsize=11,color=colors[1])


plt.axvspan(0,100, color='y', alpha=0.5, lw=0)

plt.subplot(1,2,2)
um=np.where(apo['Numax'].filled() > 100.)[0]
plt.plot(apo['FeH'][um],apo['AFe'][um],'o',color='black',ms=ss)
plt.plot(apok2['feh'],apok2['alpham'],'o',color='black',ms=ss)

plt.xlabel('[Fe/H]')
plt.ylabel(r'[$\alpha/Fe$]')

plt.plot([-2.66],[0.38],'*',ms=20,color=colors[1],label='KIC8144907')
#plt.annotate('KIC8144907',(-2.58,0.38),xycoords='data',fontsize=11,color=colors[1])

plt.plot([-1.62],[0.32],'o',color=colors[0],label='KIC7341231')
#plt.annotate('KIC7341231',(-1.57,0.3),xycoords='data',fontsize=11,color=colors[0])

plt.plot([-1.43],[0.32],'o',color=colors[2],label='$\\nu$ Indi')
#plt.annotate('$\\nu$ Indi',(-1.38,0.32),xycoords='data',fontsize=11,color=colors[2])

plt.plot([-2.23],[0.23],'o',color=colors[3],label='EPIC206443679')
#plt.annotate('EPIC206443679',(-2.18,0.23),xycoords='data',fontsize=11,color=colors[3])
plt.legend(fontsize=13)
plt.ylim([-0.1,0.45])
#plt.ylim([-0.1,0.65])
#plt.xlim([-3,0.5])

plt.tight_layout()

#plt.savefig('fig2.png',dpi=200)

