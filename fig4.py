import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sb
import match

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


# seismic exoplanet hosts from Silva Aguirre et al. 2015
exo=ascii.read('data/hosts.tsv')
x=(exo['E_Age']+exo['e_Age'])/2
um=np.where(x < 1)[0]
exo=exo[um]

# Legacy sample
legacy1=ascii.read('data/legacy-t1.tsv')
legacy2=ascii.read('data/legacy-t2.tsv')
um=np.where(legacy2['Pipe'] == 'BASTA')[0]
legacy2=legacy2[um]
ix,iy=match.match(legacy1['KIC'],legacy2['KIC'])
legacy_age=legacy2['Age'][iy]
legacy_age_err=(legacy2['E_Age'][iy]+legacy2['e_Age'][iy])/2
legacy_feh=legacy1['[Fe/H]'][ix]
legacy_feh_err=legacy1['e_[Fe/H]'][ix]

# subgiants from Li et al. 2020
subgiants=ascii.read('data/li-2020.csv')
um=np.where(subgiants['e_Age'] < 1)[0]
subgiants=subgiants[um]

## Xiang & Rix (2022) subgiants
lamost = fits.open('data/subgiant_fullparam_update.fits')

## Accreted ##
accreted = np.where((lamost[1].data['ECC'] > 0.9) | ((lamost[1].data['ECC'] > 0.8) & (lamost[1].data['LZ'] < -2000)) )[0]

## In-situ ##
insitu = np.where((lamost[1].data['ECC'] < 0.8) & (lamost[1].data['LZ'] > 0))[0]

cm=plt.get_cmap('Greys')

# use a color-blind friendly palette
# orange, red, light blue, dark blue
colors=['#FF9408','#DC4D01','#00A9E0','#016795']
alf=0.8

ss=8
alf=0.4

col1=colors[1]
col2=colors[1]

plt.ion()
plt.clf()
plt.subplot(1,2,1)
#plt.grid()
plt.hexbin(lamost[1].data['AGE'][insitu],lamost[1].data['FEH'][insitu],gridsize=100,cmap=cm,bins='log')
plt.xlim([0,15])
plt.ylim([-3,0.6])

plt.ylabel('[Fe/H]')
plt.xlabel('Age (Gyr)')

plt.errorbar(exo['Age'],exo['[Fe/H]'],yerr=exo['e_[Fe/H]'],xerr=(exo['E_Age']+exo['e_Age'])/2,fmt='.',color=col1,alpha=alf)
plt.plot(exo['Age'],exo['[Fe/H]'],'o',color=col1,ms=ss,label='Asteroseismology')
plt.plot(-1,0,'o',color='grey',ms=ss,label='Isochrone Fitting')

um=np.where(legacy_age_err < 1)[0]
plt.errorbar(legacy_age[um],legacy_feh[um],yerr=legacy_feh_err[um],xerr=legacy_age_err[um],fmt='.',color=col1,alpha=alf)
plt.plot(legacy_age[um],legacy_feh[um],'o',color=col1,ms=ss)

plt.errorbar(subgiants['Age'],subgiants['MH'],yerr=subgiants['e_MH'],xerr=subgiants['e_Age'],fmt='.',color=col1,alpha=alf)
plt.plot(subgiants['Age'],subgiants['MH'],'o',color=col1,ms=ss)

plt.annotate('In-Situ',(11.75,0.25),xycoords='data')

# Kepler-444 (Campante et al. 2015)
plt.errorbar([11.23],[-0.55],xerr=[0.95],yerr=[0.07],color=col1,alpha=alf)
plt.plot([11.23],[-0.55],'o',color=col1,ms=ss)
# nu Indi (Chaplin et al. 2020)
plt.errorbar([11.0],[-1.43],xerr=[0.8],yerr=[0.07],color=col1,alpha=alf)
plt.plot([11.0],[-1.43],'o',color=col1,ms=ss)

plt.errorbar([12.0],[-2.66],xerr=[0.68],yerr=[0.07],color=col1,alpha=alf)
plt.plot([12.0],[-2.66],'o',color=col1,ms=ss)

plt.legend(loc='lower left')

plt.subplot(1,2,2)
#plt.grid()
plt.hexbin(lamost[1].data['AGE'][accreted],lamost[1].data['FEH'][accreted],gridsize=100,cmap=cm,bins='log')
plt.xlim([0,15])
plt.ylim([-3,0.6])
plt.errorbar([12.0],[-2.66],xerr=[0.68],yerr=[0.07],color=col2,alpha=alf)
plt.plot([12.0],[-2.66],'o',color=col2,ms=ss,label='Asteroseismology')
plt.plot(-1,0,'o',color='grey',ms=ss,label='Isochrone Fitting')

plt.ylabel('[Fe/H]')
plt.xlabel('Age (Gyr)')
plt.annotate('Accreted',(11.75,0.25),xycoords='data')
plt.legend(loc='lower left')

plt.tight_layout()



plt.savefig('fig4-v2.png',dpi=200)
plt.savefig('fig4.pdf')





###




plt.ion()
plt.clf()
plt.subplot(2,2,1)
plt.hexbin(lamost[1].data['FEH'][insitu],lamost[1].data['ALPHA_FE'][insitu],gridsize=100,cmap=cm,bins='log')
plt.xlim([-3,0.5])
plt.ylim([-0.1,0.5])
#plt.errorbar([-2.66],[0.38],xerr=[0.8],yerr=[0.07],color=colors[1],alpha=alf)
plt.plot([-2.66],[0.38],'*',color=colors[1],ms=12)
plt.xlabel('[Fe/H]')
plt.ylabel(r'[$\alpha/Fe$]')

plt.subplot(2,2,2)
plt.hexbin(lamost[1].data['FEH'][accreted],lamost[1].data['ALPHA_FE'][accreted],gridsize=100,cmap=cm,bins='log')
plt.xlim([-3,0.5])
plt.ylim([-0.1,0.5])
#plt.errorbar([12.0],[-2.66],xerr=[0.8],yerr=[0.07],color=colors[1],alpha=alf)
plt.plot([-2.66],[0.38],'*',color=colors[1],ms=12)
plt.xlabel('[Fe/H]')
plt.ylabel(r'[$\alpha/Fe$]')

plt.subplot(2,2,3)
plt.hexbin(lamost[1].data['AGE'][insitu],lamost[1].data['FEH'][insitu],gridsize=100,cmap=cm,bins='log')
plt.xlim([0,15])
plt.ylim([-3,0.5])
plt.errorbar([12.0],[-2.66],xerr=[0.8],yerr=[0.07],color=colors[1],alpha=alf)
plt.plot([12.0],[-2.66],'*',color=colors[1],ms=12)
plt.ylabel('[Fe/H]')
plt.xlabel('Age (Gyr)')

# Kepler-444
plt.errorbar([11.23],[-0.55],xerr=[0.95],yerr=[0.07],color=colors[0],alpha=alf)
plt.plot([11.23],[-0.55],'o',color=colors[0],ms=ss)
# nu Indi
plt.errorbar([11.0],[-1.43],xerr=[0.8],yerr=[0.07],color=colors[0],alpha=alf)
plt.plot([11.0],[-1.43],'o',color=colors[0],ms=ss)

plt.subplot(2,2,4)
plt.hexbin(lamost[1].data['AGE'][accreted],lamost[1].data['FEH'][accreted],gridsize=100,cmap=cm,bins='log')
plt.xlim([0,15])
plt.ylim([-3,0.5])
plt.errorbar([12.0],[-2.66],xerr=[0.8],yerr=[0.07],color=colors[1],alpha=alf)
plt.plot([12.0],[-2.66],'*',color=colors[1],ms=12)
plt.ylabel('[Fe/H]')
plt.xlabel('Age (Gyr)')

# Kepler-444
plt.errorbar([11.23],[-0.55],xerr=[0.95],yerr=[0.07],color=colors[0],alpha=alf)
plt.plot([11.23],[-0.55],'o',color=colors[0],ms=ss)
# nu Indi
plt.errorbar([11.0],[-1.43],xerr=[0.8],yerr=[0.07],color=colors[0],alpha=alf)
plt.plot([11.0],[-1.43],'o',color=colors[0],ms=ss)




plt.ion()
plt.clf()
plt.subplot(1,2,1)
plt.hexbin(lamost[1].data['AGE'][insitu],lamost[1].data['FEH'][insitu],gridsize=100,cmap=cm,bins='log')
plt.xlim([0,15])
plt.ylim([-3,0.5])
plt.errorbar([12.0],[-2.66],xerr=[0.8],yerr=[0.07],color=colors[1],alpha=alf)
plt.plot([12.0],[-2.66],'*',color=colors[1],ms=12)
plt.ylabel('[Fe/H]')
plt.xlabel('Age (Gyr)')

plt.subplot(1,2,2)
plt.hexbin(lamost[1].data['AGE'][accreted],lamost[1].data['FEH'][accreted],gridsize=100,cmap=cm,bins='log')
plt.xlim([0,15])
plt.ylim([-3,0.5])
plt.errorbar([12.0],[-2.66],xerr=[0.8],yerr=[0.07],color=colors[1],alpha=alf)
plt.plot([12.0],[-2.66],'*',color=colors[1],ms=12)
plt.ylabel('[Fe/H]')
plt.xlabel('Age (Gyr)')

plt.tight_layout()


plt.plot(lamost[1].data['AGE'][insitu],lamost[1].data['FEH'][insitu],'.')
plt.plot(lamost[1].data['AGE'][accreted],lamost[1].data['FEH'][accreted],'.')


plt.ion()
plt.clf()
plt.plot(lamost[1].data['AGE'][insitu],lamost[1].data['FEH'][insitu],'.',ms=0.5,color='grey',alpha=0.1)

sb.kdeplot(lamost[1].data['AGE'][insitu],lamost[1].data['FEH'][insitu],gridsize=50,log_scale=True)

plt.plot(lamost[1].data['AGE'][accreted],lamost[1].data['FEH'][accreted],'.',ms=0.5,color='red',alpha=0.1)
sb.kdeplot(lamost[1].data['AGE'][accreted],lamost[1].data['FEH'][accreted],gridsize=50,log_scale=True,color='red')

sb.kdeplot(lamost[1].data['AGE'][accreted],lamost[1].data['FEH'][accreted])


plt.contour(lamost[1].data['FEH'][insitu],lamost[1].data['AGE'][insitu])
plt.hexbin(lamost[1].data['FEH'][accreted],lamost[1].data['AGE'][accreted],gridsize=100,cmap=cm,bins='log')

plt.ion()
plt.clf()
plt.subplot(2,2,1)
plt.hexbin(lamost[1].data['FEH'][insitu],lamost[1].data['ALPHA_FE'][insitu],gridsize=100,cmap=cm,bins='log')
plt.subplot(2,2,2)
plt.hexbin(lamost[1].data['FEH'][accreted],lamost[1].data['ALPHA_FE'][accreted],gridsize=100,cmap=cm,bins='log')

plt.subplot(2,2,3)
plt.hexbin(lamost[1].data['FEH'][insitu],lamost[1].data['AGE'][insitu],gridsize=100,cmap=cm,bins='log')
plt.subplot(2,2,4)
plt.hexbin(lamost[1].data['FEH'][accreted],lamost[1].data['AGE'][accreted],gridsize=100,cmap=cm,bins='log')


## Plotting ##

fig = plt.figure(figsize=(9, 9))

ax2a, ax3a = fig.add_subplot(221), fig.add_subplot(222)
ax2, ax3,=  fig.add_subplot(223), fig.add_subplot(224)

im3a = ax3a.scatter(lamost[1].data['FEH'][insitu],lamost[1].data['ALPHA_FE'][insitu],
                  c=lamost[1].data['ECC'][insitu], s=3, vmin=0, vmax=1)
ax3a.set_yticklabels([])
ax3a.set_xlabel('[Fe/H] (dex)')
ax3a.text(x=0.25, y=0.05, s='In-situ',fontsize=15, transform=ax3a.transAxes, ha='center')

im2a = ax2a.scatter( lamost[1].data['FEH'][accreted],lamost[1].data['ALPHA_FE'][accreted],
                  c=lamost[1].data['ECC'][accreted], s=3, vmin=0, vmax=1)
ax2a.set_xlim(ax3a.get_xlim())
ax2a.set_ylim(ax3a.get_ylim())
ax2a.set_xlabel('[Fe/H] (dex)')
ax2a.set_ylabel('[$\\alpha$/Fe] (dex)')
ax2a.text(x=0.25, y=0.05, s='Accreted',fontsize=15, transform=ax2a.transAxes, ha='center')


im2 = ax2.scatter(lamost[1].data['AGE'][accreted], lamost[1].data['FEH'][accreted],
                  c=lamost[1].data['ECC'][accreted], s=3, vmin=0, vmax=1)

ax2.set_ylim(ax1.get_ylim())
ax2.set_xlim(ax1.get_xlim())
ax2.set_xlabel('Age (Gyr)')
ax2.set_ylabel('[Fe/H]')
ax2.text(x=0.25, y=0.05, s='Accreted',fontsize=15, transform=ax2.transAxes, ha='center')

im3 = ax3.scatter(lamost[1].data['AGE'][insitu], lamost[1].data['FEH'][insitu],
                  c=lamost[1].data['ECC'][insitu], s=3, vmin=0, vmax=1)

ax3.set_ylim(ax1.get_ylim())
ax3.set_xlim(ax1.get_xlim())
ax3.text(x=0.25, y=0.05, s='In-situ',fontsize=15,
         transform=ax3.transAxes, ha='center')
ax3.set_yticklabels([])
cb = fig.colorbar(im3, ax=ax3,
                   orientation='horizontal', pad=0.15, cax = fig.add_axes([0.08, 1.0, 0.9, 0.015]))

cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')

cb.set_label('Eccentricity')
plt.tight_layout(h_pad=0.5, w_pad=0.1)
ax3.set_xlabel('Age (Gyr)')
plt.show()