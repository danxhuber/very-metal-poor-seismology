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

plt.ion()
plt.clf()

hds=ascii.read('data/KIC8144907_spec.txt')
plt.plot(hds['waveobs'],hds['flux']+0.3,lw=1.5,color=colors[3])

solar=ascii.read('data/KIC10006158_spec.txt')
plt.plot(solar['waveobs'],solar['flux']-0.3,lw=1.5,color=colors[1])

plt.xlim([516.65,519.7])
plt.ylim([-0.7,1.8])

plt.annotate('KIC8144907 ([Fe/H] = -2.6)',(516.7,1.5),xycoords='data',color=colors[3])
plt.annotate('KIC10006158 ([Fe/H] = 0.0)',(516.7,-0.5),xycoords='data',color=colors[1])

plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalized Flux')
plt.show()
plt.draw()

plt.tight_layout()

plt.savefig('fig1.png',dpi=200)




