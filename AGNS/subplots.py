

from scipy.stats import chi2
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure, show, rc
from scipy.stats import chisquare
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
from astropy.io import fits
import pylab





hdus1 = fits.open('all_mass_8_1_SFH_Mass.fits.gz')
img1 = hdus1[0].data

hdus2 = fits.open('all_mass_9_1_SFH_Mass.fits.gz')
img2 = hdus2[0].data

hdus3 = fits.open('all_mass_9_2_SFH_Mass.fits.gz')
img3 = hdus3[0].data

hdus4 = fits.open('all_mass_9_3_SFH_Mass.fits.gz')
img4 = hdus4[0].data

hdus5 = fits.open('all_mass_10_1_SFH_Mass.fits.gz')
img5 = hdus5[0].data

hdus6 = fits.open('all_mass_10_2_SFH_Mass.fits.gz')
img6 = hdus6[0].data

hdus7 = fits.open('all_mass_10_3_SFH_Mass.fits.gz')
img7 = hdus7[0].data

hdus8 = fits.open('all_mass_11_1_SFH_Mass.fits.gz')
img8 = hdus8[0].data

hdus9 = fits.open('all_mass_11_2_SFH_Mass.fits.gz')
img9 = hdus9[0].data

hdus10 = fits.open('all_mass_11_3_SFH_Mass.fits.gz')
img10 = hdus10[0].data

hdus11 = fits.open('all_mass_11_4_SFH_Mass.fits.gz')
img11 = hdus11[0].data




fig, ((ax1,ax2,ax3,ax4), (ax5,ax6,ax7,ax8),(ax9,ax10,ax11,ax12),(ax13,ax14,ax15,ax16))= plt.subplots(4, 4, sharex='col', sharey='row')

#ticksy13=np.linspace(1,2,9)
#ticksy9=np.linspace(2.1,3,9)
#ticksy5=np.linspace(3.1,4,9)
#ticksy1=np.linspace(4.1,5,9)
#ticksx13=np.linspace(8,9,9)
#ticksx14=np.linspace(9.1,10,9)
#ticksx15=np.linspace(10.1,11,9)
#ticksx16=np.linspace(11.1,12,9)

#ax13.set_xticklabels(ticksx13)
#ax14.set_xticklabels(ticksx14)
#ax15.set_xticklabels(ticksx15)
#ax16.set_xticklabels(ticksx16)
#ax13.set_yticklabels(ticksy13)
#ax9.set_yticklabels(ticksy9)
#ax5.set_yticklabels(ticksy5)
#ax1.set_yticklabels(ticksy1)

#ax13.axes.get_xaxis().set_visible(False)
#fig.text(0.5, 0.04, '$\log (M/M_{\odot})$', ha='center', va='center')
#fig.text(0.06, 0.5, '$g-r$', ha='center', va='center', rotation='vertical')
#ax13.text(0.05, 0.004, '$8<\log (M/M_{\odot})<9$', ha='center', va='center')
ax13.axes.get_xaxis().set_ticks([])
ax14.axes.get_xaxis().set_ticks([])
ax15.axes.get_xaxis().set_ticks([])
ax16.axes.get_xaxis().set_ticks([])
ax13.axes.get_yaxis().set_ticks([])
ax9.axes.get_yaxis().set_ticks([])
ax5.axes.get_yaxis().set_ticks([])
ax1.axes.get_yaxis().set_ticks([])
ax13.set_xlabel('$8<\log (M/M_{\odot})<9$')
ax14.set_xlabel('$9<\log (M/M_{\odot})<10$')
ax15.set_xlabel('$10<\log (M/M_{\odot})<11$')
ax16.set_xlabel('$11<\log (M/M_{\odot})<12$')
ax13.set_ylabel('$1<g-r<2$')
ax9.set_ylabel('$2<g-r<3$')
ax5.set_ylabel('$3<g-r<4$')
ax1.set_ylabel('$4<g-r<5$')



ax13.imshow(img1, origin = 'lower',aspect='auto') 
ax13.pcolor(img1, norm=LogNorm(),cmap=califa)
ax14.imshow(img2, origin = 'lower',aspect='auto') 
ax14.pcolor(img2, norm=LogNorm(),cmap=califa)
ax10.imshow(img3, origin = 'lower',aspect='auto') 
ax10.pcolor(img3, norm=LogNorm(),cmap=califa)
ax6.imshow(img4, origin = 'lower',aspect='auto') 
ax6.pcolor(img4, norm=LogNorm(),cmap=califa)
ax15.imshow(img5, origin = 'lower',aspect='auto') 
ax15.pcolor(img5, norm=LogNorm(),cmap=califa)
ax11.imshow(img6, origin = 'lower',aspect='auto') 
ax11.pcolor(img6, norm=LogNorm(),cmap=califa)
ax7.imshow(img7, origin = 'lower',aspect='auto') 
ax7.pcolor(img7, norm=LogNorm(),cmap=califa)
ax16.imshow(img8, origin = 'lower',aspect='auto') 
ax16.pcolor(img8, norm=LogNorm(),cmap=califa)
ax12.imshow(img9, origin = 'lower',aspect='auto') 
ax12.pcolor(img9, norm=LogNorm(),cmap=califa)
ax8.imshow(img10, origin = 'lower',aspect='auto') 
ax8.pcolor(img10, norm=LogNorm(),cmap=califa)
ax4.imshow(img11, origin = 'lower',aspect='auto') 
c1=ax4.pcolor(img11, norm=LogNorm(),cmap=califa)


fig.subplots_adjust(right=0.87)
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
cbar=fig.colorbar(c1, cax=cax)
cbar.ax.set_ylabel('Log $\Sigma_{*}$ $[M_{\odot} kpc^{-2}]$', fontsize=12)




