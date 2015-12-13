

from scipy.stats import chi2
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure, show, rc
from scipy.stats import chisquare
from scipy.stats import norm
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
from astropy.io import fits
import pylab
import seaborn as sn
import matplotlib.mlab as mlab
import math





hdus = fits.open('all_mass_8_1_SFH_Mass.fits.gz')
img = hdus[0].data

hdus1 = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)

hdus2 = fits.open('IC2095.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus3 = fits.open('MCG-02-06-016.p_e.rad_SFH_lum_Mass.fits.gz')
img3 = hdus3[0].data
img_mask3= np.power(10, img3)

hdus4 = fits.open('NGC3353.p_e.rad_SFH_lum_Mass.fits.gz')
img4 = hdus4[0].data
img_mask4= np.power(10, img4)

hdus5 = fits.open('NGC3773.p_e.rad_SFH_lum_Mass.fits.gz')
img5 = hdus5[0].data
img_mask5= np.power(10, img5)

hdus6 = fits.open('NGC3896.p_e.rad_SFH_lum_Mass.fits.gz')
img6 = hdus6[0].data
img_mask6= np.power(10, img6)

hdus7 = fits.open('UGC10650.p_e.rad_SFH_lum_Mass.fits.gz')
img7 = hdus7[0].data
img_mask7= np.power(10, img7)


nx, ny = img.shape
nx1,ny1 = img_mask1.shape
nx2,ny2 = img_mask2.shape
nx3,ny3 = img_mask3.shape
nx4,ny4 = img_mask4.shape
nx5,ny5 = img_mask5.shape
nx6,ny6 = img_mask6.shape
nx7,ny7 = img_mask7.shape


imgh = np.reshape(img, nx*ny)
imgh1 = np.reshape(img_mask1, nx1*ny1)
imgh2 = np.reshape(img_mask2, nx2*ny2)
imgh3 = np.reshape(img_mask3,nx3*ny3)
imgh4 = np.reshape(img_mask4,nx4*ny4)
imgh5 = np.reshape(img_mask5,nx5*ny5)
imgh6 = np.reshape(img_mask6,nx6*ny6)
imgh7 = np.reshape(img_mask7,nx7*ny7)

def chi2_distance(histA, histB, eps = 1e-10):
	d = 0.5 * np.sum([((a - b) ** 2) / (a + b + eps)
		for (a, b) in zip(histA, histB)])
	return d

df=36*38

dof=chi2_distance(imgh,imgh1)
chi2_distance(imgh,imgh2)
chi2_distance(imgh,imgh3)
chi2_distance(imgh,imgh4)
chi2_distance(imgh,imgh5)
chi2_distance(imgh,imgh6)
chi2_distance(imgh,imgh7)


chi_sfh=np.array([1439.1873977893379,3047.2940747101552,774.67561500830425,897.71091148083065, 517.56654522371889])
#1187.7377355082149
weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
df = 1368
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01')
ax.hist(chi_sfh/df,bins=13, normed=False,weights=weights, histtype='step', lw=3,label='Mass range $8<\log (M/M_{\odot})<9$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()

