
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


hdus = fits.open('all_agns_SFH_Mass.fits.gz')
img = hdus[0].data

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

hdus12 = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img12 = hdus12[0].data
img_masked1= np.power(10, img12)


nx, ny = img.shape
nx1,ny1 = img1.shape
nx2,ny2 = img2.shape
nx3,ny3 = img3.shape
nx4,ny4 = img4.shape
nx5,ny5 = img5.shape
nx6,ny6 = img6.shape
nx7,ny7 = img7.shape
nx8,ny8 = img8.shape
nx9,ny9 = img9.shape
nx10,ny10 = img10.shape
nx11,ny11 = img11.shape
nx12,ny12 = img_masked1.shape

imgh = np.reshape(img, nx*ny)
imgh1 = np.reshape(img1, nx1*ny1)
imgh2 = np.reshape(img2, nx2*ny2)
imgh3 = np.reshape(img3,nx3*ny3)
imgh4 = np.reshape(img4,nx4*ny4)
imgh5 = np.reshape(img5,nx5*ny5)
imgh6 = np.reshape(img6,nx6*ny6)
imgh7 = np.reshape(img7,nx7*ny7)
imgh8 = np.reshape(img8,nx8*ny8)
imgh9 = np.reshape(img9,nx9*ny9)
imgh10 = np.reshape(img10,nx10*ny10)
imgh11 = np.reshape(img11,nx11*ny11)
imgh12 = np.reshape(img_masked1,nx12*ny12)

def chi2_distance(histA, histB, eps = 1e-10):
	d = 0.5 * np.sum([((a - b) ** 2) / (a + b + eps)
		for (a, b) in zip(histA, histB)])
	return d



dof=chi2_distance(imgh12,imgh1)
dof1=chi2_distance(imgh12,imgh2)
dof2=chi2_distance(imgh12,imgh3)
dof3=chi2_distance(imgh12,imgh4)
dof4=chi2_distance(imgh12,imgh5)
dof5=chi2_distance(imgh12,imgh6)
dof6=chi2_distance(imgh12,imgh7)
dof7=chi2_distance(imgh12,imgh8)
dof8=chi2_distance(imgh12,imgh9)
dof9=chi2_distance(imgh12,imgh10)
dof10=chi2_distance(imgh12,imgh11)
dof11=chi2_distance(imgh12,imgh12)
dof12=chi2_distance(imgh12,imgh)


df = 1368


fig, ax = plt.subplots(1, 1)
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
xu1 = np.linspace(chi2.ppf(0.01, dof1), chi2.ppf(0.99, dof1), 100)
xu2 = np.linspace(chi2.ppf(0.01, dof2), chi2.ppf(0.99, dof2), 100)
xu3 = np.linspace(chi2.ppf(0.01, dof3), chi2.ppf(0.99, dof3), 100)
xu4 = np.linspace(chi2.ppf(0.01, dof4), chi2.ppf(0.99, dof4), 100)
xu5 = np.linspace(chi2.ppf(0.01, dof5), chi2.ppf(0.99, dof5), 100)
xu6 = np.linspace(chi2.ppf(0.01, dof6), chi2.ppf(0.99, dof6), 100)
xu7 = np.linspace(chi2.ppf(0.01, dof7), chi2.ppf(0.99, dof7), 100)
xu8 = np.linspace(chi2.ppf(0.01, dof8), chi2.ppf(0.99, dof8), 100)
xu9 = np.linspace(chi2.ppf(0.01, dof9), chi2.ppf(0.99, dof9), 100)
xu10 = np.linspace(chi2.ppf(0.01, dof10), chi2.ppf(0.99, dof10), 100)
xu11 = np.linspace(chi2.ppf(0.01, dof11), chi2.ppf(0.99, dof11), 100)
xu12 = np.linspace(chi2.ppf(0.01, dof12), chi2.ppf(0.99, dof12), 100)

ax.plot(xu/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $8< \log (M/M_{\odot})<9$, g-r=1')
ax.plot(xu1/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $9< \log (M/M_{\odot})<10$, g-r=1')
ax.plot(xu2/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $9< \log (M/M_{\odot})<10$, g-r=2')
ax.plot(xu3/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $9< \log (M/M_{\odot})<10$, g-r=3')
ax.plot(xu4/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $10< \log (M/M_{\odot})<11$, g-r=1')
ax.plot(xu5/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $10< \log (M/M_{\odot})<11$, g-r=2')
ax.plot(xu6/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $10< \log (M/M_{\odot})<11$, g-r=3')
ax.plot(xu7/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $11< \log (M/M_{\odot})<12$, g-r=1')
ax.plot(xu8/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $11< \log (M/M_{\odot})<12$, g-r=2')
ax.plot(xu9/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $11< \log (M/M_{\odot})<12$, g-r=3')
ax.plot(xu10/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ $11< \log (M/M_{\odot})<12$, g-r=4')
ax.plot(xu11/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ reduced UGC11680NED01 PDF')
ax.plot(xu12/df, chi2.pdf(x, df), lw=3, label='$\chi^2$ AGNs')
ax.plot(x/df, chi2.pdf(x, df), 'k--', lw=3, label='Perfect Fit')


ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()