
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



hdus = fits.open('/home/jeffrey/Documents/mass_10_2/all_mass_10_2_SFH_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/Documents/mass_11_3/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/Documents/mass_10_2/NGC5379.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_2/NGC7631.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_2/NGC6150B.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_2/UGC01918.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_2/UGC05359.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_2/UGC08778.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_2/UGC09759.p_e.rad_SFH_lum_Mass.fits.gz',]

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(7,1404)

def chi2_distance(histA, histB, eps = 1e-10):
	d = 0.5 * np.sum([((a - b) ** 2) / (a + b + eps)
		for (a, b) in zip(histA, histB)])
	return d


nx, ny = img.shape
nx1,ny1 = img_mask1.shape
imgh = np.reshape(img, nx*ny)
imgh1 = np.reshape(img_mask1, nx1*ny1)


df=36*38
dof=chi2_distance(imgh,imgh1)
chi2_distance(imgh, imght[0,:])
chi2_distance(imgh, imght[1,:])
chi2_distance(imgh, imght[2,:])
chi2_distance(imgh, imght[3,:])
chi2_distance(imgh, imght[4,:])
chi2_distance(imgh, imght[5,:])
chi2_distance(imgh, imght[6,:])

chi_sfh=np.array([502.5312995308081,580.45729191839803,204.19667370317001,518.27719309677684,1534.8645907539676,1555.9265125639893])


weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=10, normed=False,weights=weights, histtype='step', lw=3,label='Mass $10<\log (M/M_{\odot})<11$ , color $2<g-r<3$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()


