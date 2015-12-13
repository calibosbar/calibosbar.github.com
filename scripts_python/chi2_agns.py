

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



hdus = fits.open('/home/jeffrey/tesis/AGNS/all_agns_chi.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/tesis/mass_11_3/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/tesis/AGNS/IC0540.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/MCG-02-02-030.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC0426.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC0833.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC1667.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC2410.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC2623.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC2639.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC3106.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC3160.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC3303.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC5635.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC5675.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC5739.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC6323.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC6394.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC7466.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/NGC7738.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/UGC00005.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/UGC00987.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/UGC03995.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/UGC05771.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/UGC09711.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/AGNS/UGC03973.p_e.rad_SFH_lum_Mass.fits.gz']

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(24,1404)

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
chi2_distance(imgh, imght[7,:])
chi2_distance(imgh, imght[8,:])
chi2_distance(imgh, imght[9,:])
chi2_distance(imgh, imght[10,:])
chi2_distance(imgh, imght[11,:])
chi2_distance(imgh, imght[12,:])
chi2_distance(imgh, imght[13,:])
chi2_distance(imgh, imght[14,:])
chi2_distance(imgh, imght[15,:])
chi2_distance(imgh, imght[16,:])
chi2_distance(imgh, imght[17,:])
chi2_distance(imgh, imght[18,:])
chi2_distance(imgh, imght[19,:])
chi2_distance(imgh, imght[20,:])
chi2_distance(imgh, imght[21,:])
chi2_distance(imgh, imght[22,:])
chi2_distance(imgh, imght[23,:])

chi_sfh=np.array([672.47162721967879,1205.1883751002006,1827.4363180724549,2920.0802417989307,3650.1448984639828,643.25289248599643,2510.2515829297854,
9889.0137201549878,1021.6477798171122,545.62641216777536,560.62492118763612,1024.0544918973637,200.72236833204934,1263.0929406712014,902.31479427960994,
1221.1222834800806,886.15540777789283,913.80484648649815,746.32569818753075,2895.3268688997964,803.02163381714956,623.09225434519669,1061.6856004611554,
4385.7822934832575])





weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=30, normed=False,weights=weights, histtype='step', lw=3,label='Mass $9<\log (M/M_{\odot})<10$ , color $1<g-r<2$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()


