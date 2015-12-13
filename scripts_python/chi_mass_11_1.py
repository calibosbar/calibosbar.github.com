




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



hdus = fits.open('/home/jeffrey/Documents/mass_11_1/all_mass_11_1_SFH_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/Documents/mass_11_2/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/Documents/mass_11_1/IC1078.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0023.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0165.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0214.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0234.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0257.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0309.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0768.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0776.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0835.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC0976.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC1094.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC1142.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC2347.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC4185.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC5000.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC5056.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC5947.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC6301.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC7321.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC7489.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/NGC7653.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/UGC00005.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/UGC01659.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/UGC02311.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/UGC12767.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/UGC12810.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_1/MCG-01-09-006.p_e.rad_SFH_lum_Mass.fits.gz']

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(28,1404)

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
chi2_distance(imgh, imght[24,:])
chi2_distance(imgh, imght[25,:])
chi2_distance(imgh, imght[26,:])
chi2_distance(imgh, imght[27,:])



chi_sfh=np.array([1606.9905047423656,1178.5242856122097,844.51447527036623,711.80202671126574,571.40591313717778,914.77256010030146,
250.75376679608681,315.23482588599967,3885.0959829745429,2254.0395467089406,799.70556773379599,2180.0159454836926,1843.0874877664567,
822.97951604021193,265.52059810138906,297.07942916209117,766.47342705889571,1899.8640333375133,387.9778117934905,1773.5168553321578,
480.80792417678316,539.75778477550523,1256.9019718434627,438.59937127693973,1602.3900820984384,991.53901508915101,321.61274723238915])



weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=30, normed=False,weights=weights, histtype='step', lw=3,label='Mass $11<\log (M/M_{\odot})<12$ , color $1<g-r<2$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()


