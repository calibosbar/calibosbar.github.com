

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



hdus = fits.open('/home/jeffrey/Documents/mass_11_2/all_mass_11_2_SFH_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/Documents/mass_11_2/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/Documents/mass_11_2/ARP220.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/IC0307.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/IC0485.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/IC0674.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/IC4566.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0036.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0169.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0180.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0192.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0426.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0477.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0508.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC0787.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC1070.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC1324.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC2410.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC2449.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC2565.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC2572.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC2639.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC2916.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC4003.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5406.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5533.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5616.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5635.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5675.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5720.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5772.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC5888.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC6060.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC6154.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC6394.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC6478.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC6497.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC6941.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC7364.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC7466.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/NGC7591.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC00036.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC01368.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC02099.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC02367.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC03038.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC03973.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC03995.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC04132.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC04262.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC04455.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC05108.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC05111.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC08107.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC08234.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC08781.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC09401.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC09537.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC10205.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC10337.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC10710.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC10811.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC11717.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC12185.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC12250.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/UGC12348.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_2/MCG-02-03-015.p_e.rad_SFH_lum_Mass.fits.gz']

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(65,1404)

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
chi2_distance(imgh, imght[28,:])
chi2_distance(imgh, imght[29,:])
chi2_distance(imgh, imght[30,:])
chi2_distance(imgh, imght[31,:])
chi2_distance(imgh, imght[32,:])
chi2_distance(imgh, imght[33,:])
chi2_distance(imgh, imght[34,:])
chi2_distance(imgh, imght[35,:])
chi2_distance(imgh, imght[36,:])
chi2_distance(imgh, imght[37,:])
chi2_distance(imgh, imght[38,:])
chi2_distance(imgh, imght[39,:])
chi2_distance(imgh, imght[40,:])
chi2_distance(imgh, imght[41,:])
chi2_distance(imgh, imght[42,:])
chi2_distance(imgh, imght[43,:])
chi2_distance(imgh, imght[44,:])
chi2_distance(imgh, imght[45,:])
chi2_distance(imgh, imght[46,:])
chi2_distance(imgh, imght[47,:])
chi2_distance(imgh, imght[48,:])
chi2_distance(imgh, imght[49,:])
chi2_distance(imgh, imght[50,:])
chi2_distance(imgh, imght[51,:])
chi2_distance(imgh, imght[52,:])
chi2_distance(imgh, imght[53,:])
chi2_distance(imgh, imght[54,:])
chi2_distance(imgh, imght[55,:])
chi2_distance(imgh, imght[56,:])
chi2_distance(imgh, imght[57,:])
chi2_distance(imgh, imght[58,:])
chi2_distance(imgh, imght[59,:])
chi2_distance(imgh, imght[60,:])
chi2_distance(imgh, imght[61,:])
chi2_distance(imgh, imght[62,:])
chi2_distance(imgh, imght[63,:])
chi2_distance(imgh, imght[64,:])


chi_sfh=np.array([1408.7685050027608,742.37324080294366,3298.8785984177998,509.77318750300486,447.98469618697789,608.53294857065703,
765.71300357798452,2407.9704555837138,1572.6331865356726,1995.0150711918734,2115.6193207958918,2273.9658062483354,662.18427568214395,
1008.2421095215208,780.40271175550981,802.14406485862878,2341.6365050871964,2883.0782654885011,990.95219930923724,10469.904824242429,
504.21531910030768,496.72631072101927,778.23938895418758,687.51455575433818,437.5841666301589,1016.5822120144751,153.01790401110313,
610.35053564104362,831.83798012081991,1107.4072632366231,802.93316712523688,286.0817378541135,1272.38151605743,1141.1708062575901,
1452.8909070358307,363.47371055491874,1368.4992072572229,673.2623399080959,317.40075115448462,324.88208928081605,562.02964258915847,
1205.0491578406281,294.37870067555946,763.53586858567394,4737.085537105253,881.16720116890417,1098.2745712597127,861.46807904661796,
835.88647705821347,919.41142539365569,3009.6092381203571,1132.4922604358196,3368.0117044145104,785.20621189382905,1643.857429558971,
388.25294891907606,800.55966187691979,1324.5429238393051, 695.38863677639745,1524.1279867774956,2962.347935350851,665.51000016284979,
515.66096567015029,3877.1166153683312,662.07881985471045])


weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=30, normed=False,weights=weights, histtype='step', lw=3,label='Mass $11<\log (M/M_{\odot})<12$ , color $2<g-r<3$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()


