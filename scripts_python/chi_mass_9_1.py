
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



hdus = fits.open('/home/jeffrey/Documents/AGNS/all_mass_9_1_SFH_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/Documents/mass_11_3/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/Documents/mass_9_1/IC758.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/IC0776.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/IC0995.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/IC2604.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/IC3631.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC0014.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC0216.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC0755.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC0941.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC1677.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC2480.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC3057.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC3220.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC3395.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC3600.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC3913.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC4630.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5402.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5425.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5439.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5630.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5682.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5731.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5951.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC5954.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC6168.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/NGC7800.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC6930.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC9663.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGCA021.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC00312.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC00809.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC01014.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC02443.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC03899.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC04054.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC04258.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC04659.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC05187.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC05244.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC05326.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC05358.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC05377.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC05990.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC07012.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC07129.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC08231.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC08662.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC08733.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09056.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09071.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09080.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09356.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09448.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09849.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09901.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC09919.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC10796.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC12054.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC12308.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/UGC12494.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_9_1/MCG-01-54-016.p_e.rad_SFH_lum_Mass.fits.gz']

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(62,1404)

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

chi_sfh=np.array([987.70656759827409,701.34150243154431,389.49233343926437,98.575210745184833,223.61116743822788,281.33117087370999,509.34447514571258,
551.52259827133355,552.57815793380212,600.36053156113508,648.15015361698931,455.67490516288251,701.18885833398713,625.5387511745497,1165.4132015895016,
1723.4834391709721,483.5281504065652,2285.7749545137472,952.47412815591667,208.14802087058862,628.27879499637038,687.08023339503325,1906.8307791929012,
1241.1338465998351,3143.5889856780955,482.64935266456132,710.65150030346979,1225.0045247087303,4662.4964099949775,1496.0706655456122,576.67822955364454,
282.68167988269278,551.35438240395013,440.17937628074401,294.70734581518246,664.85928105609571,2276.0876847571785,111.42045062684967,2906.1443543091086,
3423.6716289107835,499.79682963358186,411.56041474956379,582.19511876668082,233.53703564678767,884.25571754550697,751.65513050897448,4186.2793067272314,
1929.9291099582674,239.08941419645308,321.56077984788857,907.29977202874079,176.84766623361691,586.59105518578167,3224.3537301483607,213.98408701908781,
173.12580470081946,1079.8251650801674,384.45707404352095,1754.9402190141755,633.81333074988447,716.40343228326503,655.88097121835256,165.10963901712532,
218.82204100887481,55.67548717020128,2640.7642065544765,329.34000024295898,811.89414524051665,1675.0500206655211,462.49066293881873,1664.6401571722347,
1925.9187596495317,288.51447214681303,309.02001337988446,735.27344451990791,392.52432582560863,2725.6153028996869,494.59395899001044,373.45374223565301,
537.59186814080692,242.52445296793783,1526.5797039649121,955.21932151451244,291.2583498350333,878.05631588094207,899.60389348121271,765.0617841438941,
318.69583094051359,263.37372267650841,699.19362431267734,848.02014698106575,344.00889345526502,371.05976015502722,1213.1301317327025,202.71807413198295,
885.85921503036502,484.35142019016308,909.80209146235961,931.8525855045782,874.59417477021873,362.17061005862666,913.25046449622891,865.51542078087584,
280.617082100446,743.32661321620458,1030.3994625248936,905.12273988768015,1869.2606235098181,531.59404008710226,1052.3396438960956,2183.6929730521229,
1034.9566247792231,1247.8185322479417,1573.7435963632029])



weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=50, normed=False,weights=weights, histtype='step', lw=3,label='Mass $10<\log (M/M_{\odot})<11$ , color $1<g-r<2$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()


