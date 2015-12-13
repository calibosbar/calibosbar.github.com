
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



hdus = fits.open('/home/jeffrey/Documents/mass_10_1/all_mass_10_1_SFH_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/Documents/mass_11_3/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/Documents/mass_10_1/IC0159.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/IC0208.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/IC1151.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/IC1256.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/IC1528.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/IC2101.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/IC5309.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0197.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0237.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0444.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0496.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0551.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0716.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0873.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC0991.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC1056.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC1093.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC1656.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC2530.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC2540.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC2604.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC2623.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC2730.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC2805.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC2906.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC3381.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC3614.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC3815.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC3994.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC4047.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC4470.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC4711.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC4961.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5016.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5320.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5480.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5520.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5622.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5633.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5656.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5665.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5730.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5732.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5735.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5829.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5950.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5953.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC5980.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC6004.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC6063.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC6090.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC6132.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC6155.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC6186.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC7536.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC7549.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC7691.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/NGC7819.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC9837.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC00139.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC00148.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC00386.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC01057.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC01938.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC02134.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC02405.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC02628.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC02690.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC03539.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC03944.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC03969.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC04140.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC04195.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC04308.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC04375.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC04461.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC04730.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC05396.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC05598.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC06517.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC07145.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC08004.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC08250.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09067.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09253.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09262.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09291.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09476.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09542.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09665.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09708.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09777.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09842.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09873.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC09892.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC10257.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC10331.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC10384.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC11262.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC12224.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC12519.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC12688.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC12816.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC12857.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC12864.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/CGCG163-062.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/CGCG536-030.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/ESO539-G014.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/MCG-01-10-015.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/MCG-01-10-019.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/MCG-02-02-030.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/MCG-02-02-040.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/MCG-02-51-004.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_1/UGC11680NED02.p_e.rad_SFH_lum_Mass.fits.gz',]

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(114,1404)

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
chi2_distance(imgh, imght[65,:])
chi2_distance(imgh, imght[66,:])
chi2_distance(imgh, imght[67,:])
chi2_distance(imgh, imght[68,:])
chi2_distance(imgh, imght[69,:])
chi2_distance(imgh, imght[70,:])
chi2_distance(imgh, imght[71,:])
chi2_distance(imgh, imght[72,:])
chi2_distance(imgh, imght[73,:])
chi2_distance(imgh, imght[74,:])
chi2_distance(imgh, imght[75,:])
chi2_distance(imgh, imght[76,:])
chi2_distance(imgh, imght[77,:])
chi2_distance(imgh, imght[78,:])
chi2_distance(imgh, imght[79,:])
chi2_distance(imgh, imght[80,:])
chi2_distance(imgh, imght[81,:])
chi2_distance(imgh, imght[82,:])
chi2_distance(imgh, imght[83,:])
chi2_distance(imgh, imght[84,:])
chi2_distance(imgh, imght[85,:])
chi2_distance(imgh, imght[86,:])
chi2_distance(imgh, imght[87,:])
chi2_distance(imgh, imght[88,:])
chi2_distance(imgh, imght[89,:])
chi2_distance(imgh, imght[90,:])
chi2_distance(imgh, imght[91,:])
chi2_distance(imgh, imght[92,:])
chi2_distance(imgh, imght[93,:])
chi2_distance(imgh, imght[94,:])
chi2_distance(imgh, imght[95,:])
chi2_distance(imgh, imght[96,:])
chi2_distance(imgh, imght[97,:])
chi2_distance(imgh, imght[98,:])
chi2_distance(imgh, imght[99,:])
chi2_distance(imgh, imght[100,:])
chi2_distance(imgh, imght[101,:])
chi2_distance(imgh, imght[102,:])
chi2_distance(imgh, imght[103,:])
chi2_distance(imgh, imght[104,:])
chi2_distance(imgh, imght[105,:])
chi2_distance(imgh, imght[106,:])
chi2_distance(imgh, imght[107,:])
chi2_distance(imgh, imght[108,:])
chi2_distance(imgh, imght[109,:])
chi2_distance(imgh, imght[110,:])
chi2_distance(imgh, imght[111,:])
chi2_distance(imgh, imght[112,:])
chi2_distance(imgh, imght[113,:])

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


