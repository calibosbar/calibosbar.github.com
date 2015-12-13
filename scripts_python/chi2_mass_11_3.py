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



hdus = fits.open('/home/jeffrey/Documents/mass_11_3/all_mass_11_3_SFH_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/Documents/mass_11_3/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/Documents/mass_11_3/ARP118.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/IC0994.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/IC1079.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/IC1755.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/IC3598.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0155.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0160.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0171.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0217.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0364.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0447.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0499.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0507.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0515.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0528.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0570.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0774.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0833.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0842.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC0924.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC1026.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC1041.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC1060.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC1132.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC1167.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC1349.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC2486.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC2487.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC2507.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC2522.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC2554.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC3106.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC3303.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC4816.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC4874.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC4956.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5157.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5267.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5513.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5532.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5549.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5598.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5614.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5739.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5797.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5908.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5928.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC5987.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6023.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6081.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6146.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6150.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6314.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6338.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6977.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6978.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7025.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7236.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7311.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7550.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7563.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7671.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7684.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7711.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7722.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7738.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7782.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC7824.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6150B.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC01062.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC01271.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC01274.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC01749.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC02018.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC02229.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC03151.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC04136.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC05113.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC05771.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC06036.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC06312.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC08322.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC09492.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC09629.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC10097.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC10380.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC10695.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC10905.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC11228.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC11694.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC11958.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/UGC12274.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/VIIZw700.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_11_3/NGC6166NED01.p_e.rad_SFH_lum_Mass.fits.gz',]

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(94,1404)

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


chi_sfh=np.array([3337.5338616435602,482.32851460998756,1665.5876806016786,1031.6177142692582,1142.4705146765598,759.66900297487359,1260.2054019795296,
1535.6552639346291,3085.8014371812692,581.45847757945364,1421.5321717656993,2682.0681723747739,1976.5406738256913,1579.4350420199719,3542.1520906793471,1350.7002407837742,
1205.3242134287375,1540.4932984147561,992.47788832493688,1426.9466770436979,1051.1893073060214,1809.2006293205536,2069.4488702225462,1926.2640431556938,1057.5643715406786,
606.67294110328248,1876.1596351743501,2614.8379825149259,429.4107627057906,4744.5751493440712,1219.8741668464447,1257.7411849072739,1143.6265270554277,2324.9205280500623,
2835.7180221711651,2438.0576446798882,1134.0479971615489,864.21778915149787,407.71110503923455,958.79731534262112,982.18626886218203,510.12548650028867,881.60724246862333,
539.51538245537552,1744.1115969877219,2844.6746355288178,1248.040795460146,701.69854193125911,3629.4242533194829,711.87486283876456,1255.7293461543604,1306.7677556033605,
775.61062946004631,2601.9736624448442,906.73557663111183,654.03315585407063,1610.4628813550621,1028.6293957406469,7096.0466289239948,1008.004686537206,1450.5753716327069,
5620.1545666518177,5279.4287483440394,1123.0323831187904,1076.6617373127262,1536.1679894376657,558.5356214534695,925.91706295879828,2479.6974087151966,3536.8457430330432,
731.48650564743014,4740.7108168141776,1978.7995451555605,614.1875374343291,807.28376595412362,612.5407299331855,1811.4610626328713,2795.8917945571607,580.48101692904152,
5433.7541395662829,2339.4364845925261,3473.6502467993446,815.39502585751165,1374.3488306170473,1765.8749621512288,4424.4733527323742,1472.9931358896677,498.36004942627835,
798.329257797129,2869.226730889377,2687.674411268018,867.84664273107728,3566.5688896878041,3390.1422908711279])

weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=30, normed=False,weights=weights, histtype='step', lw=3,label='Mass $11<\log (M/M_{\odot})<12$ , color $3<g-r<4$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()











