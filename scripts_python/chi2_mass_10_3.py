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



hdus = fits.open('/home/jeffrey/Documents/mass_10_3/all_mass_10_3_SFH_Mass.fits.gz')
img = hdus[0].data
hdus1 = fits.open('/home/jeffrey/Documents/mass_11_3/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)



image_list = ['/home/jeffrey/Documents/mass_10_3/IC1652.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/IC2341.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/IC4215.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/IC4534.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0429.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0495.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0504.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0517.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0548.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0681.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0781.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0825.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0932.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC0955.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC1211.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC1666.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC2481.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC2553.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC2577.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC2880.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC3160.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC3300.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC3619.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5358.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5378.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5473.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5475.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5481.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5485.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5580.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5602.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5611.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5631.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5684.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5687.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5689.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5794.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5876.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC5935.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC6278.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC6427.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC6762.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC6945.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC7611.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/NGC7623.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC01123.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC02222.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC03960.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC08984.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC09539.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC09711.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC09937.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC10388.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC12518.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/UGC12653.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/CGCG429-012.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/SDSSJ015424.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/MCG-01-52-012.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/mass_10_3/MCG-02-08-014.p_e.rad_SFH_lum_Mass.fits.gz',]

image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= np.power(10, image_concat)
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

imght=img_concat_lin.reshape(59,1404)

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

chi_sfh=np.array([6554.4632793743131,1284.219024060897, 1856.8665038967686,636.96296298623577,1299.1325002105134,1580.7576903573281,
1662.9082881504246,1566.1409534301167,2270.3959853582628,2389.1803615714607,5410.8247268641353,1026.5255617672915,1255.2092589556105,
10926.840451035703,2295.0780400285389,1751.509362622402,24322.870163261825,1704.1712227013193, 1869.176666350646,2120.5180344505166,
1850.0741454074346,1020.1634825242691,2363.548981012525,1074.1625044190866,2361.6087077846755,2860.3347962972157,5308.8083334279054,
2404.5282521808067,1147.3797887929593,1480.6041039862916,3151.4666209599632,3943.4016736396984,2322.6651376937461,1668.3074034162944,
1388.171704256576,1461.4098965937512,1939.5471309571067,1616.2406232927001,2729.9228647384925,2640.3608235556412,3358.5434609436811,
728.50839049199499,8070.8725079093119,5109.6687217450972,1026.2468328583618,1457.2248604201527,1809.3160632562831,5015.1091396012598,
3322.0647640138723,2333.787989058635,3824.8737727714943,4638.0808883143618,1015.2765572158016,2434.5747736826211,5472.0183578092319,
1437.5049128337555,2903.0658720267484,1284.6521847658066,5320.6706635334922])


weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=35, normed=False,weights=weights, histtype='step', lw=3,label='Mass $10<\log (M/M_{\odot})<11$ , color $3<g-r<4$')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()











