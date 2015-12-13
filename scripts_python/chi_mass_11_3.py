

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

hdus = fits.open('all_mass_11_3_SFH_Mass.fits.gz')
img = hdus[0].data
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)

hdus2 = fits.open('ARP118.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus3 = fits.open('IC0944.p_e.rad_SFH_lum_Mass.fits.gz')
img3 = hdus3[0].data
img_mask3= np.power(10, img3)

hdus4 = fits.open('IC0994.p_e.rad_SFH_lum_Mass.fits.gz')
img4 = hdus4[0].data
img_mask4= np.power(10, img4)

hdus5 = fits.open('IC1079.p_e.rad_SFH_lum_Mass.fits.gz')
img5 = hdus5[0].data
img_mask5= np.power(10, img5)

hdus6 = fits.open('IC1755.p_e.rad_SFH_lum_Mass.fits.gz')
img6 = hdus6[0].data
img_mask6= np.power(10, img6)

hdus7 = fits.open('IC3598.p_e.rad_SFH_lum_Mass.fits.gz')
img7 = hdus7[0].data
img_mask7= np.power(10, img7)

hdus8 = fits.open('NGC0155.p_e.rad_SFH_lum_Mass.fits.gz')
img8 = hdus8[0].data
img_mask8= np.power(10, img8)

hdus9 = fits.open('NGC0160.p_e.rad_SFH_lum_Mass.fits.gz')
img9 = hdus9[0].data
img_mask9= np.power(10, img9)

hdus10 = fits.open('NGC0171.p_e.rad_SFH_lum_Mass.fits.gz')
img10 = hdus10[0].data
img_mask10= np.power(10, img10)

hdus11 = fits.open('NGC0217.p_e.rad_SFH_lum_Mass.fits.gz')
img11 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus12 = fits.open('NGC0364.p_e.rad_SFH_lum_Mass.fits.gz')
img12 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus13 = fits.open('NGC0447.p_e.rad_SFH_lum_Mass.fits.gz')
img13 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus14 = fits.open('NGC0499.p_e.rad_SFH_lum_Mass.fits.gz')
img14= hdus2[0].data
img_mask2= np.power(10, img2)

hdus15 = fits.open('NGC0507.p_e.rad_SFH_lum_Mass.fits.gz')
img15 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus16 = fits.open('NGC0515.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus17 = fits.open('NGC0528.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus18 = fits.open('NGC0570.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus19 = fits.open('NGC0774.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus20 = fits.open('NGC0833.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus21 = fits.open('NGC0842.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus22 = fits.open('NGC0924.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus23 = fits.open('NGC1026.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus24 = fits.open('NGC1041.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus25 = fits.open('NGC1060.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus26 = fits.open('NGC1132.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus27 = fits.open('NGC1167.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus28 = fits.open('NGC1349.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus29 = fits.open('NGC2486.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus30 = fits.open('NGC2487.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus31 = fits.open('NGC2507.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus32 = fits.open('NGC2522.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus33 = fits.open('NGC2554.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus34 = fits.open('NGC3106.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus35 = fits.open('NGC3303.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus36 = fits.open('NGC4816.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus37 = fits.open('NGC4874.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus38 = fits.open('NGC4956.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus39 = fits.open('NGC5157.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus40 = fits.open('NGC5267.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus41 = fits.open('NGC5513.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus42 = fits.open('NGC5532.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus43 = fits.open('NGC5549.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus44 = fits.open('NGC5598.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus45 = fits.open('NGC5614.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus46 = fits.open('NGC5739.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus47 = fits.open('NGC5797.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus48 = fits.open('NGC5908.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus49 = fits.open('NGC5928.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus50 = fits.open('NGC5987.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus51 = fits.open('NGC6023.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus52 = fits.open('NGC6081.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus53 = fits.open('NGC6146.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus54 = fits.open('NGC6150.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus55 = fits.open('NGC6314.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus56 = fits.open('NGC6338.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus57 = fits.open('NGC6977.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus58 = fits.open('NGC6978.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus59 = fits.open('NGC7025.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus60 = fits.open('NGC7236.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus61 = fits.open('NGC7311.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus62 = fits.open('NGC7550.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus63 = fits.open('NGC7563.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus64 = fits.open('NGC7671.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus65 = fits.open('NGC7684.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus66 = fits.open('NGC7711.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus67 = fits.open('NGC7722.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus68 = fits.open('NGC7738.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus69 = fits.open('NGC7782.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus70 = fits.open('NGC7824.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus71= fits.open('NGC6150B.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus72 = fits.open('UGC01062.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus73 = fits.open('UGC01271.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus74 = fits.open('UGC01274.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus75 = fits.open('UGC01749.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus76 = fits.open('UGC02018.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus77 = fits.open('UGC02229.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus78 = fits.open('UGC03151.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus79 = fits.open('UGC04136.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus80 = fits.open('UGC05113.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus81 = fits.open('UGC05771.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus82 = fits.open('UGC06036.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus83 = fits.open('UGC06312.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus84 = fits.open('UGC08322.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus85 = fits.open('UGC09492.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus86 = fits.open('UGC09629.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus87 = fits.open('UGC10097.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus88 = fits.open('UGC10380.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus89 = fits.open('UGC10695.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus90 = fits.open('UGC10905.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus100 = fits.open('UGC11228.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus101 = fits.open('UGC11694.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus102 = fits.open('UGC11958.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus103 = fits.open('UGC12274.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus104 = fits.open('VIIZw700.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus105 = fits.open('NGC6166NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)








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

