
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

hdus1 = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus1[0].data
img_mask1= np.power(10, img1)

hdus2 = fits.open('IC0540.p_e.rad_SFH_lum_Mass.fits.gz')
img2 = hdus2[0].data
img_mask2= np.power(10, img2)

hdus3 = fits.open('NGC0426.p_e.rad_SFH_lum_Mass.fits.gz')
img3 = hdus3[0].data
img_mask3= np.power(10, img3)

hdus4 = fits.open('NGC0833.p_e.rad_SFH_lum_Mass.fits.gz')
img4 = hdus4[0].data
img_mask4= np.power(10, img4)

hdus5 = fits.open('NGC1667.p_e.rad_SFH_lum_Mass.fits.gz')
img5 = hdus5[0].data
img_mask5= np.power(10, img5)

hdus6 = fits.open('NGC2410.p_e.rad_SFH_lum_Mass.fits.gz')
img6 = hdus6[0].data
img_mask6= np.power(10, img6)

hdus7 = fits.open('NGC2623.p_e.rad_SFH_lum_Mass.fits.gz')
img7 = hdus7[0].data
img_mask7= np.power(10, img7)

hdus8= fits.open('NGC2639.p_e.rad_SFH_lum_Mass.fits.gz')
img8 = hdus8[0].data
img_mask8= np.power(10, img8)

hdus9 = fits.open('NGC3106.p_e.rad_SFH_lum_Mass.fits.gz')
img9= hdus9[0].data
img_mask9= np.power(10, img9)

hdus10 = fits.open('NGC3160.p_e.rad_SFH_lum_Mass.fits.gz')
img10 = hdus10[0].data
img_mask10= np.power(10, img10)

hdus11 = fits.open('NGC3303.p_e.rad_SFH_lum_Mass.fits.gz')
img11= hdus11[0].data
img_mask11= np.power(10, img11)

hdus12 = fits.open('NGC5635.p_e.rad_SFH_lum_Mass.fits.gz')
img12 = hdus12[0].data
img_mask12= np.power(10, img12)

hdus13 = fits.open('NGC5675.p_e.rad_SFH_lum_Mass.fits.gz')
img13 = hdus13[0].data
img_mask13= np.power(10, img13)

hdus14 = fits.open('NGC5739.p_e.rad_SFH_lum_Mass.fits.gz')
img14 = hdus14[0].data
img_mask14= np.power(10, img14)

hdus15 = fits.open('NGC6323.p_e.rad_SFH_lum_Mass.fits.gz')
img15 = hdus15[0].data
img_mask15= np.power(10, img15)

hdus16 = fits.open('NGC6394.p_e.rad_SFH_lum_Mass.fits.gz')
img16 = hdus16[0].data
img_mask16= np.power(10, img16)

hdus17 = fits.open('NGC7466.p_e.rad_SFH_lum_Mass.fits.gz')
img17 = hdus17[0].data
img_mask17= np.power(10, img17)

hdus18 = fits.open('NGC7738.p_e.rad_SFH_lum_Mass.fits.gz')
img18 = hdus18[0].data
img_mask18= np.power(10, img18)

hdus19 = fits.open('UGC00005.p_e.rad_SFH_lum_Mass.fits.gz')
img19 = hdus19[0].data
img_mask19= np.power(10, img19)

hdus20 = fits.open('UGC00987.p_e.rad_SFH_lum_Mass.fits.gz')
img20 = hdus20[0].data
img_mask20= np.power(10, img20)

hdus21 = fits.open('UGC03995.p_e.rad_SFH_lum_Mass.fits.gz')
img21 = hdus21[0].data
img_mask21= np.power(10, img21)

hdus22 = fits.open('UGC05771.p_e.rad_SFH_lum_Mass.fits.gz')
img22 = hdus22[0].data
img_mask22= np.power(10, img22)

hdus23 = fits.open('UGC09711.p_e.rad_SFH_lum_Mass.fits.gz')
img23 = hdus23[0].data
img_mask23= np.power(10, img23)

hdus24 = fits.open('MCG-02-02-030.p_e.rad_SFH_lum_Mass.fits.gz')
img24 = hdus24[0].data
img_mask24= np.power(10, img24)

nx, ny = img.shape
nx1,ny1 = img_mask1.shape
nx2,ny2 = img_mask2.shape
nx3,ny3 = img_mask3.shape
nx4,ny4 = img_mask4.shape
nx5,ny5 = img_mask5.shape
nx6,ny6 = img_mask6.shape
nx7,ny7 = img_mask7.shape
nx8,ny8 = img_mask8.shape
nx9,ny9 = img_mask9.shape
nx10,ny10 = img_mask10.shape
nx11,ny11 = img_mask11.shape
nx12,ny12 = img_mask12.shape
nx13,ny13 = img_mask13.shape
nx14,ny14 = img_mask14.shape
nx15,ny15 = img_mask15.shape
nx16,ny16 = img_mask16.shape
nx17,ny17 = img_mask17.shape
nx18,ny18 = img_mask18.shape
nx19,ny19 = img_mask19.shape
nx20,ny20 = img_mask20.shape
nx21,ny21 = img_mask21.shape
nx22,ny22 = img_mask22.shape
nx23,ny23 = img_mask23.shape
nx24,ny24 = img_mask24.shape

imgh = np.reshape(img, nx*ny)
imgh1 = np.reshape(img_mask1, nx1*ny1)
imgh2 = np.reshape(img_mask2, nx2*ny2)
imgh3 = np.reshape(img_mask3,nx3*ny3)
imgh4 = np.reshape(img_mask4,nx4*ny4)
imgh5 = np.reshape(img_mask5,nx5*ny5)
imgh6 = np.reshape(img_mask6,nx6*ny6)
imgh7 = np.reshape(img_mask7,nx7*ny7)
imgh8 = np.reshape(img_mask8,nx8*ny8)
imgh9 = np.reshape(img_mask9,nx9*ny9)
imgh10 = np.reshape(img_mask10,nx10*ny10)
imgh11 = np.reshape(img_mask11,nx11*ny11)
imgh12 = np.reshape(img_mask12,nx12*ny12)
imgh13 = np.reshape(img_mask13,nx13*ny13)
imgh14 = np.reshape(img_mask14,nx14*ny14)
imgh15 = np.reshape(img_mask15,nx15*ny15)
imgh16 = np.reshape(img_mask16,nx16*ny16)
imgh17 = np.reshape(img_mask17,nx17*ny17)
imgh18 = np.reshape(img_mask18,nx18*ny18)
imgh19 = np.reshape(img_mask19,nx19*ny19)
imgh20 = np.reshape(img_mask20,nx20*ny20)
imgh21 = np.reshape(img_mask21,nx21*ny21)
imgh22 = np.reshape(img_mask22,nx22*ny22)
imgh23 = np.reshape(img_mask23,nx23*ny23)
imgh24 = np.reshape(img_mask24,nx24*ny24)





def chi2_distance(histA, histB, eps = 1e-10):
	d = 0.5 * np.sum([((a - b) ** 2) / (a + b + eps)
		for (a, b) in zip(histA, histB)])
	return d

df=36*38

dof=chi2_distance(imgh,imgh1)
chi2_distance(imgh,imgh2)
chi2_distance(imgh,imgh3)
chi2_distance(imgh,imgh4)
chi2_distance(imgh,imgh5)
chi2_distance(imgh,imgh6)
chi2_distance(imgh,imgh7)
chi2_distance(imgh,imgh8)
chi2_distance(imgh,imgh9)
chi2_distance(imgh,imgh10)
chi2_distance(imgh,imgh11)
chi2_distance(imgh,imgh12)
chi2_distance(imgh,imgh13)
chi2_distance(imgh,imgh14)
chi2_distance(imgh,imgh15)
chi2_distance(imgh,imgh16)
chi2_distance(imgh,imgh17)
chi2_distance(imgh,imgh18)
chi2_distance(imgh,imgh19)
chi2_distance(imgh,imgh20)
chi2_distance(imgh,imgh21)
chi2_distance(imgh,imgh22)
chi2_distance(imgh,imgh23)
chi2_distance(imgh,imgh24)

chi_sfh=np.array([1163.7951354172355,1174.1648382707642,2101.6737019327734,2490.5802151243556,229.31686185335579,3702.9850404366912,8010.8976785589657,1426.7212688049581,731.80523340894786,633.8375833521842,668.22543662170722,461.59385929530822,1075.2450141562779,1038.3811265072629, 2069.479220105387, 1649.2137587619873,844.8898560351239,1402.7387981047777,2353.4017982222886,1152.548420343079,590.47101548905209,1849.3699518480601,1052.9411549059328])
#1187.7377355082149
weights = np.ones_like(chi_sfh)/float(len(chi_sfh))
fig, ax = plt.subplots(1, 1)
df = 1368
dof= 1187.7377355082149
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
xu = np.linspace(chi2.ppf(0.01, dof), chi2.ppf(0.99, dof), 100)
ax.axvline(x=dof/df,color='k', linestyle='dashed',lw=4, label='UGC11680NED01 $\chi^2$')
ax.hist(chi_sfh/df,bins=30, normed=False,weights=weights, histtype='step', lw=3,label='All AGNs $\chi^2$ reduced distribution')
ax.legend(loc='best', frameon=False)
ax.set_ylabel('Probability density $\chi ^2$')
ax.set_xlabel('$x$')
plt.show()

1+ 1./2. + 1./3. + 1./4. +1./5. + 1./6. +1./7. +1./8. + 1./9. + 1./10. + 1./11. +1./12. +1./13. + 1./14. + 1./15. +1./16. +1./17. +1./18. +1./19. +1./20. +1./21.+ 1./22. +1./23. +1./24.

3.7759581777535067

2(1368) = 2736


n=1368.9986218003824*(chi2.pdf(x, df))


2736+3.7759581777535067

varianza= 52.34286921995921
media= 1368


mass 8, color 1

2736+2.45
varianza= 52.33020160480943
media= 1368

mass 9 color 1
2736+5.11
varianza= 52.35561096959905

mass 9 color 2
2736+2.08
varianza=52.326666241984114


mass 9 color 3
2736+2.283
varianza=52.32860594359456

mass 10 color 1
2736+5.317
varianza= 52.3575877977586

mass 10 color2 
2736+ 3.815
varianza = 52.34324216171559

mass 10 color 3
2736+ 4.663
varianza= 52.35134191212294

