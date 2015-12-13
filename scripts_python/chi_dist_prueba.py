
from matplotlib.pyplot import figure, show, rc
from kapteyn import kmpfit
from scipy.special import gammainc, chdtrc
from scipy.optimize import fminbound
from scipy.stats import chisquare
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
import pylab

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
imgh = np.reshape(img, nx*ny)
imgh1 = np.reshape(img_mask1, nx1*ny1)
imgh2 = np.reshape(img_mask2, nx2*ny2)
imgh3 = np.reshape(img_mask3,nx3*ny3)

def func(p, x):
   A, mu, sigma, zerolev = p
   return( A * np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + zerolev )

def residuals(p, data):
   x, y, err = data
   return (y-func(p,x))/ err


x = imgh
truepars = [741.01421963870212 ,9.3883286887501178,54.147681988464271, 0.0]
p0 = [741.01421963870212 ,9.3883286887501178,54.147681988464271, 0.0]
y = func(truepars, x) + imgh1
print "A max=", y.max()
N = len(x)
err = np.random.normal(9.3883286887501178, 54.147681988464271, N)

fitter = kmpfit.Fitter(residuals=residuals, data=(x,y,err))
fitter.parinfo = [{}, {}, {}, {'fixed':True}]  # Take zero level fixed in fit
fitter.fit(params0=p0)

if (fitter.status <= 0): 
   print "Status:  ", fitter.status
   print 'error message = ', fitter.errmsg
   raise SystemExit 

# Rescale the errors to force a reasonable result:
err[:] *= np.sqrt(1.187*fitter.rchi2_min)
fitter.fit()

print "======== Fit results =========="
print "Initial params:", fitter.params0
print "Params:        ", fitter.params
print "Iterations:    ", fitter.niter
print "Function ev:   ", fitter.nfev 
print "Uncertainties: ", fitter.xerror
print "dof:           ", fitter.dof
print "chi^2, rchi2:  ", fitter.chi2_min, fitter.rchi2_min
print "stderr:        ", fitter.stderr   
print "Status:        ", fitter.status

print "\n======== Statistics ========"

from scipy.stats import chi2
rv = chi2(fitter.dof)
print "Three methods to calculate the right tail cumulative probability:"
print "1. with gammainc(dof/2,chi2/2):  ", 1-gammainc(0.5*fitter.dof, 0.5*fitter.chi2_min)
print "2. with scipy's chdtrc(dof,chi2):", chdtrc(fitter.dof,fitter.chi2_min)
print "3. with scipy's chi2.cdf(chi2):  ", 1-rv.cdf(fitter.chi2_min)
print ""


xc = fitter.chi2_min
print "Threshold chi-squared at alpha=0.05: ", rv.ppf(1-0.05)
print "Threshold chi-squared at alpha=0.01: ", rv.ppf(1-0.01)

f = lambda x: -rv.pdf(x)
x_max = fminbound(f,1,200)
print """For %d degrees of freedom, the maximum probability in the distribution is
at chi-squared=%g """%(fitter.dof, x_max)

alpha = 0.05           # Select a p-value
chi2max = max(3*x_max, fitter.chi2_min)
chi2_threshold = rv.ppf(1-alpha)

print "For a p-value alpha=%g, we found a threshold chi-squared of %g"%(alpha, chi2_threshold)
print "The chi-squared of the fit was %g. Therefore: "%fitter.chi2_min 
if fitter.chi2_min <= chi2_threshold:
   print "we do NOT reject the hypothesis that the data is consistent with the model"
else:
   print "we REJECT the hypothesis that the data is consistent with the model"


# Plot the result
rc('legend', fontsize=8)
fig = figure(figsize=(7.2,9.5))
fig.subplots_adjust(hspace=0.5)
frame3 = fig.add_subplot(3,1,3)
xchi = np.linspace(1250, chi2max, 100)
ychi = rv.pdf(xchi)
delta = (xchi.max()-xchi.min())/40.0
frame3.plot(xchi, ychi, label="Degrees of freedom = %d"%(fitter.dof))
frame3.set_xlabel("$\chi^2$")
frame3.set_ylabel("$\mathrm{Probability}$")
frame3.set_title("$\chi^2 \mathrm{Probability\ density\ function\ for\, } \\nu=%d$"%fitter.dof)
frame3.plot((xc,xc),(0,ychi.max()), 'g', label="chi square (fit) = %g"%fitter.chi2_min)
frame3.plot((chi2_threshold,chi2_threshold),(0,ychi.max()), 'r', label="chi square threshold = %g"%chi2_threshold)

bbox_props = dict(boxstyle="larrow,pad=0.2", fc="cyan", ec="b", lw=2)
t = frame3.text(chi2_threshold-delta, ychi.max()/2, "Accept", ha="right", va="center", size=12,
               bbox=bbox_props)
bbox_props = dict(boxstyle="rarrow,pad=0.2", fc="red", ec="b", lw=2)
t = frame3.text(chi2_threshold+delta, ychi.max()/2, "Reject", ha="left", va="center", size=12,
               bbox=bbox_props)
leg = frame3.legend(loc=1)

ychi = rv.cdf(xchi)
frame2 = fig.add_subplot(3,1,2)
frame2.plot(xchi, ychi, label="Degrees of freedom = %d"%(fitter.dof))
frame2.set_xlabel("$\chi^2$")
frame2.set_ylabel("$\mathrm{Cumulative\ probability}$")
frame2.plot((xc,xc),(0,ychi.max()), 'g', label="chi square (fit) = %g"%fitter.chi2_min)
frame2.plot((chi2_threshold,chi2_threshold),(0,ychi.max()), 'r', label="chi square threshold = %g"%chi2_threshold)
frame2.plot((0,chi2_threshold),(1-alpha,1-alpha), 'r--', label="threshold for alpha = %g (=1-%g)"%(alpha,1-alpha))
frame2.set_title("$\chi^2 \mathrm{Cumulative\ distribution\ function}$")
leg = frame2.legend(loc=4)

frame1 = fig.add_subplot(3,1,1)
xd = np.linspace(x.min(), x.max(), 200)
label = "fit model: $y = A\ \exp\\left( \\frac{-(x-\mu)^2}{\sigma^2}\\right) + 0$"
frame1.plot(xd, func(fitter.params,xd), 'g', label=label)
frame1.errorbar(x, y, yerr=err,  fmt='bo', label="data")
frame1.set_xlabel("$x$")
frame1.set_ylabel("$y$")
vals = (fitter.chi2_min, fitter.rchi2_min, fitter.dof)
title = r"$\mathrm{Fit\ with\ } \chi^2=%g \mathrm{\ and\ } \chi^2_{\nu}=%g \,(\nu=%d)$"%vals
frame1.set_title(title, y=1.05)
leg = frame1.legend(loc=2)

show()


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
imgh21 = np.reshape(img_mask21,n21*ny21)
imgh22 = np.reshape(img_mask22,nx22*ny22)
imgh23 = np.reshape(img_mask23,nx23*ny23)
imgh24 = np.reshape(img_mask24,nx24*ny24)





def chi2_distance(histA, histB, eps = 1e-10):
	d = 0.5 * np.sum([((a - b) ** 2) / (a + b + eps)
		for (a, b) in zip(histA, histB)])
	return d

dof=36*38

chi2_distance(imgh,imgh1)/dof
chi2_distance(imgh1,img1)/dof


def chi2_distance(histA, histB, eps = 1e-10):
	d = 0.5 * np.sum([((a - b) ** 2) / (np.std(a) + eps)
		for (a, b) in zip(histA, histB)])
	return d

fig, ax = plt.subplots(1, 1)
