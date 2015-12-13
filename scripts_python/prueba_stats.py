
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
from scipy.stats import chisquare
import pylab

#plt.xkcd()
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img1 = hdus[0].data

fig=plt.figure()
ax1 = fig.add_subplot(111)
img_mask= np.power(10, img1)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0
ticksy=np.linspace(0,2,9)
ticksx=np.linspace(0,13,9)
ax1.set_xticklabels(ticksx)
ax1.set_yticklabels(ticksy)
plt.imshow(img_masked, origin = 'lower',aspect='auto') 
plt.pcolor(img_masked, norm=LogNorm())
plt.set_cmap('seismic')
cbar=plt.colorbar() #ticks=[0,1,2,3,4,5]
#cbar.ax.set_yticklabels(np.array([0.01,0.1,1,10,100]))
cbar.ax.set_ylabel('Log $\Sigma_{*}$ $[M_{sun}/pc^{2}]$', rotation=270, fontsize=10, verticalalignment='top')
plt.ylabel('$R/R_e$')
plt.xlabel('Look Back in time (Gigayears)')
plt.title('UGC11680NED01') 

hdus = fits.open('all_agns_SFH_Mass.fits.gz')
img2 = hdus[0].data

fig=plt.figure()
ax1 = fig.add_subplot(111)
img_agn= img2
where_are_NaNs = isnan(img_agn)
img_agn[where_are_NaNs] = 0
ticksy=np.linspace(0,2,9)
ticksx=np.linspace(0,13,9)
ax1.set_xticklabels(ticksx)
ax1.set_yticklabels(ticksy)
plt.imshow(img_agn, origin = 'lower',aspect='auto') 
plt.pcolor(img_agn, norm=LogNorm())
plt.set_cmap('seismic')
cbar=plt.colorbar() #ticks=[0,1,2,3,4,5]
#cbar.ax.set_yticklabels(np.array([0.01,0.1,1,10,100]))
cbar.ax.set_ylabel('Log $\Sigma_{*}$ $[M_{sun}/pc^{2}]$', rotation=270, fontsize=10, verticalalignment='top')
plt.ylabel('$R/R_e$')
plt.xlabel('Look Back in time (Gigayears)')
plt.title('All AGNs') 

mean_masked=np.mean(img_masked)
std_masked=np.std(img_masked)
mean_agn=np.mean(img_agn)
std_agn=np.std(img_agn)

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

t38= img_masked[:,38]
t37= img_masked[:,37]
t36= img_masked[:,36]
t35= img_masked[:,35]
t34= img_masked[:,34]
t33= img_masked[:,33]
t32= img_masked[:,32]
t31= img_masked[:,31]
t30= img_masked[:,30]
t29= img_masked[:,29]
t28= img_masked[:,28]
t27= img_masked[:,27]
t26= img_masked[:,26]
t25= img_masked[:,25]
t24= img_masked[:,24]
t23= img_masked[:,23]
t22= img_masked[:,22]
t21= img_masked[:,21]
t20= img_masked[:,20]
t19= img_masked[:,19]
t18= img_masked[:,18]
t17= img_masked[:,17]
t16= img_masked[:,16]
t15= img_masked[:,15]
t14= img_masked[:,14]
t13= img_masked[:,13]
t12= img_masked[:,12]
t11= img_masked[:,11]
t10= img_masked[:,10]
t9= img_masked[:,9]
t8= img_masked[:,8]
t7= img_masked[:,7]
t6= img_masked[:,6]
t5= img_masked[:,5]
t4= img_masked[:,4]
t3= img_masked[:,3]
t2= img_masked[:,2]
t1= img_masked[:,1]
t0= img_masked[:,0]

t0= img_masked[0,:]
t1= img_masked[1,:]
t2= img_masked[2,:]
t3= img_masked[3,:]
t4= img_masked[4,:]
t5= img_masked[5,:]
t6= img_masked[6,:]
t7= img_masked[7,:]
ts1= (t0+t1+t2+t3+t4+t5+t6+t7)/8

a0= img_agn[0,:]
a1= img_agn[1,:]
a2= img_agn[2,:]
a3= img_agn[3,:]
a4= img_agn[4,:]
a5= img_agn[5,:]
a6= img_agn[6,:]
a7= img_agn[7,:]
as1= (a0+a1+a2+a3+a4+a5+a6+a7)/8

bins=np.arange(0,14)
indt = np.digitize(img_masked, bins)
inda = np.digitize(img_agn, bins)

for n in range(ts1.size):
	print bins[indt[n]-1], "<=", ts1[n], "<", bins[indt[n]]

#ts1.sort()
tmean = np.mean(ts1)
tstd = np.std(ts1)
pdf = stats.norm.pdf(ts1, tmean, tstd)
plt.plot(ts1, pdf)


#as1.sort()
amean = np.mean(as1)
astd = np.std(as1)
pdf = stats.norm.pdf(as1, amean, astd)
plt.plot(as1, pdf)


np.mean(t38)
np.std(t38)

mu, kappa = 70.377539192790579, 152.339097359177
import scipy.special as sps
count, bins, ignored = plt.hist(t38,38)
x = np.arange(-2, 2, 38)
y = -np.exp(kappa*np.cos(x-mu))/(2*np.pi*sps.jn(0,kappa))
plt.plot(x, y/max(y), linewidth=2, color='r')
plt.show()




t38.sort()
hmean = np.mean(t38)
hstd = np.std(t38)
pdf = stats.norm.pdf(t38, hmean, hstd)
plt.plot(t38, pdf)











