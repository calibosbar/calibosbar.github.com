
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
import pylab

#plt.xkcd()
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

plt.figure()
x =np.linspace(0, 14, 38)
plt.ylabel(' $SFR(t)$ ')
plt.xlabel('look back in time ($Gyrs$)')
#plt.yscale('log')
#plt.xscale('log')
plt.title('UGC11680NED01')


t0= img_masked[0,:]
t1= img_masked[1,:]
t2= img_masked[2,:]
t3= img_masked[3,:]
t4= img_masked[4,:]
t5= img_masked[5,:]
t6= img_masked[6,:]
t7= img_masked[7,:]

d0= np.diff(t0)/0.37837838
d1= np.diff(t1)/0.37837838
d2= np.diff(t2)/0.37837838
d3= np.diff(t3)/0.37837838
d4= np.diff(t4)/0.37837838
d5= np.diff(t5)/0.37837838
d6= np.diff(t6)/0.37837838
d7= np.diff(t7)/0.37837838
sd1= (d4+d5+d6+d7)/4

plt.plot(x,sd1, label= '0.0 < R/Re < 0.5 ')
plt.legend(loc=2)

t8= img_masked[8,:]
t9= img_masked[9,:]
t10= img_masked[10,:]
t11= img_masked[11,:]
t12= img_masked[12,:]
t13= img_masked[13,:]
t14= img_masked[14,:]
t15= img_masked[15,:]
t16= img_masked[16,:]
t17= img_masked[17,:]
t18= img_masked[18,:]
t19= img_masked[19,:]
t20= img_masked[20,:]
t21= img_masked[21,:]
t22= img_masked[22,:]
t23= img_masked[23,:]
t24= img_masked[24,:]
t25= img_masked[25,:]

d8= np.diff(t8)/0.37837838
d9= np.diff(t9)/0.37837838
d10= np.diff(t10)/0.37837838
d11= np.diff(t11)/0.37837838
d12= np.diff(t12)/0.37837838
d13= np.diff(t13)/0.37837838
d14= np.diff(t14)/0.37837838
d15= np.diff(t15)/0.37837838
d16= np.diff(t16)/0.37837838
d17= np.diff(t17)/0.37837838
d18= np.diff(t18)/0.37837838
d19= np.diff(t19)/0.37837838
d20= np.diff(t20)/0.37837838
d21= np.diff(t21)/0.37837838
d22= np.diff(t22)/0.37837838
d23= np.diff(t23)/0.37837838
d24= np.diff(t24)/0.37837838
d25= np.diff(t25)/0.37837838

sd2= (d8+d9+d10+d11+d12+d13+d14+d15+d16+d17+d18+d19+d20+d21+d22+d23+d24+d25)/18
plt.plot(x,sd2, label= '0.5 < R/Re < 1 ')
plt.legend(loc=2)

t26= img_masked[26,:]
t27= img_masked[27,:]
t28= img_masked[28,:]
t29= img_masked[29,:]
t30= img_masked[30,:]
t31= img_masked[31,:]
t32= img_masked[32,:]
t33= img_masked[33,:]
t34= img_masked[34,:]
t35= img_masked[35,:]

d26= np.diff(t20)/0.37837838
d27= np.diff(t21)/0.37837838
d28= np.diff(t22)/0.37837838
d29= np.diff(t23)/0.37837838
d30= np.diff(t24)/0.37837838
d31= np.diff(t25)/0.37837838
d32= np.diff(t20)/0.37837838
d33= np.diff(t21)/0.37837838
d34= np.diff(t22)/0.37837838
d35= np.diff(t23)/0.37837838

sd3= (d26+d27+d28+d29+d30+d31+d32+d33+d34+d35)/10
plt.plot(x,sd3, label= '1 < R/Re < 2 ')
plt.legend(loc=2)





