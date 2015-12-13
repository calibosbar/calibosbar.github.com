


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
from astropy.io import fits
import pylab

hdus = fits.open('/home/jeffrey/Documents/dataproducts/UGC11680NED01.p_e.rad_SFH_lum.fits.gz')
img = hdus[0].data
fig=plt.figure()
ax1 = fig.add_subplot(111)
#img_mask= np.power(10, img)
img_masked= img
where_are_NaNs = np.isnan(img_masked)
img_masked[where_are_NaNs] = 0
ticksy=np.linspace(0,2,9)
ticksx=np.linspace(0,10,9)
ax1.set_xticklabels(ticksx)
ax1.set_yticklabels(ticksy)
plt.imshow(img_masked, origin = 'lower', aspect='auto', interpolation='bicubic', vmin=0.0001, vmax= 1) 
plt.pcolor(img_masked, cmap=califa,) 
#plt.set_cmap('califa')
cbar=plt.colorbar() #ticks=[0,1,2,3,4,5]
#cbar.ax.set_yticklabels(np.array([0.01,0.1,1,10,100]))
cbar.ax.set_ylabel('Log $L_{*}$ $[L_{\odot} pc^{-2}]$', fontsize=10)
plt.ylabel('$R/R_e$')
plt.xlabel('Log time (yr)')
plt.title('UGC11680NED01') 
plt.show()


plt.figure()
x =np.linspace(0, 2, 36)
plt.ylabel('Log Surface Brigthness Intensity $[L_{sun}/kpc^{2}]$ ')
plt.xlabel('R/Re')
#plt.yscale('log')
#plt.xscale('log')
plt.title('UGC11680NED01 L(r) in the V-Band history radial')


t38= img_masked[:,38]
s38= t38/np.max(t38)
plt.plot(x,t38 , label= '10 Gyrs' )
plt.legend()


t37= img_masked[:,37]
t37n= t37+t38
s37= t37n/np.max(t37n)
#plt.plot(x,s37)

t36= img_masked[:,36]
t36n= t36+t37n
s36= t36n/np.max(t36n)
#plt.plot(x,s36)

t35= img_masked[:,35]
t35n= t35+t36n
s35= t35n/np.max(t35n)
plt.plot(x,t35n, label= '7.8 Gyrs')
plt.legend()


t34= img_masked[:,34]
t34n= t34+t35n
s34= t34n/np.max(t34n)
#plt.plot(x,s34)



t33= img_masked[:,33]
t33n= t33+t34n
s33= t33n/np.max(t33n)
#plt.plot(x,t33n, label= '7.8 Gyrs')
#plt.legend()

t32= img_masked[:,32]
t32n= t32+t33n
s32= t32n/np.max(t32n)
#plt.plot(x,s32)

t31= img_masked[:,31]
t31n= t31+t32n
s31= t31n/np.max(t31n)
#plt.plot(x,s31)

t30= img_masked[:,30]
t30n= t30+t31n
s30= t30n/np.max(t30n)
#plt.plot(x,t30n, label= '7.8 Gyrs')
#plt.legend()

t29= img_masked[:,29]
t29n= t29+t30n
s29= t29n/np.max(t29n)
#plt.plot(x,s29)

t28= img_masked[:,28]
t28n= t28+t29n
s28= t28n/np.max(t28n)
#plt.plot(x,s28)


t27= img_masked[:,27]
t27n= t27+t28n
s27= t27n/np.max(t27n)
#plt.plot(x,s27)

t26= img_masked[:,26]
t26n= t26+t27n
s26= t26n/np.max(t26n)
#plt.plot(x,s26)


t25= img_masked[:,25]
t25n= t25+t26n
s25= t25n/np.max(t25n)
#plt.plot(x,s25)

t24= img_masked[:,24]
t24n= t24+t25n
s24= t24n/np.max(t24n)
#plt.plot(x,s24)


t23= img_masked[:,23]
t23n= t23+t24n
s23= t23n/np.max(t23n)
#plt.plot(x,s23)

t22= img_masked[:,22]
t22n= t22+t23n
s22= t22n/np.max(t22n)
#plt.plot(x,s22)

t21= img_masked[:,21]
t21n= t21+t22n
s21= t21n/np.max(t21n)
#plt.plot(x,s21)

t20= img_masked[:,20]
t20n= t20+t21n
s20= t20n/np.max(t20n)
#plt.plot(x,s20)


t19= img_masked[:,19]
t19n= t19+t20n
s19= t19n/np.max(t19n)
#plt.plot(x,s19)


t18= img_masked[:,18]
t18n= t18+t19n
s18= t18n/np.max(t18n)
#plt.plot(x,s18)

t17= img_masked[:,17]
t17n= t17+t18n
s17= t17n/np.max(t17n)
#plt.plot(x,s17)


t16= img_masked[:,16]
t16n= t16+t17n
s16= t16n/np.max(t16n)
#plt.plot(x,s16)

t15= img_masked[:,15]
t15n= t15+t16n
s15= t15n/np.max(t15n)
#plt.plot(x,s15)

t14= img_masked[:,14]
t14n= t14+t15n
s14= t14n/np.max(t14n)
#plt.plot(x,s14)

t13= img_masked[:,13]
t13n= t13+t14n
s13= t13n/np.max(t13n)
#plt.plot(x,s13)

t12= img_masked[:,12]
t12n= t12+t13n
s12= t12n/np.max(t12n)
#plt.plot(x,s12)

t11= img_masked[:,11]
t11n= t11+t12n
s11= t11n/np.max(t11n)
#plt.plot(x,s11)

t10= img_masked[:,10]
t10n= t10+t11n
s10= t10n/np.max(t10n)
#plt.plot(x,s10)

t9= img_masked[:,9]
t9n= t9+t10n
s9= t9n/np.max(t9n)
#plt.plot(x,s9)

t8= img_masked[:,8]
t8n= t8+t9n
s8= t8n/np.max(t8n)
#plt.plot(x,s8)


t7= img_masked[:,7]
t7n= t7+t8n
s7= t7n/np.max(t7n)
#plt.plot(x,s7)

t6= img_masked[:,6]
t6n= t6+t7n
s6= t6n/np.max(t6n)
#plt.plot(x,s6)

t5= img_masked[:,5]
t5n= t5+t6n
s5= t5n/np.max(t5n)
#plt.plot(x,s5)

t4= img_masked[:,4]
t4n= t4+t5n
s4= t4n/np.max(t4n)
#plt.plot(x,s4)

t3= img_masked[:,3]
t3n= t3+t4n
s3= t3n/np.max(t3n)
#plt.plot(x,s3)

t2= img_masked[:,2]
t2n= t2+t3n
s2= t2n/np.max(t2n)
#plt.plot(x,s2)

t1= img_masked[:,1]
t1n= t1+t2n
s1= t1n/np.max(t1n)
#plt.plot(x,s1)

t0= img_masked[:,0]
t0n= t0+t1n
s0= t0n/np.max(t0n)
plt.plot(x,t0n, label= 'Actual')
plt.legend()



