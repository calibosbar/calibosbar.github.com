
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
import pylab






hdus = fits.open('/home/jeffrey/Documents/dataproducts/NGC4047.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
time=np.linspace(0,14,38)
time1=np.diff(time)
diff= np.diff(img_mask, axis=1)
diff1= diff/0.37837838
img_masked= diff1
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

plt.figure()
num_plots = 6
colormap = plt.cm.brg
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])

x =np.linspace(0, 2, 36)
plt.ylabel(' $Log <SFR> (M_{sun} yr^{-1} kpc^{-2})$ ')
plt.xlabel('$R/R_e$')
plt.yscale('log')
#plt.xscale('log')
plt.title('All AGNs')

t37= img_masked[:,37]
t36= img_masked[:,36]
t35= img_masked[:,35]
t34= img_masked[:,34]
st4=(t37+t36+t35+t34)/4
plt.plot(x,st4, label= '14<Gyrs<12.5' )
plt.legend()


t33= img_masked[:,33]
plt.plot(x,t33, label= '11.5 Gyrs' )
plt.legend()

t32= img_masked[:,32]
plt.plot(x,t32, label= '8 Gyrs' )
plt.legend()


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

st1=(t31+t30+t29+t28+t27+t26+t25+t24+t23+t22+t21+20+t19+t18)/14
plt.plot(x,st1, label= '8<Gyrs<4' )
plt.legend()


t18= img_masked[:,18]
t17= img_masked[:,17]
t16= img_masked[:,16]
t15= img_masked[:,15]
t14= img_masked[:,14]
t13= img_masked[:,13]
t12= img_masked[:,12]
t11= img_masked[:,11]
t10= img_masked[:,10]
t09= img_masked[:,9]
st3=(t18+t17+t16+t15+t14+t13+t12+t11+t10+t09)/10
plt.plot(x,st3, label= '4< Gyrs<1' )
plt.legend()


t13= img_masked[:,13]
t12= img_masked[:,12]
t11= img_masked[:,11]
t10= img_masked[:,10]
t09= img_masked[:,9]
t08= img_masked[:,8]
t07= img_masked[:,7]
t06= img_masked[:,6]
t05= img_masked[:,5]
t04= img_masked[:,4]
t03= img_masked[:,3]
t02= img_masked[:,2]
st4= (t13+t12+t11+t10+t09+t08+t07+t06+t05+t04+t03+t02)/12
plt.plot(x,st4, label= '1<Gyrs<Actual' )
plt.legend()



outfile = 'NGC4047.p_e.rad_SFR_lum_Mass.fits.gz'

hdu = fits.PrimaryHDU(img_masked)
hdu.writeto(outfile, clobber=True)















