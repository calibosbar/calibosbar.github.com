

from scipy.stats import chi2
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure, show, rc
from scipy.stats import chisquare
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
from astropy.io import fits
import pylab





hdus1 = fits.open('all_mass_8_1_SFH_Mass.fits.gz')
img1 = hdus1[0].data

hdus2 = fits.open('all_mass_9_1_SFH_Mass.fits.gz')
img2 = hdus2[0].data

hdus3 = fits.open('all_mass_9_2_SFH_Mass.fits.gz')
img3 = hdus3[0].data

hdus4 = fits.open('all_mass_9_3_SFH_Mass.fits.gz')
img4 = hdus4[0].data

hdus5 = fits.open('all_mass_10_1_SFH_Mass.fits.gz')
img5 = hdus5[0].data

hdus6 = fits.open('all_mass_10_2_SFH_Mass.fits.gz')
img6 = hdus6[0].data

hdus7 = fits.open('all_mass_10_3_SFH_Mass.fits.gz')
img7 = hdus7[0].data

hdus8 = fits.open('all_mass_11_1_SFH_Mass.fits.gz')
img8 = hdus8[0].data

hdus9 = fits.open('all_mass_11_2_SFH_Mass.fits.gz')
img9 = hdus9[0].data

hdus10 = fits.open('all_mass_11_3_SFH_Mass.fits.gz')
img10 = hdus10[0].data

hdus11 = fits.open('all_mass_11_4_SFH_Mass.fits.gz')
img11 = hdus11[0].data




fig, ((ax1,ax2,ax3,ax4), (ax5,ax6,ax7,ax8),(ax9,ax10,ax11,ax12),(ax13,ax14,ax15,ax16))= plt.subplots(4, 4, sharex='col', sharey='row')


ax13.axes.get_xaxis().set_ticks([])
ax14.axes.get_xaxis().set_ticks([])
ax15.axes.get_xaxis().set_ticks([])
ax16.axes.get_xaxis().set_ticks([])
ax13.axes.get_yaxis().set_ticks([])
ax9.axes.get_yaxis().set_ticks([])
ax5.axes.get_yaxis().set_ticks([])
ax1.axes.get_yaxis().set_ticks([])
ax13.set_xlabel('$8<\log (M/M_{\odot})<9$')
ax14.set_xlabel('$9<\log (M/M_{\odot})<10$')
ax15.set_xlabel('$10<\log (M/M_{\odot})<11$')
ax16.set_xlabel('$11<\log (M/M_{\odot})<12$')
ax13.set_ylabel('$1<g-r<2$')
ax9.set_ylabel('$2<g-r<3$')
ax5.set_ylabel('$3<g-r<4$')
ax1.set_ylabel('$4<g-r<5$')


x =np.linspace(0, 2, 36)



ax13.set_yscale('log')
diffa= (np.diff(img1, axis=1))/0.37837838

t37= diffa[:,37]
t36= diffa[:,36]
t35= diffa[:,35]
t34= diffa[:,34]
t33= diffa[:,33]
t32= diffa[:,32]
t31= diffa[:,31]
t30= diffa[:,30]
t29= diffa[:,29]
t28= diffa[:,28]
t27= diffa[:,27]
t26= diffa[:,26]
t25= diffa[:,25]
t24= diffa[:,24]
t23= diffa[:,23]
t22= diffa[:,22]
t21= diffa[:,21]
t20= diffa[:,20]
t19= diffa[:,19]
t18= diffa[:,18]
t17= diffa[:,17]
t16= diffa[:,16]
t15= diffa[:,15]
t14= diffa[:,14]
t13= diffa[:,13]
t12= diffa[:,12]
t11= diffa[:,11]
t10= diffa[:,10]
t09= diffa[:,9]
t08= diffa[:,8]
t07= diffa[:,7]
t06= diffa[:,6]
t05= diffa[:,5]
t04= diffa[:,4]
t03= diffa[:,3]
t02= diffa[:,2]

st1=(t37+t36+t35+t34+t33+t32+t31+t30+t29+t28+t27+t26+t25+t24+t23+t22+t21+t20+t19+t18+t17+t16+t15+t14+t13+t12+t11+t10+t09+t13+t12+t11+t10+t09+t08+t07+t06+t05+t04+t03+t02)/36
ax13.plot(x,st1, label= '8<Gyrs<4' )
#plt.legend()


ax14.set_yscale('log')
diffb= (np.diff(img2, axis=1))/0.37837838
a37= diffb[:,37]
a36= diffb[:,36]
a35= diffb[:,35]
a34= diffb[:,34]
a33= diffb[:,33]
a32= diffb[:,32]
a31= diffb[:,31]
a30= diffb[:,30]
a29= diffb[:,29]
a28= diffb[:,28]
a27= diffb[:,27]
a26= diffb[:,26]
a25= diffb[:,25]
a24= diffb[:,24]
a23= diffb[:,23]
a22= diffb[:,22]
a21= diffb[:,21]
a20= diffb[:,20]
a19= diffb[:,19]
a18= diffb[:,18]
a17= diffb[:,17]
a16= diffb[:,16]
a15= diffb[:,15]
a14= diffb[:,14]
a13= diffb[:,13]
a12= diffb[:,12]
a11= diffb[:,11]
a10= diffb[:,10]
a09= diffb[:,9]
a08= diffb[:,8]
a07= diffb[:,7]
a06= diffb[:,6]
a05= diffb[:,5]
a04= diffb[:,4]
a03= diffb[:,3]
a02= diffb[:,2]

ta1=(a37+a36+a35+a34+a33+a32+a31+a30+a29+a28+a27+a26+a25+a24+a23+a22+a21+a20+a19+a18+a17+a16+a15+a14+a13+a12+a11+a10+a09+a13+a12+a11+a10+a09+a08+a07+a06+a05+a04+a03+a02)/36
ax14.plot(x,ta1, label= '8<Gyrs<4' )


ax10.set_yscale('log')
diffb= (np.diff(img3, axis=1))/0.37837838
b37= diffb[:,37]
b36= diffb[:,36]
b35= diffb[:,35]
b34= diffb[:,34]
b33= diffb[:,33]
b32= diffb[:,32]
b31= diffb[:,31]
b30= diffb[:,30]
b29= diffb[:,29]
b28= diffb[:,28]
b27= diffb[:,27]
b26= diffb[:,26]
b25= diffb[:,25]
b24= diffb[:,24]
b23= diffb[:,23]
b22= diffb[:,22]
b21= diffb[:,21]
b20= diffb[:,20]
b19= diffb[:,19]
b18= diffb[:,18]
b17= diffb[:,17]
b16= diffb[:,16]
b15= diffb[:,15]
b14= diffb[:,14]
b13= diffb[:,13]
b12= diffb[:,12]
b11= diffb[:,11]
b10= diffb[:,10]
b09= diffb[:,9]
b08= diffb[:,8]
b07= diffb[:,7]
b06= diffb[:,6]
b05= diffb[:,5]
b04= diffb[:,4]
b03= diffb[:,3]
b02= diffb[:,2]

ub1=(b37+b36+b35+b34+b33+b32+b31+b30+b29+b28+b27+b26+b25+b24+b23+b22+b21+b20+b19+b18+b17+b16+b15+b14+b13+b12+b11+b10+b09+b13+b12+b11+b10+b09+b08+b07+b06+b05+b04+b03+b02)/36
ax10.plot(x,ub1, label= '8<Gyrs<4' )
#plb.legend()




#pla.legend()




