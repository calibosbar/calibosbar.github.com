

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
#import seaborn as sn





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
#ax14.axes.get_xaxis().set_ticks([])
ax15.axes.get_xaxis().set_ticks([])
ax16.axes.get_xaxis().set_ticks([])
#ax13.axes.get_yaxis().set_ticks([])
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

ax10.set_ylabel('$M(t)/M_0$')
ax10.set_xlabel('Look back in time (Gyrs)')
ticksy=np.linspace(0,1,9)
ticksx=np.array([0,3.5,7,10,12,14])
ax14.set_xticklabels(ticksx)
ax13.set_yticklabels(ticksy)



x =np.linspace(0, 14, 39)

t1= img1[3,:]
t1[:] = t1[::-1].cumsum()
t1[:] = t1[::-1]
#s1=t1
s1= t1/np.max(t1)
ax13.plot(x ,s1, label= '0.0 < R/Re < 0.5' )
#ax13.legend(loc=3)

t20= img1[20,:]
t20[:] = t20[::-1].cumsum()
t20[:] = t20[::-1]
s20= t20/np.max(t20)
ax13.plot(x,s20, label= '0.5 < R/Re < 1.0' )
#ax13.legend(loc=3)


t31= img1[31,:]
t31[:] = t31[::-1].cumsum()
t31[:] = t31[::-1]
s31= t31/np.max(t31)
ax13.plot(x,s31, label= '1 < R/Re < 2' )
#ax13.legend(loc=3)


a1= img2[3,:]
a1[:] = a1[::-1].cumsum()
a1[:] = a1[::-1]
#s1=t1
r1= a1/np.max(a1)
ax14.plot(x ,s1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

a20= img2[20,:]
a20[:] = a20[::-1].cumsum()
a20[:] = a20[::-1]
r20= a20/np.max(a20)
ax14.plot(x,r20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


a31= img2[31,:]
a31[:] = a31[::-1].cumsum()
a31[:] = a31[::-1]
r31= a31/np.max(a31)
ax14.plot(x,r31, label= '1 < R/Re < 2' )
#plt.legend(loc=3)

b1= img3[3,:]
b1[:] = b1[::-1].cumsum()
b1[:] = b1[::-1]
#s1=t1
s1= b1/np.max(b1)
ax10.plot(x ,r1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

b20= img3[20,:]
b20[:] = b20[::-1].cumsum()
b20[:] = b20[::-1]
s20= b20/np.max(b20)
ax10.plot(x,s20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


b31= img3[31,:]
b31[:] = b31[::-1].cumsum()
b31[:] = b31[::-1]
s31= b31/np.max(b31)
ax10.plot(x,s31, label= '1 < R/Re < 2' )
#plt.legend(loc=3)

c1= img4[3,:]
c1[:] = c1[::-1].cumsum()
c1[:] = c1[::-1]
#s1=t1
u1= c1/np.max(c1)
ax6.plot(x ,u1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

c20= img4[20,:]
c20[:] = c20[::-1].cumsum()
c20[:] = c20[::-1]
u20= c20/np.max(c20)
ax6.plot(x,u20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


c31= img4[31,:]
c31[:] = c31[::-1].cumsum()
c31[:] = c31[::-1]
u31= c31/np.max(c31)
ax6.plot(x,u31, label= '1 < R/Re < 2' )
#plt.legend(loc=3)

d1= img5[3,:]
d1[:] = d1[::-1].cumsum()
d1[:] = d1[::-1]
#s1=t1
v1= d1/np.max(d1)
ax15.plot(x ,v1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

d20= img5[20,:]
d20[:] = d20[::-1].cumsum()
d20[:] = d20[::-1]
v20= d20/np.max(d20)
ax15.plot(x,v20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


d31= img5[31,:]
d31[:] = d31[::-1].cumsum()
d31[:] = d31[::-1]
v31= d31/np.max(d31)
ax15.plot(x,v31, label= '1 < R/Re < 2' )
#plt.legend(loc=3)

e1= img6[3,:]
e1[:] = e1[::-1].cumsum()
e1[:] = e1[::-1]
#s1=t1
w1= e1/np.max(e1)
ax11.plot(x ,w1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

e20= img6[20,:]
e20[:] = e20[::-1].cumsum()
e20[:] = e20[::-1]
w20= e20/np.max(e20)
ax11.plot(x,w20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


e31= img6[31,:]
e31[:] = e31[::-1].cumsum()
e31[:] = e31[::-1]
w31= e31/np.max(e31)
ax11.plot(x,w31, label= '1 < R/Re < 2' )
#plt.legend(loc=3)

f1= img7[3,:]
f1[:] = f1[::-1].cumsum()
f1[:] = f1[::-1]
#s1=t1
x1= f1/np.max(f1)
ax7.plot(x ,x1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

f20= img7[20,:]
f20[:] = f20[::-1].cumsum()
f20[:] = f20[::-1]
x20= f20/np.max(f20)
ax7.plot(x,x20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


f31= img7[31,:]
f31[:] = f31[::-1].cumsum()
f31[:] = f31[::-1]
x31= f31/np.max(f31)
ax7.plot(x,x31, label= '1 < R/Re < 2' )
#plt.legend(loc=3)

g1= img8[3,:]
g1[:] = g1[::-1].cumsum()
g1[:] = g1[::-1]
#s1=t1
y1= g1/np.max(g1)
ax16.plot(x ,y1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

g20= img8[20,:]
g20[:] = g20[::-1].cumsum()
g20[:] = g20[::-1]
y20= g20/np.max(g20)
ax16.plot(x,y20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


g31= img8[31,:]
g31[:] = g31[::-1].cumsum()
g31[:] = g31[::-1]
y31= g31/np.max(g31)
ax16.plot(x,y31, label= '1 < R/Re < 2' )
#plt.legend(loc=3

h1= img9[3,:]
h1[:] = h1[::-1].cumsum()
h1[:] = h1[::-1]
#s1=t1
z1= h1/np.max(h1)
ax12.plot(x ,z1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

h20= img9[20,:]
h20[:] = h20[::-1].cumsum()
h20[:] = h20[::-1]
z20= h20/np.max(h20)
ax12.plot(x,z20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


h31= img9[31,:]
h31[:] = h31[::-1].cumsum()
h31[:] = h31[::-1]
z31= h31/np.max(h31)
ax12.plot(x,z31, label= '1 < R/Re < 2' )
#plt.legend(loc=3

i1= img10[3,:]
i1[:] = i1[::-1].cumsum()
i1[:] = i1[::-1]
#s1=t1
p1= i1/np.max(i1)
ax8.plot(x ,p1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

i20= img10[20,:]
i20[:] = i20[::-1].cumsum()
i20[:] = i20[::-1]
p20= i20/np.max(i20)
ax8.plot(x,p20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


i31= img10[31,:]
i31[:] = i31[::-1].cumsum()
i31[:] = i31[::-1]
p31= i31/np.max(i31)
ax8.plot(x,p31, label= '1 < R/Re < 2' )
#plt.legend(loc=3

j1= img11[3,:]
j1[:] = j1[::-1].cumsum()
j1[:] = j1[::-1]
#s1=t1
q1= j1/np.max(j1)
ax4.plot(x ,q1, label= '0.0 < R/Re < 0.5' )
#plt.legend(loc=3)

j20= img11[20,:]
j20[:] = j20[::-1].cumsum()
j20[:] = j20[::-1]
q20= j20/np.max(j20)
ax4.plot(x,q20, label= '0.5 < R/Re < 1.0' )
#plt.legend(loc=3)


j31= img11[31,:]
j31[:] = j31[::-1].cumsum()
j31[:] = j31[::-1]
q31= j31/np.max(j31)
ax4.plot(x,q31, label= '1 < R/Re < 2' )
#plt.legend(loc=3

#ax16.legend(loc=9)
