# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
plt.imshow(img)
img = hdus[0].data
plt.imshow(img)
plt.clf()
plt.imshow(img, origin = 'lower')
img.shape
img.min()
t38= img[:,38].cumsum()
plt.figure()
plt.plot(t38)
t38= img[:,38]
plt.plot(t38)
plt.plot(t38)
t38= img[:,38].cumsum()
plt.plot(t38)
t37= img[:,37].cumsum()
s37= t37+t38
plt.plot(s37)
t36= img[:,36].cumsum()
s36= t36+s37
plt.plot(s36)
t35= img[:,35].cumsum()
s35= t35+s36
plt.plot(s35)
t34= img[:,34].cumsum()
s34= t34+s35
plt.plot(s34)
t33= img[:,33].cumsum()
s33= t33+s34
plt.plot(s33)
t32= img[:,32].cumsum()
s32= t32+s33
plt.plot(s32)
t31= img[:,31].cumsum()
s31= t31+s32
plt.plot(s31)
t30= img[:,30].cumsum()
s30= t30+s31
plt.plot(s30)
t29= img[:,29].cumsum()
s29= t29+s30
plt.plot(s29)
t28= img[:,28].cumsum()
s28= t28+s29
plt.plot(s28)
t27= img[:,27].cumsum()
s27= t27+s28
plt.plot(s27)
t26= img[:,26].cumsum()
s26= t26+s27
plt.plot(s26)
t25= img[:,25].cumsum()
s25= t25+s26
plt.plot(s25)
t24= img[:,24].cumsum()
s24= t24+s25
plt.plot(s24)
t23= img[:,23].cumsum()
s23= t23+s24
plt.plot(s23)
t22= img[:,22].cumsum()
s22= t22+s23
plt.plot(s22)
t21= img[:,21].cumsum()
s21= t21+s22
t21= img[:,21].cumsum()
plt.plot(s21)
t20= img[:,20].cumsum()
s20= t20+s21
plt.plot(s20)
t19= img[:,19].cumsum()
s19= t19+s20
plt.plot(s19)
t18= img[:,18].cumsum()
s18= t18+s19
plt.plot(s18)
t17= img[:,17].cumsum()
s17= t17+s18
plt.plot(s17)
t17= img[:,17].cumsum()
plt.plot(t17)
t16= img[:,16].cumsum()
plt.plot(t16)
t15= img[:,15].cumsum()
plt.plot(t15)
s15= t15+s18
plt.plot(s15)
t14= img[:,14].cumsum()
s14= t14+s17
plt.plot(s14)
plt.plot(s15)
t13= img[:,13].cumsum()
s13= t13+s17
plt.plot(s13)
plt.plot(s15)
t12= img[:,12].cumsum()
plt.plot(t12)
plt.plot(s15)
s12= t12+s17
plt.plot(s12)
plt.plot(s17)
s13= t13+s15
plt.plot(s13)
s12= t12+ s13
plt.plot(s12)
t11= img[:,11].cumsum()
s11= t11+s12
plt.plot(s11)
