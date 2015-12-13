


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
import pylab

hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
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

x =np.linspace(0, 14, 38)
plt.ylabel(' $SFR (M_{sun} yr^{-1})$ ')
plt.xlabel('$R/R_e$')
#plt.yscale('log')
#plt.xscale('log')
plt.title('All AGNs')

t36= img_masked[35,:]
t35= img_masked[34,:]
t34= img_masked[33,:]
t33= img_masked[32,:]
t32= img_masked[31,:]
t31= img_masked[30,:]
t30= img_masked[29,:]
t29= img_masked[28,:]
t28= img_masked[27,:]
t27= img_masked[26,:]
t26= img_masked[25,:]
t25= img_masked[24,:]
t24= img_masked[23,:]
t23= img_masked[22,:]
t22= img_masked[21,:]
t21= img_masked[20,:]
t20= img_masked[19,:]
t19= img_masked[18,:]
t18= img_masked[17,:]
t17= img_masked[16,:]
t16= img_masked[15,:]
t15= img_masked[14,:]
t14= img_masked[13,:]
t13= img_masked[12,:]
t12= img_masked[11,:]
t11= img_masked[10,:]
t10= img_masked[9,:]
t9= img_masked[8,:]
t8= img_masked[7,:]
t7= img_masked[6,:]
t6= img_masked[5,:]
t5= img_masked[4,:]
t4= img_masked[3,:]
t3= img_masked[2,:]
t2= img_masked[1,:]
st4=(t36+t35+t34+t33+t32+t31+t30+t29+t28+t27+t26+t25+t24+t23+t22+t21+t20+t19+t18+t17+t16+t15+t14+t13+t12+t11+t10+t9+t8+t7+t6+t5+t4+t3+t2)/35
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




python starpy.py 2.472543716 0.5 21.79722 0.19 0.025761 8733 316.92 3.87

Best fit [t, tau] values found by starpy for input parameters are : [ 7.6012507518 0.322877696627 ]

Best fit [t, tau] values found by starpy for input parameters are : [ 7.60632485936 0.325843373704 ]

Best fit [t, tau] values found by starpy for input parameters are : [ 7.61801495164 0.306381749944 ]

####################################################
# UGC11680NED01
# Centroid 36.840103428583,31.7533234547072
# Redshift = 0.0258870075940476
# DA = 0.496342681361508
# Re_kpc = 8.33855704687333
# Mass St uncorrected = 11.3379781183796
# Ha/Hb = 4.05258832920205 +- 175.892528409773
# mean Av = 1.09253465255246 +- 43.4025156570507 (mag)
# weighted Av = 1.67441667512879  (mag)
# mean mu_SFR = 0.000459608188054147 +- 0.00460587585555642  (Msun/yr^2/spaxel)
# SFR = 4.14393641113 (Msun/yr^2)
# log_SFR = 0.61741308252211 log10(Msun/yr^2)
##########################################
# IONIZATION 
# -1 -> No class 
#  0 -> No Gas 
#  1 -> SF  <Kewley and abs_EW(Ha)>3
#  2 -> wAGN  >Kewley and 3>abs_EW(Ha)>6
#  3 -> sAGN  >Kewley and abs_EW(Ha)>6
#  4 -> pAGB  abs_EW(Ha)<3
#  5 -> Shocks?
# CENTRAL IONIZATION = 3
# AREA_SFR_pAGB = 851 0.110490825677263
# OH_O3N2 (M13) = 8.51418296177332 +- 0.0699017597810228
# Mass_corr = 11.037303973818 vs. 11.3379781183796
# VEL_SYS_GAS = 7771.84063657731 km/s
# VEL_SYS_SSP = 7749.61863719108 km/s
# Av_ssp = 0.2952967146091
./proc_elines.pl.20150531.pl UGC11680NED01 175 0.46 20 1/xs
proc_elines.pl UGC11680NED01 106.78 0.05 33.6 1/xs
 perl ./proc_elines.pl.20150531.pl UGC11680NED01 175 0.46 20 1/xs


# NGC4047
# Centroid 36.3928552589237,32.4056359348288
# Redshift = 0.011356738508533
# DA = 0.22216311659558
# Re_kpc = 1.43295210204149
# Mass St uncorrected = 11.2199110114187
# Ha/Hb = 3.55998236693923 +- 58.8365667271996
# mean Av = 0.686282673783659 +- 16.5272073461941 (mag)
# weighted Av = 1.07200725328495  (mag)
# mean mu_SFR = 0.000403174012116852 +- 0.00123591131584436  (Msun/yr^2/spaxel)
# SFR = 3.64225628449595 (Msun/yr^2)
# log_SFR = 0.561370501323378 log10(Msun/yr^2)
##########################################
# IONIZATION 
# -1 -> No class 
#  0 -> No Gas 
#  1 -> SF  <Kewley and abs_EW(Ha)>3
#  2 -> wAGN  >Kewley and 3>abs_EW(Ha)>6
#  3 -> sAGN  >Kewley and abs_EW(Ha)>6
#  4 -> pAGB  abs_EW(Ha)<3
#  5 -> Shocks?
# CENTRAL IONIZATION = 1
# AREA_SFR_pAGB = 3194 0.404132535566878
# OH_O3N2 (M13) = 8.5241556895091 +- 0.0665782901044763
# Mass_corr = 10.5465874583847 vs. 11.2199110114187
# VEL_SYS_GAS = 3387.55190537336 km/s
# VEL_SYS_SSP = 3421.77719929937 km/s
# Av_ssp = 0.34784838797332


/usr/local/texlive/2015/bin/i386-linux