

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


sum_chi_1= ((img_mask1-img)**2)
norm_1=0.5*(img_mask1+img)
sum_tot_1= sum_chi_1/norm_1
chi_j_1= np.sum(sum_tot_1, axis=1)
chi_ugc11680= np.sum(chi_j_1)
norm_df= 36*38
chi_red_ugc11680= chi_ugc11680/norm_df

sum_chi_2= ((img_mask2-img)**2)
norm_2=0.5*(img_mask2+img)
sum_tot_2= sum_chi_2/norm_2
chi_j_2= np.sum(sum_tot_2, axis=1)
chi_ic0540= np.sum(chi_j_2)
norm_df= 36*38
chi_red_ic0540= chi_ic0540/norm_df

sum_chi_3= ((img_mask3-img)**2)
norm_3=0.5*(img_mask3+img)
sum_tot_3= sum_chi_3/norm_3
chi_j_3= np.sum(sum_tot_3, axis=1)
chi_ngc0426= np.sum(chi_j_3)
norm_df= 36*38
chi_red_ngc0426= chi_ngc0426/norm_df

sum_chi_4= ((img_mask4-img)**2)
norm_4=0.5*(img_mask4+img)
sum_tot_4= sum_chi_4/norm_4
chi_j_4= np.sum(sum_tot_4, axis=1)
chi_ngc0833= np.sum(chi_j_4)
norm_df= 36*38
chi_red_ngc0833= chi_ngc0833/norm_df

sum_chi_5= ((img_mask5-img)**2)
norm_5=0.5*(img_mask5+img)
sum_tot_5= sum_chi_5/norm_5
chi_j_5= np.sum(sum_tot_5, axis=1)
chi_ngc1667= np.sum(chi_j_5)
norm_df= 36*38
chi_red_ngc1667= chi_ngc1667/norm_df

sum_chi_6= ((img_mask6-img)**2)
norm_6=0.5*(img_mask6+img)
sum_tot_6= sum_chi_6/norm_6
chi_j_6= np.sum(sum_tot_6, axis=1)
chi_ngc2410= np.sum(chi_j_6)
norm_df= 36*38
chi_red_ngc2410= chi_ngc2410/norm_df

sum_chi_7= ((img_mask7-img)**2)
norm_7=0.5*(img_mask7+img)
sum_tot_7= sum_chi_7/norm_7
chi_j_7= np.sum(sum_tot_7, axis=1)
chi_ngc2623= np.sum(chi_j_7)
norm_df= 36*38
chi_red_ngc2623= chi_ngc2623/norm_df

sum_chi_8= ((img_mask8-img)**2)
norm_8=0.5*(img_mask8+img)
sum_tot_8= sum_chi_8/norm_8
chi_j_8= np.sum(sum_tot_8, axis=1)
chi_ngc2639= np.sum(chi_j_8)
norm_df= 36*38
chi_red_ngc2639= chi_ngc2639/norm_df

sum_chi_9= ((img_mask9-img)**2)
norm_9=0.5*(img_mask9+img)
sum_tot_9= sum_chi_9/norm_9
chi_j_9= np.sum(sum_tot_9, axis=1)
chi_ngc3106= np.sum(chi_j_9)
norm_df= 36*38
chi_red_ngc3106= chi_ngc3106/norm_df

sum_chi_10= ((img_mask10-img)**2)
norm_10=0.5*(img_mask10+img)
sum_tot_10= sum_chi_10/norm_10
chi_j_10= np.sum(sum_tot_10, axis=1)
chi_ngc3160= np.sum(chi_j_10)
norm_df= 36*38
chi_red_ngc3160= chi_ngc3160/norm_df

sum_chi_11= ((img_mask11-img)**2)
norm_11=0.5*(img_mask11+img)
sum_tot_11= sum_chi_11/norm_11
chi_j_11= np.sum(sum_tot_11, axis=1)
chi_ngc3303= np.sum(chi_j_11)
norm_df= 36*38
chi_red_ngc3303= chi_ngc3303/norm_df

sum_chi_12= ((img_mask12-img)**2)
norm_12=0.5*(img_mask12+img)
sum_tot_12= sum_chi_12/norm_12
chi_j_12= np.sum(sum_tot_12, axis=1)
chi_ngc5635= np.sum(chi_j_12)
norm_df= 36*38
chi_red_ngc5635= chi_ngc5635/norm_df

sum_chi_13= ((img_mask13-img)**2)
norm_13=0.5*(img_mask13+img)
sum_tot_13= sum_chi_13/norm_13
chi_j_13= np.sum(sum_tot_13, axis=1)
chi_ngc5675= np.sum(chi_j_13)
norm_df= 36*38
chi_red_ngc5675= chi_ngc5675/norm_df

sum_chi_14= ((img_mask14-img)**2)
norm_14=0.5*(img_mask14+img)
sum_tot_14= sum_chi_14/norm_14
chi_j_14= np.sum(sum_tot_14, axis=1)
chi_ngc5739= np.sum(chi_j_14)
norm_df= 36*38
chi_red_ngc5739= chi_ngc5739/norm_df

sum_chi_15= ((img_mask15-img)**2)
norm_15=0.5*(img_mask15+img)
sum_tot_15= sum_chi_15/norm_15
chi_j_15= np.sum(sum_tot_15, axis=1)
chi_ngc6323= np.sum(chi_j_15)
norm_df= 36*38
chi_red_ngc6323= chi_ngc6323/norm_df

sum_chi_16= ((img_mask16-img)**2)
norm_16=0.5*(img_mask16+img)
sum_tot_16= sum_chi_16/norm_16
chi_j_16= np.sum(sum_tot_16, axis=1)
chi_ngc6394= np.sum(chi_j_16)
norm_df= 36*38
chi_red_ngc6394= chi_ngc6394/norm_df

sum_chi_17= ((img_mask17-img)**2)
norm_17=0.5*(img_mask17+img)
sum_tot_17= sum_chi_17/norm_17
chi_j_17= np.sum(sum_tot_17, axis=1)
chi_ngc7466= np.sum(chi_j_17)
norm_df= 36*38
chi_red_ngc7466= chi_ngc7466/norm_df

sum_chi_18= ((img_mask18-img)**2)
norm_18=0.5*(img_mask18+img)
sum_tot_18= sum_chi_18/norm_18
chi_j_18= np.sum(sum_tot_18, axis=1)
chi_ngc7738= np.sum(chi_j_18)
norm_df= 36*38
chi_red_ngc7738= chi_ngc7738/norm_df

sum_chi_19= ((img_mask19-img)**2)
norm_19=0.5*(img_mask19+img)
sum_tot_19= sum_chi_19/norm_19
chi_j_19= np.sum(sum_tot_19, axis=1)
chi_ugc00005= np.sum(chi_j_19)
norm_df= 36*38
chi_red_ugc00005= chi_ugc00005/norm_df

sum_chi_20= ((img_mask20-img)**2)
norm_20=0.5*(img_mask20+img)
sum_tot_20= sum_chi_20/norm_20
chi_j_20= np.sum(sum_tot_20, axis=1)
chi_ugc00987= np.sum(chi_j_20)
norm_df= 36*38
chi_red_ugc00987= chi_ugc00987/norm_df

sum_chi_21= ((img_mask21-img)**2)
norm_21=0.5*(img_mask21+img)
sum_tot_21= sum_chi_21/norm_21
chi_j_21= np.sum(sum_tot_21, axis=1)
chi_ugc03995= np.sum(chi_j_21)
norm_df= 36*38
chi_red_ugc03995= chi_ugc03995/norm_df

sum_chi_22= ((img_mask22-img)**2)
norm_22=0.5*(img_mask22+img)
sum_tot_22= sum_chi_22/norm_22
chi_j_22= np.sum(sum_tot_22, axis=1)
chi_ugc05771= np.sum(chi_j_22)
norm_df= 36*38
chi_red_ugc05771= chi_ugc05771/norm_df

sum_chi_23= ((img_mask23-img)**2)
norm_23=0.5*(img_mask23+img)
sum_tot_23= sum_chi_23/norm_23
chi_j_23= np.sum(sum_tot_23, axis=1)
chi_ugc09711= np.sum(chi_j_23)
norm_df= 36*38
chi_red_ugc09711= chi_ugc09711/norm_df

sum_chi_24= ((img_mask24-img)**2)
norm_24=0.5*(img_mask24+img)
sum_tot_24= sum_chi_24/norm_24
chi_j_24= np.sum(sum_tot_24, axis=1)
chi_mcg_02_02_030= np.sum(chi_j_24)
norm_df= 36*38
chi_red_MCG_02_02_030= chi_mcg_02_02_030/norm_df

chi_red_ugc11680
chi_red_ugc09711
chi_red_ugc05771
chi_red_ugc03995
chi_red_ugc00987
chi_red_ugc00005
chi_red_ngc7738
chi_red_ngc7466
chi_red_ngc6394
chi_red_ngc6323
chi_red_ngc5739
chi_red_ngc5675
chi_red_ngc5635
chi_red_ngc3303
chi_red_ngc3160
chi_red_ngc3106
chi_red_ngc2639
chi_red_ngc2623
chi_red_ngc2410
chi_red_ngc1667
chi_red_ngc0833
chi_red_ngc0426
chi_red_ic0540
chi_red_MCG_02_02_030

chi_tot=np.array([3.4729173553279491,5.4075144791484773,1.7265234372027134,3.3700246209748657,6.8812918076930689,4.1015754330938874,2.4704381755833476,4.8222624525863349,6.0511088308201177,3.0362021243289381,3.1439912696072607,1.3496896471099378,1.9538755457513528,1.8533262672127828,2.139781384290095,4.1716996165374916,23.42367742274438,10.82744163883533,0.67051713995132334,7.2823982899103017,6.1452447425749588,3.4332305213587033,3.4029097527128211,3.078775306774113])


plt.hist(np.random.noncentral_chisquare(3, 20, 100000), bins=200, normed=True)



