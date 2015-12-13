#!/usr/bin/perl
#
use Statistics::OLS;
use Math::FFT;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);
use Math::Approx;
use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Slatec;
use PDL::Image2D;
use PDL::Graphics::LUT;
use Carp;
use PDL::Graphics::PGPLOT;

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";


require "/home/jeffrey/Documents/AGNS/my.pl"; # "/home/sanchez/sda2/code/R3D/my.pl";

$speed_of_light=299792.458;
if ($#ARGV<4) {
    print "USE: proc_elines.pl NAME PA[degress] Ellipticity Re DEV\n";
    exit;
}

#################################################
# Stellar Template!
#$ssp_template="/disk-a/sanchez/ppak/legacy/gsd61_156.fits";

#
# Mstar -> Mass Remaining after mass loss.
#
$base_info="/home/jeffrey/Documents/AGNS/BASE.gsd01";
open(FH,"<$base_info");
$line=<FH>; chop($line); $nbase=$line;
for ($i=0;$i<$nbase;$i++) {
    $line=<FH>;    chop($line); my @data=split(" ",$line);
    $mstar[$i]=$data[4];    
}
close(FH);



$ssp_template="/home/jeffrey/Documents/AGNS/gsd61_156.fits";
$pdl_flux_c_ini=rfits($ssp_template);
($nxf,$nf)=$pdl_flux_c_ini->dims();
for ($iii=0;$iii<$nf;$iii++) {
    $header="NAME".$iii;
    $name[$iii]=$pdl_flux_c_ini->hdr->{$header};
    $name_min=$name[$iii];
    $name_min =~ s/spec_ssp_//;
    $name_min =~ s/.spec//;    
    $name_min =~ s/.dat//;
    ($AGE,$MET)=split("_",$name_min);
    if ($AGE =~ "Myr") {
	$age=$AGE;
	$age =~ s/Myr//;
	$age=$age/1000;
    } else {
	$age=$AGE;
	$age =~ s/Gyr//;
    }
    $met=$MET;
    $met =~ s/z/0\./;
    $age_mod[$iii]=$age;
    $met_mod[$iii]=$met;
    $header="NORM".$iii;

    $val_ml=$pdl_flux_c_ini->hdr->{$header};
    if ($val_ml!=0) {
	$ml[$iii]=1/$val_ml;
    } else {
	$ml[$iii]=1;
    }
    $ML{$age}=$ML{$age}+$ml[$iii]/4;
    $age_met=$age."-".$met;
    $INDEX{$age_met}=$iii;
   #print "$age,$met,$ML{$age}\n";
}

#exit;

###############################################
# CONFIG!
#
$spax_scale=1;

$n_V_seg=3; $label[$n_V_seg][1]="V-band intensity Segmented";

$n_V=0; $label[$n_V][1]="V-band intensity";
$n_vel_Ha=96; $label[$n_vel_Ha][0]="H\ga velocity map";
$n_F_Ha=45; $label[$n_F_Ha][0]="H\ga intensity map";
$n_eF_Ha=249; $label[$n_eF_Ha][0]="H\ga error intensity map";
$n_EW_Ha=198; $label[$n_EW_Ha][0]="EW(H\ga)";

$n_vel_Hb=79; $label[$n_vel_Hb][0]="H\gb velocity map";
$n_F_Hb=28; $label[$n_F_Hb][0]="H\gb intensity map";
$n_eF_Hb=232; $label[$n_eF_Hb][0]="H\gb error intensity map";
$n_EW_Hb=181; $label[$n_EW_Hb][0]="EW(H\gb)";

$n_vel_OIII=77; $label[$n_vel_OIII][0]="[OIII] velocity map";
$n_F_OIII=26; $label[$n_F_OIII][0]="[OIII] intensity map";
$n_eF_OIII=230; $label[$n_eF_OIII][0]="[OIII] error intensity map";
$n_EW_OIII=179; $label[$n_EW_OIII][0]="EW([OIII])";

$n_vel_OIIIw=78; $label[$n_vel_OIIIw][0]="[OIII]4959 velocity map";
$n_F_OIIIw=27; $label[$n_F_OIIIw][0]="[OIII]4959 intensity map";
$n_eF_OIIIw=231; $label[$n_eF_OIIIw][0]="[OIII]4959 error intensity map";
$n_EW_OIIIw=180; $label[$n_EW_OIIIw][0]="EW([OIII]4959)";

$n_vel_OII=51; $label[$n_vel_OII][0]="[OII] velocity map";
$n_F_OII=0; $label[$n_F_OII][0]="[OII] intensity map";
$n_eF_OII=204; $label[$n_eF_OII][0]="[OII] error intensity map";
$n_EW_OII=153; $label[$n_EW_OII][0]="EW([OII])";

$n_vel_NII=97; $label[$n_vel_NII][0]="[NII] velocity map";
$n_F_NII=46; $label[$n_F_NII][0]="[NII] intensity map";
$n_eF_NII=250; $label[$n_eF_NII][0]="[NII] error intensity map";
$n_EW_NII=199; $label[$n_EW_NII][0]="EW([NII])";

$n_vel_NIIw=98; $label[$n_vel_NII][0]="[NII] velocity map";
$n_F_NIIw=47; $label[$n_F_NII][0]="[NII] intensity map";
$n_eF_NIIw=251; $label[$n_eF_NII][0]="[NII] error intensity map";
$n_EW_NIIw=200; $label[$n_EW_NII][0]="EW([NII])";

$n_vel_SII6717=100; $label[$n_vel_SII6717][0]="[SII6717] velocity map";
$n_F_SII6717=49; $label[$n_F_SII6717][0]="[SII6717] intensity map";
$n_eF_SII6717=253; $label[$n_eF_SII6717][0]="[SII6717] error intensity map";
$n_EW_SII6717=202; $label[$n_EW_SII6717][0]="EW([SII6717])";

$n_vel_SII6731=101; $label[$n_vel_SII6731][0]="[SII6731] velocity map";
$n_F_SII6731=50; $label[$n_F_SII6731][0]="[SII6731] intensity map";
$n_eF_SII6731=254; $label[$n_eF_SII6731][0]="[SII6731] error intensity map";
$n_EW_SII6731=203; $label[$n_EW_SII6731][0]="EW([SII6731])";





$n_vel_ssp=13; $label[$n_vel_ssp][1]="Stellar velocity map";
$n_disp_ssp=15; $label[$n_disp_ssp][1]="Stellar vel. dispersion map";
$n_scale_ssp=2; $label[$n_scale_ssp][1]="Deszonification map";
$n_Mass_ssp=18;  $label[$n_Mass_ssp][1]="Stellar Mass (Msun/arcsec\u2\d)";
$n_Av_ssp=11;
$n_eAv_ssp=12;

###############################################
#
###############################################

$name=$ARGV[0];
$elines_file="flux_elines.".$ARGV[0].".cube.fits.gz";
$ssp_file=$ARGV[0].".SSP.cube.fits.gz";
$sfh_file=$ARGV[0].".SFH.cube.fits.gz";
$pa=$ARGV[1];
$elip=$ARGV[2];
$Re=$ARGV[3]; #spaxels
$DEV=$ARGV[4];
if (($DEV !~ "XS")&&($DEV !~ "xs")) {
    $DEV=$name."_proc_elines".$DEV;
}

if ($#ARGV==6) {
    $XC=$ARGV[5];
    $YC=$ARGV[6];
}




#print "$DEV\n"; exit;

$pa_rad=((90-$pa)*(3.1416/180));
if ($elip<0.99) {
    $median_incl_PDL=acos(sqrt((1-$elip**2-0.13**2)/(1-0.13**2)));
    $median_incl=$median_incl_PDL->at(0);
    $median_incl_deg=180*($median_incl/3.1416);
} else {
    $median_incl_deg=0;
}
#print "$pa_rad,$median_incl_deg\n";

#
# Reading Elines
#
print "####################################################\n";
print "# $name\n";
#print "# Reading elines\n";
$pdl_elines=rfits($elines_file);
#print "# Reading ssp\n";
$pdl_ssp=rfits($ssp_file);
#print "# Reading sfh\n";
$pdl_sfh=rfits($sfh_file);
#print "# Reading sfh\n";

my $h=$pdl_sfh->gethdr();

$pdl_sfh_lum_age=$pdl_sfh->slice(":,:,156:194");
$pdl_sfh_lum_age->sethdr($h);
$pdl_sfh_y_200=$pdl_sfh->slice(":,:,156:165");

my $pdl_junk=sumover($pdl_sfh_y_200->xchg(0,2));
$pdl_map_y_200=$pdl_junk->xchg(0,1);


$V_img_seg=$pdl_ssp->slice(":,:,($n_V_seg)");




#
# Centroid
#
$V_img=$pdl_ssp->slice(":,:,($n_V)");
($nx,$ny)=$V_img->dims;    
$FoV=sqrt((0.5*$nx)**2+(0.5*$ny)**2)/$Re;
$flux_V=sum($V_img);
# Zero-points, assuming
# 10^-16 erg/s/cm^2
$mag_g=-2.5*log10($flux_V/(4.705*10**(7)));#*10**(-0.4*($g));

#$mag_V=$f[1]=4.705*10**(7)*10**(-0.4*($g));

$img_mask=zeroes($nx,$ny);
$V_mask="/home/jeffrey/Documents/AGNS/".$NAME.".V500.mask.fits";
$nx1=$nx; $ny1=$ny;
eval {
    $img_mask=rfits($V_mask);
    ($nx1,$ny1)=$img_mask->dims;    
};

$img_mask_new=ones($nx,$ny);;


if ($nx1>$nx) {
    $nx1=$nx;
} 

if ($ny1>$ny) {
    $ny1=$ny;
} 
 
for ($ii=0;$ii<$nx1;$ii++) {                                                     
    for ($jj=0;$jj<$ny1;$jj++) {                                                 
	$val_mask=$img_mask->at($ii,$jj);                                                 
	if ($val_mask==2) {
	    set($img_mask_new,$ii,$jj,0);
	} else {
	    set($img_mask_new,$ii,$jj,1);
	}
	$val=$V_img->at($ii,$jj);                                                 
	if ($val<0.01) {
	    set($img_mask_new,$ii,$jj,0);
	}	

    }
}
$val_max=-1e12;  

if (($XC==0)&&($YC==0)) {                                                  
$XC=0;                                                                         
$YC=0;                        
$SUM=0;  
$nx0=0;#int($nx*0.33);
$nx1=$nx;#int($nx*2*0.33);
$ny0=0;#int($nx*0.33);
$ny1=$ny;#int($nx*2*0.33);
for ($ii=$nx0;$ii<$nx1;$ii++) {
    for ($jj=$ny0;$jj<$ny1;$jj++) {
        $val=$V_img->at($ii,$jj);
        $val_mask=$img_mask->at($ii,$jj);
#       print "$ii/$nx0 $jj/$ny0 $val $val_mask $nx $ny\n";                                                                                             
        if ($val_mask!=2) {
            $XC=$XC+$ii*($val)**5;
            $YC=$YC+$jj*($val)**5;
            $SUM=$SUM+($val)**5;
        }
    }
}
$xc=$XC/$SUM;
$yc=$YC/$SUM;
} else {
    $xc=$XC;
    $yc=$YC;

}

my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*0.5;
$$h{CDELT1}=0.5;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*0.5;
$$h{CDELT2}=0.5;
$$h{CRPIX2}=1;
$V_img->sethdr($h);



print "# Centroid $xc,$yc\n";
#
# END CENTROID!!!!!
#


$pdl_r=zeroes($nx,$ny);
$pdl_r_arcsec=zeroes($nx,$ny);
$pdl_x_rot=zeroes($nx,$ny);
$pdl_y_rot=zeroes($nx,$ny);

for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$x_S=$xc-$i;
	$y_S=$j-$yc;       
	$x_ROT=$x_S*cos($pa_rad)+$y_S*sin($pa_rad);
	$y_ROT=$x_S*sin($pa_rad)-$y_S*cos($pa_rad);    
	if (abs($median_incl_deg)>35) {
	    $x_ROT=$x_ROT/cos($median_incl);
	}
	$r=sqrt($x_ROT**2+$y_ROT**2)/$Re;
	$r_arcsec=sqrt($x_ROT**2+$y_ROT**2);
	set($pdl_x_rot,$i,$j,$x_ROT);
	set($pdl_y_rot,$i,$j,$y_ROT);
	set($pdl_r,$i,$j,$r);
	set($pdl_r_arcsec,$i,$j,$r_arcsec);
    }
}

$outfile=$name.".p_e.pdl_r.fits";
$pdl_r->wfits($outfile);
#$pdl_x_rot->wfits("pdl_x_rot.fits");
#$pdl_y_rot->wfits("pdl_y_rot.fits");

################################################
# Redshift
################################################






$pdl_vel_Ha=$pdl_elines->slice(":,:,($n_vel_Ha)");
$pdl_F_Ha=$pdl_elines->slice(":,:,($n_F_Ha)");
$pdl_eF_Ha=$pdl_elines->slice(":,:,($n_eF_Ha)");
$pdl_EW_Ha=$pdl_elines->slice(":,:,($n_EW_Ha)");

$pdl_vel_Hb=$pdl_elines->slice(":,:,($n_vel_Hb)");
$pdl_F_Hb=$pdl_elines->slice(":,:,($n_F_Hb)");
$pdl_eF_Hb=$pdl_elines->slice(":,:,($n_eF_Hb)");
$pdl_EW_Hb=$pdl_elines->slice(":,:,($n_EW_Hb)");

$pdl_vel_OIII=$pdl_elines->slice(":,:,($n_vel_OIII)");
$pdl_F_OIII=$pdl_elines->slice(":,:,($n_F_OIII)");
$pdl_eF_OIII=$pdl_elines->slice(":,:,($n_eF_OIII)");
$pdl_EW_OIII=$pdl_elines->slice(":,:,($n_EW_OIII)");


$pdl_vel_OIIIw=$pdl_elines->slice(":,:,($n_vel_OIIIw)");
$pdl_F_OIIIw=$pdl_elines->slice(":,:,($n_F_OIIIw)");
$pdl_eF_OIIIw=$pdl_elines->slice(":,:,($n_eF_OIIIw)");
$pdl_EW_OIIIw=$pdl_elines->slice(":,:,($n_EW_OIIIw)");


$pdl_vel_OII=$pdl_elines->slice(":,:,($n_vel_OII)");
$pdl_F_OII=$pdl_elines->slice(":,:,($n_F_OII)");
$pdl_eF_OII=$pdl_elines->slice(":,:,($n_eF_OII)");
$pdl_EW_OII=$pdl_elines->slice(":,:,($n_EW_OII)");

$pdl_vel_NII=$pdl_elines->slice(":,:,($n_vel_NII)");
$pdl_F_NII=$pdl_elines->slice(":,:,($n_F_NII)");
$pdl_eF_NII=$pdl_elines->slice(":,:,($n_eF_NII)");
$pdl_EW_NII=$pdl_elines->slice(":,:,($n_EW_NII)");

$pdl_vel_NIIw=$pdl_elines->slice(":,:,($n_vel_NIIw)");
$pdl_F_NIIw=$pdl_elines->slice(":,:,($n_F_NIIw)");
$pdl_eF_NIIw=$pdl_elines->slice(":,:,($n_eF_NIIw)");
$pdl_EW_NIIw=$pdl_elines->slice(":,:,($n_EW_NIIw)");

$pdl_vel_SII6717=$pdl_elines->slice(":,:,($n_vel_SII6717)");
$pdl_F_SII6717=$pdl_elines->slice(":,:,($n_F_SII6717)");
$pdl_eF_SII6717=$pdl_elines->slice(":,:,($n_eF_SII6717)");
$pdl_EW_SII6717=$pdl_elines->slice(":,:,($n_EW_SII6717)");

$pdl_vel_SII6731=$pdl_elines->slice(":,:,($n_vel_SII6731)");
$pdl_F_SII6731=$pdl_elines->slice(":,:,($n_F_SII6731)");
$pdl_eF_SII6731=$pdl_elines->slice(":,:,($n_eF_SII6731)");
$pdl_EW_SII6731=$pdl_elines->slice(":,:,($n_EW_SII6731)");

$pdl_F_SII=$pdl_F_SII6717+$pdl_F_SII6731;
$pdl_eF_SII=sqrt($pdl_eF_SII6717**2+$pdl_eF_SII6731**2);
$pdl_EW_SII=$pdl_EW_SII6717+$pdl_EW_SII6731;

$pdl_V=$pdl_ssp->slice(":,:,($n_V)");
$pdl_vel_ssp=$pdl_ssp->slice(":,:,($n_vel_ssp)");
$pdl_disp_ssp=$pdl_ssp->slice(":,:,($n_disp_ssp)");

$pdl_V->inplace->setvaltobad(0);
@V_sum=stats($pdl_V);
$pdl_V_rat=$pdl_V/$V_sum[0];

$i0=int($xc-0.5);
$i1=int($xc+1);
$j0=int($yc-0.5);
$j1=int($yc+1);
@a_vel_Ha_cen=stats($pdl_vel_Ha->slice("$i0:$i1,$j0:$j1"));
@a_vel_ssp_cen=stats($pdl_vel_ssp->slice("$i0:$i1,$j0:$j1"));
#$vel_center_Ha=$pdl_vel_Ha->at($xy,$yc);



#print "$a_vel_Ha_cen[0] $a_vel_ssp_cen[0]\n";

if (abs($a_vel_Ha_cen[0]-$a_vel_ssp_cen[0])<150) {
    $redshift=0.5*($a_vel_Ha_cen[0]+$a_vel_ssp_cen[0])/$speed_of_light;
} else {
    $redshift=$a_vel_ssp_cen[0]/$speed_of_light;
}

if ($redshift<0.0001) {
    exit;
}

@cosmo=cosm($redshift,71,0.27,0.73);

$modz=$cosmo[4]; # Distance Modulus
$ratio=3.08567758e24;
$DL=10**(($modz-25)/5)*$ratio;
$L=4*3.1416*($DL**2)*1/(1+$redshift);
$DA=$DL/(1+$redshift)**2/206.264806/$ratio;  # Angular Scale (Kpc/");

$pdl_r_kpc=$pdl_r_arcsec*$DA;
$Re_kpc=$Re*$DA/2; # If Re is in spaxels for MaNGA!!!!


print "# Redshift = $redshift\n";
print "# DA = $DA\n";
print "# Re_kpc = $Re_kpc\n";

########################################
# Mass
#$n_scale_ssp=2; $label[$n_scale_ssp][1]="Deszonification map";
$pdl_scale=$pdl_ssp->slice(":,:,($n_scale_ssp)");
$pdl_log_Mass=$pdl_ssp->slice(":,:,($n_Mass_ssp)");
$pdl_Av_ssp=$pdl_ssp->slice(":,:,($n_Av_ssp)");
$pdl_e_Av_ssp=$pdl_ssp->slice(":,:,($n_eAv_ssp)");
$pdl_log_Mass=$pdl_log_Mass+0.4*$pdl_Av_ssp;
$pdl_Mass=$pdl_scale*10**($pdl_log_Mass);
$pdl_Mass=$pdl_Mass*$img_mask_new;
$pdl_Mass->inplace->setnantobad;
$Mass=sum($pdl_Mass);
$log_Mass=log10($Mass);

$pdl_map_y_200=$pdl_map_y_200*$pdl_scale;
#$pdl_map_y_200->wfits("pdl.fits"); exit;


$pdl_log_Mass=$pdl_log_Mass*$img_mask_new;
#$pdl_Sigma_Mass=$pdl_log_Mass+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2

#
# Old correction:
#


print "# Mass St uncorrected = $log_Mass\n";

#
# SFR
#
#$pdl_Ha_Hb=$pdl_F_Ha/$pdl_F_Hb;
#$pdl_e_Ha_Hb=$pdl_eF_Ha/$pdl_F_Hb;#+$pdl_F_Ha*($pdl_eF_Hb)/($pdl_F_Hb**2);

$pdl_e_Ha_rat=$pdl_eF_Ha/$pdl_F_Ha;
($pdl_e_Ha_clean,$pdl_mask_e_Ha,$stats_e_Ha)=create_mask($pdl_e_Ha_rat,-1,0.3);

$pdl_e_Hb_rat=$pdl_eF_Hb/$pdl_F_Hb;
($pdl_e_Hb_clean,$pdl_mask_e_Hb,$stats_e_Hb)=create_mask($pdl_e_Hb_rat,-1,3);

($pdl_EW_Ha_SFR,$pdl_mask_EW_Ha_SFR,$stats_EW_Ha)=create_mask($pdl_EW_Ha,-1e10,-10);

$pdl_Ha_Hb=($pdl_F_Ha/$pdl_F_Hb);
$pdl_Ha_Hb_clean=$pdl_Ha_Hb*($pdl_mask_e_Ha);#$pdl_mask_EW_Ha_SFR;
$pdl_Ha_Hb_clean->inplace->setnantobad;
$pdl_Ha_Hb_clean->inplace->setvaltobad(0);
@stats_Ha_Hb=stats($pdl_Ha_Hb_clean);
print "# Ha/Hb = $stats_Ha_Hb[2] +- $stats_Ha_Hb[1]\n";
$pdl_Ha_Hb->wfits("pdl_Ha_Hb.fits");


$rat_Ha_Hb=$stats_Ha_Hb[2];#/$stats_Ha[2];
$e_rat_Ha_Hb=$stats_Ha_Hb[1];#/$stats_Ha[2];
$Rv=3.1;
$a_1=A_l($Rv,4861);
$a_2=A_l($Rv,6562);
if ($rat_Ha_Hb>2.86) {
    $Av=(2.5*log10($rat_Ha_Hb/2.86))/($a_1-$a_2);
    $e_Av=$e_rat_Ha_Hb/$rat_Ha_Hb;
} else {
    $Av=0;
    $e_Av=0;
}

print "# mean Av = $Av +- $e_Av (mag)\n";	    

$pdl_Ha_Hb_clean->inplace->setbadtoval(0);
$pdl_Ha_Hb_clean->inplace->clip(2.86,1000);
$pdl_log_Ha_Hb_clean=$pdl_Ha_Hb_clean/2.86;
$pdl_log_Ha_Hb_clean->inplace->log10; 
$pdl_Av=(2.5*$pdl_log_Ha_Hb_clean)/($a_1-$a_2);

my @a_Av_w;
my $n_Av_w=0;
my $sum_Ha_flux=0;
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val_Ha=$pdl_F_Ha->at($i,$j);
	$val_EW_Ha=$pdl_EW_Ha->at($i,$j);
	$val_e_Ha=$pdl_eF_Ha->at($i,$j);
	$val_Av=$pdl_Av->at($i,$j);
	if (($val_Ha>(3*abs($val_e_Ha)))&&($val_EW_Ha<-10)) {
	    $Av_w=$Av_w+$val_Ha*$val_Av;
	    $sum_Ha_flux=$sum_Ha_flux+$val_Ha;
	    $n_Av_w++;
	}
    }
}

#@stats_Av_w=stats(pdl(@a_Av_w));
if ($sum_Ha_flux>0) {
    $Av_w=$Av_w/$sum_Ha_flux;
} else {
    $Av_w=0;
}
print "# weighted Av = $Av_w  (mag)\n";	    

$pdl_Av->inplace->clip(0,30);
$pdl_Av->inplace->setvaltobad(0);
$pdl_Av_org=$pdl_Av->copy;
$pdl_Av->inplace->setbadtoval($Av_w);


$outfile=$name.".p_e.Av_gas.fits";
$pdl_Av->wfits($outfile);
#$pdl_Av->wfits("pdl_Av.fits");

$Al=A_l($Rv,6562);
#
$pdl_F_Ha_cor=$pdl_F_Ha*10**(0.4*$Av_w*$Al);
#$pdl_F_Ha_cor=$pdl_F_Ha*10**(0.4*$pdl_Av*$Al);
$pdl_L_Ha_cor=$pdl_F_Ha_cor*$L*(1e-16);
$pdl_SFR=$pdl_L_Ha_cor*0.8*(1e-41);
$pdl_e_Ha_rat=$pdl_eF_Ha/$pdl_F_Ha;
($pdl_e_Ha_clean,$pdl_mask_e_Ha,$stats_e_Ha)=create_mask($pdl_e_Ha_rat,-1,1);
$pdl_SFR=$pdl_SFR*$pdl_mask_e_Ha;
$pdl_SFR->inplace->setvaltobad(0);
$pdl_SFR->inplace->setnantobad;
@stats_SFR=stats($pdl_SFR);
$SFR=sum($pdl_SFR);
print "# mean mu_SFR = $stats_SFR[2] +- $stats_SFR[1]  (Msun/yr^2/spaxel)\n";
print "# SFR = $SFR (Msun/yr^2)\n";
$lSFR=log10($SFR);
print "# log_SFR = $lSFR log10(Msun/yr^2)\n";

#($pdl_EW_Ha_SFR,$pdl_mask_EW_Ha_SFR,$stats_EW_Ha)=create_mask($pdl_EW_Ha,-1e10,-1);

$pdl_log_SFR=$pdl_SFR->log10;
$pdl_Sigma_SFR=$pdl_log_SFR-log10(($DA)**2)-6; # Msun/yr/pc^2

$outfile=$name.".p_e.Sigma_SFR.fits";
$pdl_Sigma_SFR->wfits($outfile);

#
# BPT and ionization classification
#

($pdl_e_Ha_clean,$pdl_mask_e_Ha,$stats_e_Ha)=create_mask($pdl_e_Ha_rat,-1,2);

$pdl_e_OIII_rat=$pdl_eF_OIII/$pdl_F_OIII;
($pdl_e_OIII_clean,$pdl_mask_e_OIII,$stats_e_OIII)=create_mask($pdl_e_OIII_rat,-1,3);
$pdl_OIII_Hb=($pdl_F_OIII/$pdl_F_Hb);
$pdl_OIII_Hb_clean=$pdl_OIII_Hb*($pdl_mask_e_Ha)*($pdl_mask_e_Hb)*($pdl_mask_e_OIII);
$pdl_OIII_Hb_clean->inplace->log10;
$pdl_OIII_Hb_clean->inplace->setnantobad;

$pdl_e_OII_rat=$pdl_eF_OII/$pdl_F_OII;
($pdl_e_OII_clean,$pdl_mask_e_OII,$stats_e_OII)=create_mask($pdl_e_OII_rat,-1,3);
$pdl_OII_Hb=($pdl_F_OII/$pdl_F_Hb);
$pdl_OII_Hb_clean=$pdl_OII_Hb*($pdl_mask_e_Ha)*($pdl_mask_e_Hb)*($pdl_mask_e_OII);
$pdl_OII_Hb_clean->inplace->log10;
$pdl_OII_Hb_clean->inplace->setnantobad;

$pdl_e_NII_rat=$pdl_eF_NII/$pdl_F_NII;
($pdl_e_NII_clean,$pdl_mask_e_NII,$stats_e_NII)=create_mask($pdl_e_NII_rat,-1,3);
$pdl_NII_Ha=($pdl_F_NII/$pdl_F_Ha);
$pdl_NII_Ha_clean=$pdl_NII_Ha*($pdl_mask_e_Ha)*($pdl_mask_e_Ha)*($pdl_mask_e_NII);
$pdl_NII_Ha_clean->inplace->log10;
$pdl_NII_Ha_clean->inplace->setnantobad;

$pdl_e_SII_rat=$pdl_eF_SII/$pdl_F_SII;
($pdl_e_SII_clean,$pdl_mask_e_SII,$stats_e_SII)=create_mask($pdl_e_SII_rat,-1,3);
$pdl_SII_Ha=($pdl_F_SII/$pdl_F_Ha);
$pdl_SII_Ha_clean=$pdl_SII_Ha*($pdl_mask_e_Ha)*($pdl_mask_e_Ha)*($pdl_mask_e_SII);
$pdl_SII_Ha_clean->inplace->log10;
$pdl_SII_Ha_clean->inplace->setnantobad;




$pdl_OIII_OII=($pdl_F_OIII/$pdl_F_OII);
$pdl_OIII_OII_clean=$pdl_OIII_OII*($pdl_mask_e_OIII)*($pdl_mask_e_OII);#*($pdl_mask_e_SII);
$pdl_OIII_OII_clean->inplace->log10;
$pdl_OIII_OII_clean->inplace->setnantobad;


#$pdl_OIII_Hb_clean->wfits("pdl_OIII_Hb.fits");
#$pdl_OII_Hb_clean->wfits("pdl_OII_Hb.fits");
#$pdl_NII_Ha_clean->wfits("pdl_NII_Ha.fits");

$outfile=$name.".p_e.OIII_Hb.fits";
$pdl_OIII_Hb_clean->wfits($outfile);
$outfile=$name.".p_e.OII_Hb.fits";
$pdl_OII_Hb_clean->wfits($outfile);
$outfile=$name.".p_e.NII_Ha.fits";
$pdl_NII_Ha_clean->wfits($outfile);

#
# We classify the ionization regions
#
#  0 -> No GAS
# -1 -> No class
#  1 -> SF  <Kewley and abs_EW(Ha)>3
#  2 -> wAGN  >Kewley and 3>abs_EW(Ha)>6
#  3 -> sAGN  >Kewley and abs_EW(Ha)>6
#  4 -> pAGB  abs_EW(Ha)<3
#  5 -> Shocks?
$pdl_ion_class=zeroes($nx,$ny);
$pdl_OH_O3N2=zeroes($nx,$ny);
$pdl_OH_N2=zeroes($nx,$ny);
$pdl_OH_ONS=zeroes($nx,$ny);
$pdl_OH_R23=zeroes($nx,$ny);
$pdl_log_U_iter=zeroes($nx,$ny);

#$pdl_F_Ha_cor=$pdl_F_Ha*10**(0.4*$pdl_Av*$Al);
$pdl_F_Hb_cor=$pdl_F_Hb*10**(0.4*$pdl_Av*$Al);
$pdl_F_OII_cor=$pdl_F_OII*10**(0.4*$pdl_Av*$Al);
$pdl_F_OIII_cor=$pdl_F_OIII*10**(0.4*$pdl_Av*$Al);
$pdl_F_SII_cor=$pdl_F_SII*10**(0.4*$pdl_Av*$Al);
$pdl_F_NII_cor=$pdl_F_NII*10**(0.4*$pdl_Av*$Al);

$pdl_R2=$pdl_F_OII_cor/$pdl_F_Hb_cor;
$pdl_N2=1.333*$pdl_F_NII_cor/$pdl_F_Hb_cor;
$pdl_S2=$pdl_F_SII_cor/$pdl_F_Hb_cor;
$pdl_R3=1.333*$pdl_F_OIII_cor/$pdl_F_Hb_cor;
$pdl_P=$pdl_R3/($pdl_R3+$pdl_R2);
$pdl_O23=1.333*$pdl_F_OIII_cor/$pdl_F_OII_cor;
$pdl_log_U_P=1.22*($pdl_O23->log10)-2.25;


$pdl_EW_Ha_cor=$pdl_EW_Ha*$pdl_V_rat;
$sum_SFR=0;
$sum_SFR_pure=0;
$sum_GAS=0;
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
#
# If there is gas emission
#
	if (($pdl_mask_e_Ha->at($i,$j))==1) {
	    $sum_GAS=$sum_GAS+1;
	    $X=$pdl_NII_Ha_clean->at($i,$j);
	    $Y=$pdl_OIII_Hb_clean->at($i,$j);
	    $EW=abs($pdl_EW_Ha_cor->at($i,$j));
	    $disp=$pdl_disp_ssp->at($i,$j);
	    $r_kpc=$pdl_r_kpc->at($i,$j);
	    if (($X ne "BAD")&&($Y ne "BAD")) {
		$OH_O3N2=8.505-0.221*($Y-$X)+0.04;
		$OH_N2=8.743+0.462*$X-0.04;
		my $N2=$pdl_N2->at($i,$j);
		my $R2=$pdl_R2->at($i,$j);
		my $R3=$pdl_R3->at($i,$j);
		my $S2=$pdl_S2->at($i,$j);
		my $O23=$pdl_O23->at($i,$j);
		my $P=$pdl_P->at($i,$j);
		if (($R3>0)&&($R2>0)&&($N2>0)&&($S2>0)) {
		    if (log10($N2)>=-1) {
			$OH_ONS=8.277+0.657*$P-0.399*log10($R3)-0.061*log10($N2/$R2)+0.005*log10($S2/$R2);
		    } else {
			if (log10($S2/$R2)>=-0.25) {
			    $OH_ONS=8.816-0.733*$P+0.454*log10($R3)+0.710*log10($N2/$R2)-0.337*log10($S2/$R2);
			} else {
			    $OH_ONS=8.774-1.855*$P+1.517*log10($R3)+0.304*log10($N2/$R2)+0.328*log10($S2/$R2);
			}
		    }
		} else {
		    $OH_ONS=nan;
		}
		
		$OH_iter=$OH_O3N2;
		$n_iter=0;
		$Delta_iter=1e12;
		$converg=0.001;
		$break=0;
		my $log_q;
		do {
		    my $y=log10($O23);
		    my $z=$OH_iter;
		    my $x=log10($R2+$R3);
#		    print "$n_iter $OH_iter\n";
		    $log_q=(32.81-1.153*($y**2)+$z*(-3.396-0.025*$y+0.1444*($y**2)))/(4.603-0.3119*$y-0.163*$y**2+$z*(-0.48+0.0271*$y+0.02037*($y**2)))+0.3;
		    if ($z>8.3) {
			$OH_iter=9.72-0.777*$x-0.951*($x**2)-0.072*($x**3)-0.811*($x**4)-$log_q*(0.0737-0.0713*$x-0.141*($x**2)+0.0373*($x**3)-0.058*($x**4))-0.6;
		    } else {
			$OH_iter=9.40+4.65*$x-3.17*($x**2)-$log_q*(0.272-0.547*$x-0.513*($x**2))-0.6;
		    }
		    $Delta_iter=abs($OH_iter-$z);
		    if ($n_iter<4) {
			$Delta_iter=100;
		    }-log10(3e10);
		    $n_iter++;
		} while (($Delta_iter>$converg)&&($n_iter<20));
		set($pdl_log_U_iter,$i,$j,($log_q-log10(3e10)));
		


		$Y_KEWLEY=0.61/($X-0.47)+1.19;
		# pAGBs
		if ($EW<3) {
		    $sum_SFR=$sum_SFR+1;
		    set($pdl_OH_O3N2,$i,$j,$OH_O3N2);
		    set($pdl_OH_N2,$i,$j,$OH_N2);
		    set($pdl_OH_ONS,$i,$j,$OH_ONS);
		    set($pdl_OH_R23,$i,$j,$OH_iter);
		    set($pdl_ion_class,$i,$j,4);
		} else {
		    # SF
		    if ($Y<=$Y_KEWLEY) {
			set($pdl_ion_class,$i,$j,1);
#			if ($EW>20) {
			$sum_SFR=$sum_SFR+1;
			$sum_SFR_pure=$sum_SFR_pure+1;
			set($pdl_OH_O3N2,$i,$j,$OH_O3N2);
			set($pdl_OH_N2,$i,$j,$OH_N2);
			set($pdl_OH_ONS,$i,$j,$OH_ONS);
			set($pdl_OH_R23,$i,$j,$OH_iter);
#			}
		    } else {
			# AGNs/LINER
			if ($Y>1.5*$X) {
			    if ($EW>6) {
				# sAGN	
				set($pdl_ion_class,$i,$j,3);
			    } else {
				# wAGN	
				set($pdl_ion_class,$i,$j,2);
			    }						    
			} else {
			    if (($r_kpc>3)&&($disp>10)) {
				# Shocks?
				set($pdl_ion_class,$i,$j,5);
			    } else {
				# pAGBs
				set($pdl_OH_O3N2,$i,$j,$OH_O3N2);
				set($pdl_OH_N2,$i,$j,$OH_N2);
				set($pdl_OH_ONS,$i,$j,$OH_ONS);
				set($pdl_ion_class,$i,$j,4);
				set($pdl_OH_R23,$i,$j,$OH_iter);
			    }
			}
		    }
		}
	    } else {
		set($pdl_ion_class,$i,$j,-1);
	    }
	}
    }
}





my $h = {NAXIS=>2, NAXIS1=>$nx, NAXIS2=>$ny, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="diagnostic BTP+EW_Ha";
$$h{VAL0}="No GAS 0";
$$h{VAL1}="SF     1";
$$h{VAL2}="wAGN   2";
$$h{VAL3}="sAGN   3";
$$h{VAL4}="pAGB   4";
$$h{VAL5}="Shocks 5";
$pdl_ion_class->sethdr($h);
$outfile=$name.".p_e.ion_class.fits";
$pdl_ion_class->wfits($outfile);

$i0=int($xc-1);
$i1=int($xc+1);
$j0=int($yc-1);
$j1=int($yc+1);
@a_ion_cen=stats($pdl_ion_class->slice("$i0:$i1,$j0:$j1"));
#print "@a_ion_cen\n";

#$ion_cen_type=$pdl_ion_class->($xc,$yc);





#$pdl_ion_class->wfits("pdl_ion_class.fits");
#$pdl_EW_Ha->wfits("pdl_EW_Ha.fits");
#$pdl_EW_Ha_cor->wfits("pdl_EW_Ha_cor.fits");
$pdl_OH_O3N2->inplace->setvaltobad(0);
#$pdl_OH_O3N2->wfits("pdl_OH.fits");
$outfile=$name.".p_e.OH_O3N2.fits";
$pdl_OH_O3N2->wfits($outfile);



@stats_OH=stats($pdl_OH_O3N2);
$frac_area_SFR=$sum_SFR/(6*($nx/2)*sqrt(($nx/2)**2-($nx/4)**2));
$frac_area_SFR_pure=$sum_SFR_pure/(6*($nx/2)*sqrt(($nx/2)**2-($nx/4)**2));
$frac_area_GAS=$sum_GAS/(6*($nx/2)*sqrt(($nx/2)**2-($nx/4)**2));


print "##########################################\n";
print "# IONIZATION \n";
print "# -1 -> No class \n";
print "#  0 -> No Gas \n";
print "#  1 -> SF  <Kewley and abs_EW(Ha)>3\n";
print "#  2 -> wAGN  >Kewley and 3>abs_EW(Ha)>6\n";
print "#  3 -> sAGN  >Kewley and abs_EW(Ha)>6\n";
print "#  4 -> pAGB  abs_EW(Ha)<3\n";
print "#  5 -> Shocks?\n";
print "# CENTRAL IONIZATION = $a_ion_cen[2]\n";
print "# AREA_SFR_pAGB = $sum_SFR $frac_area_SFR\n";

$OH_O3N2=$stats_OH[0];
$e_OH_O3N2=$stats_OH[1];
if ($frac_area_SFR<0.05) {
    $OH_O3N2="nan";
    $e_OH_O3N2="nan";
}
print "# OH_O3N2 (M13) = $OH_O3N2 +- $e_OH_O3N2\n";

$outfile=$name.".p_e.OH_O3N2_norm.fits";
$pdl_OH_O3N2_norm=$pdl_OH_O3N2-$OH_O3N2+8.69;
$pdl_OH_O3N2_norm->wfits($outfile);




#
# Add colors
#
my $nc=254;
my @r,@l,@g,@b;
while($cmap=<DATA>) {
    chop($cmap);
    @data=split(" ",$cmap);
    $nc=$data[0];
    $r[$nc-1]=$data[1]/255;
    $g[$nc-1]=$data[2]/255;
    $b[$nc-1]=$data[3]/255;
    $l[$nc]=$nc/255;
}
my $bright=1.0; 
my $contrast=0.5;
$r[0]=1.0;
$g[0]=1.0;
$b[0]=1.0;
$pdl_ctable=pdl(pdl(@r),pdl(@g),pdl(@b),pdl(@l));

my $nc=254;
my @r,@l,@g,@b;
open(VDATA,"</disk-b/sanchez/ppak/legacy/VDATA.cmap.txt");
while($cmap=<VDATA>) {
    chop($cmap);
    @data=split(" ",$cmap);
    $nc=$data[0];
    $r[$nc-1]=$data[1]/255;
    $g[$nc-1]=$data[2]/255;
    $b[$nc-1]=$data[3]/255;
    $l[$nc]=$nc/255;
}
close(VDATA);
my $bright=1.0; 
my $contrast=0.5;
$r[0]=1.0;
$g[0]=1.0;
$b[0]=1.0;
my $pdl_ctable_vel=pdl(pdl(@r),pdl(@g),pdl(@b),pdl(@l));






#plot_diag($pdl_NII_Ha_clean,$pdl_OIII_Hb_clean,5,"55/xs");
#plot_diag_color($pdl_NII_Ha_clean,$pdl_OIII_Hb_clean,$pdl_r,$pdl_ctable,"56/xs","$name","R/Re");
#plot_diag_color($pdl_NII_Ha_clean,$pdl_OIII_Hb_clean,$pdl_EW_Ha,$pdl_ctable,"55/xs","$name","EW(H\\ga)");

#
# Uncomment!
#

$dev="1".$DEV;
plot_diag_color_scale($pdl_NII_Ha_clean,$pdl_OIII_Hb_clean,$pdl_r,$pdl_ctable,$dev,"$name","R/Re",0,3);


#$pdl_abs_EW_Ha=abs($pdl_EW_Ha);


#
$pdl_abs_EW_Ha=abs($pdl_EW_Ha_cor);
$pdl_log_EW_Ha=$pdl_abs_EW_Ha->log10;

#
$dev="2".$DEV;
($pdl_MAP_diag,$pdl_fMAP_diag)=plot_diag_color_scale($pdl_NII_Ha_clean,$pdl_OIII_Hb_clean,$pdl_log_EW_Ha,$pdl_ctable,$dev,"$name","log10 |EW(H\\ga)|",0,2.5);

$outfile=$name.".p_e.MAP_diag.fits";
$pdl_MAP_diag->wfits($outfile);
$outfile=$name.".p_e.fMAP_diag_EW.fits";
$pdl_fMAP_diag->wfits($outfile);





#
$dev="3".$DEV;
($pdl_MAP_r_OH,$pdl_fMAP_r_OH)=plot_XY_color_scale($pdl_r,$pdl_OH_O3N2,$pdl_log_EW_Ha,$pdl_ctable,$dev,0,2.5,8.11,8.79,0,2.5,"R/Re","12+log(O/H)","log10 |EW(H\\ga)|","$name");

$outfile=$name.".p_e.MAP_r_OH.fits";
$pdl_MAP_r_OH->wfits($outfile);
$outfile=$name.".p_e.fMAP_r_OH_EW.fits";
$pdl_fMAP_r_OH->wfits($outfile);




#

$dev="7".$DEV;
$pdl_map_y_200=$pdl_map_y_200*100;
($pdl_MAP_diag,$pdl_fMAP_diag)=plot_diag_color_scale($pdl_NII_Ha_clean,$pdl_OIII_Hb_clean,$pdl_map_y_200,$pdl_ctable,$dev,"$name","%",0,49);
$outfile=$name.".p_e.MAP_diag.fits";
$pdl_MAP_diag->wfits($outfile);
$outfile=$name.".p_e.fMAP_diag_young.fits";
$pdl_fMAP_diag->wfits($outfile);





$pdl_abs_EW_Ha_no_cor=abs($pdl_EW_Ha);
$pdl_log_EW_Ha_no_cor=$pdl_abs_EW_Ha_no_cor->log10;

#
$dev="9".$DEV;
($pdl_MAP_diag,$pdl_fMAP_diag)=plot_diag_color_scale($pdl_NII_Ha_clean,$pdl_OIII_Hb_clean,$pdl_log_EW_Ha_no_cor,$pdl_ctable,$dev,"$name","log10 |EW(H\\ga)|",0,2.5);

$outfile=$name.".p_e.MAP_diag.fits";
$pdl_MAP_diag->wfits($outfile);
$outfile=$name.".p_e.fMAP_diag_EW_no_cor.fits";
$pdl_fMAP_diag->wfits($outfile);



my $h=$pdl_sfh_lum_age->gethdr();
$nx0=156;
$nx1=194;
for ($i=$nx0;$i<$nx1+1;$i++) {
    $key="FILE_".$i;
    $val=$h->{$key};
    ($junk,$age_now,$junk)=split(/\_/,$val);
    $age[$i-$nx0]=$age_now;
}





my $h=$pdl_sfh_lum_age->gethdr();
$pdl_sfh_lum_age_V=$pdl_sfh_lum_age*$V_img;
$pdl_sfh_lum_age_Mass=$pdl_sfh_lum_age_V->copy;
$pdl_sfh_lum_age_Mass_Time=$pdl_sfh_lum_age_V->copy;
$nx0=156;
$nx1=194;

my ($nx_now,$ny_now,$nz_now)=$pdl_sfh_lum_age->dims;
$pdl_sfh_met=zeroes($nx_now,$ny_now,$nz_now);
$pdl_sfh_met_Mass=zeroes($nx_now,$ny_now,$nz_now);
$pdl_Lum_coadd=zeroes($nx_now,$ny_now,$nz_now);
$pdl_Lum_coadd_Mass=zeroes($nx_now,$ny_now,$nz_now);

for ($i=$nx0;$i<$nx1+1;$i++) {
    $key="FILE_".$i;
    $val=$h->{$key};
    $age_now=$val;
    $cut="map.CS.".$name."_";
    $age_now =~ s/$cut//g;
    $cut="_NORM_age.fits.gz";
    $age_now =~ s/$cut//g; 
    $age[$i-$nx0]=$age_now;
    my $I=$i-$nx0;
    my $name="NAME".$I;
    $$h{$name}=apr(9+log10($age_now));
    my $t=$pdl_sfh_lum_age_Mass->slice(":,:,($I)");
    my $t_Time=$pdl_sfh_lum_age_Mass_Time->slice(":,:,($I)");
    my $tt=$pdl_sfh_lum_age_V->slice(":,:,($I)");
    $age_now=1.0*$age_now;
    my $ML_now=0;
    my $dist_min=1e12;
    for ($iii=0;$iii<$nf;$iii++) {
	my $dist=abs($age_now-1.*$age_mod[$iii]);
	if ($dist<$dist_min) {
	    $ML_mod=$ML{$age_mod[$iii]};
	    $dist_min=$dist;
	    $ML_now=$ML_mod;
	}
    }

    $ML_final{$age_now}=$ML_now;
    #print "$I $mstar[$I]\n";
    $t .= $tt*$ML_now*$L*(1e-16)/3.826e33*$mstar[$I]; # Observed Mass
    $t_Time .= $tt*$ML_now*$L*(1e-16)/3.826e33; # Mass at the time, corrected by mass-loses

#
# Run for metals!
#
    my $pdl_Lum_coadd_t=$pdl_Lum_coadd->slice(":,:,($I)");
    my $pdl_Lum_coadd_Mass_t=$pdl_Lum_coadd_Mass->slice(":,:,($I)");
    my $ml_coadd=0;
    for ($j=0;$j<4;$j++) {
	$met_now=$met_mod[$j];
	my $age_met=$age_now."-".$met_now;
	$II=$INDEX{$age_met};
	my $lmet_now=log10($met_now);
	my $t=$pdl_sfh_met->slice(":,:,($I)");
	my $Lum=$pdl_sfh->slice(":,:,($II)");
	$t .=$t+$lmet_now*$Lum;#*$ml[$II];
	$pdl_Lum_coadd_t .= $pdl_Lum_coadd_t+$Lum;
	my $t=$pdl_sfh_met_Mass->slice(":,:,($I)");
	$t .=$t+$lmet_now*$Lum*$ml[$II];
	my @stats=stats($Lum);
	
	$pdl_Lum_coadd_Mass_t .= $pdl_Lum_coadd_Mass_t+$Lum*$ml[$II];
	$ml_coadd=$ml_coadd+$ml[$II];

    }
#    $pdl_Lum_coadd->inplace->setvaltobad(0);
#    $pdl_Lum_coadd_Mass->inplace->setvaltobad(0);

#
# ************
#
#    my $t=$pdl_sfh_met->slice(":,:,($I)");
#
#    $t .= 10**($t/$pdl_Lum_coadd_t);
#    $t .= $t/$pdl_Lum_coadd_t;
#    my $t=$pdl_sfh_met_Mass->slice(":,:,($I)");
#
#    $t .= 10**($t/$pdl_Lum_coadd_Mass_t);
#    $t .= $t/$pdl_Lum_coadd_Mass_t;
    #$t .= $t/$ml_coadd;

#
# 
#



}


$pdl_Lum_coadd_C=$pdl_Lum_coadd->copy;
$pdl_Lum_coadd_Mass_C=$pdl_Lum_coadd_Mass->copy;

$pdl_Lum_coadd->inplace->setvaltobad(0);
$pdl_Lum_coadd_Mass->inplace->setvaltobad(0);

$pdl_sfh_met_C=$pdl_sfh_met->copy;
$pdl_sfh_met_Mass_C=$pdl_sfh_met_Mass->copy;

$pdl_sfh_met=10**($pdl_sfh_met/$pdl_Lum_coadd);
$pdl_sfh_met_Mass=10**($pdl_sfh_met_Mass/$pdl_Lum_coadd_Mass);

#$pdl_sfh_met->inplace->setvaltobad(0);
#$pdl_sfh_met_Mass->inplace->setvaltobad(0);
#$pdl_sfh_met=10**($pdl_sfh_met);
#$pdl_sfh_met_Mass=10**($pdl_sfh_met_Mass);
$pdl_sfh_met->inplace->setnantobad;
$pdl_sfh_met_Mass->inplace->setnantobad;
$pdl_sfh_met->inplace->setbadtoval(0);
$pdl_sfh_met_Mass->inplace->setbadtoval(0);
$pdl_sfh_met=$pdl_sfh_met*$img_mask_new;
$pdl_sfh_met_Mass=$pdl_sfh_met_Mass*$img_mask_new;
#$pdl_sfh_met->wfits("pdl.fits");
#exit;


$pdl_sfh_lum_age_V_C=$pdl_sfh_lum_age*$V_img;
$pdl_sfh_lum_age_Mass_C=$pdl_sfh_lum_age_Mass->copy;
#$pdl_sfh_met_C=$pdl_sfh_met->copy;
$pdl_SFR_Time=$pdl_sfh_lum_age_Mass_Time->copy;
#=$pdl_sfh_lum_age_V->copy;



for ($I=($nx1-$nx0);$I>-1;$I--) {
    my $name="NAME".$I;
    $$h{$name}=$age[$I];

#    if ($I!=($nx1-$nx0)) {
#	my $t=$pdl_SFR_Time->slice(":,:,($I)");
#	my $t2=$pdl_sfh_lum_age_Mass_Time->slice(":,:,($I)");
#	my $D_time=$age[$I+1]-$age[$I];
#	$t .= $t2/$D_time;	
 #   } else {

    if ($I!=0) {	
	my $t=$pdl_SFR_Time->slice(":,:,($I)");
	my $t1=$pdl_sfh_lum_age_Mass_Time->slice(":,:,($I)");
	my $I2=$I-1;
	my $t2=$pdl_sfh_lum_age_Mass_Time->slice(":,:,($I2)");
	my $D_time=$age[$I]-$age[$I-1];
	$D_time=$D_time*1e9;
	$t .= 0.5*($t2+$t1)/$D_time;
    } else {
	my $D_time=$age[$I+1]-$date[$I];
	$D_time=$D_time*1e9;
	my $t=$pdl_SFR_Time->slice(":,:,($I)");
	my $t1=$pdl_sfh_lum_age_Mass_Time->slice(":,:,($I)");
        $t .= $t1/$D_time;
    }

  #  }


    if ($I==($nx1-$nx0)) {
	my $t=$pdl_sfh_lum_age_Mass_C->slice(":,:,($I)");
	my $tt=$pdl_sfh_lum_age_Mass->slice(":,:,($I)");
	$t .= $tt;
	my $t=$pdl_sfh_lum_age_V_C->slice(":,:,($I)");
	my $tt=$pdl_sfh_lum_age_V->slice(":,:,($I)");
	$t .= $tt;
	my $t=$pdl_sfh_met_C->slice(":,:,($I)");
	my $tt=$pdl_sfh_met_C->slice(":,:,($I)");
	$t .= $tt;
	my $t=$pdl_sfh_met_Mass_C->slice(":,:,($I)");
	my $tt=$pdl_sfh_met_Mass_C->slice(":,:,($I)");
	$t .= $tt;
	my $t=$pdl_Lum_coadd_C->slice(":,:,($I)");
	my $tt=$pdl_Lum_coadd_C->slice(":,:,($I)");
	$t .= $tt;
	my $t=$pdl_Lum_coadd_Mass_C->slice(":,:,($I)");
	my $tt=$pdl_Lum_coadd_Mass_C->slice(":,:,($I)");
	$t .= $tt;
    } else {
	my $I0=$I+1;
	my $t=$pdl_sfh_lum_age_Mass_C->slice(":,:,($I)");
	my $t1=$pdl_sfh_lum_age_Mass_C->slice(":,:,($I0)");
	$t .= $t+$t1;
	my $t=$pdl_sfh_lum_age_V_C->slice(":,:,($I)");
	my $t1=$pdl_sfh_lum_age_V_C->slice(":,:,($I0)");
	$t .= $t+$t1;
	my $t=$pdl_sfh_met_C->slice(":,:,($I)");
	my $t1=$pdl_sfh_met_C->slice(":,:,($I0)");
	$t .= $t+$t1;
	my $t=$pdl_sfh_met_Mass_C->slice(":,:,($I)");
	my $t1=$pdl_sfh_met_Mass_C->slice(":,:,($I0)");
	$t .= $t+$t1;
	my $t=$pdl_Lum_coadd_C->slice(":,:,($I)");
	my $t1=$pdl_Lum_coadd_C->slice(":,:,($I0)");
	$t .= $t+$t1;
	my $t=$pdl_Lum_coadd_Mass_C->slice(":,:,($I)");
	my $t1=$pdl_Lum_coadd_Mass_C->slice(":,:,($I0)");
	$t .= $t+$t1;
    }

}

$pdl_Lum_coadd_C->inplace->setvaltobad(0);
$pdl_Lum_coadd_Mass_C->inplace->setvaltobad(0);

$pdl_sfh_met_C=10**($pdl_sfh_met_C/$pdl_Lum_coadd_C);
$pdl_sfh_met_Mass_C=10**($pdl_sfh_met_Mass_C/$pdl_Lum_coadd_Mass_C);

$pdl_sfh_met_C=$pdl_sfh_met_C*$img_mask_new;
$pdl_sfh_met_Mass_C=$pdl_sfh_met_Mass_C*$img_mask_new;

#
#
#

#$pdl_SFR_Time->wfits("pdl.fits");
#exit;

#$pdl_sfh_lum_age_Mass->wfits("pdl.fits");
#$pdl_sfh_lum_age_Mass_C->wfits("pdl_C.fits");

#exit;

#$pdl_sfh_lum_age_Mass->inplace->log10();
$pdl_sfh_lum_age_Mass->sethdr($h);
#$pdl_sfh_lum_age_Mass->wfits("pdl.fits");
$pdl_sfh_lum_age->sethdr($h);
$pdl_sfh_lum_age_V->sethdr($h);


my $pdl_Mass_corr=$pdl_sfh_lum_age_Mass_C->slice(":,:,(0)")+0.4*$pdl_Av_ssp;
$pdl_Mass_corr=$pdl_Mass_corr*$img_mask_new;
$pdl_Mass_corr->inplace->setnantobad;
$pdl_log_Mass_corr=$pdl_Mass_corr->log10;
$Mass_corr=sum($pdl_Mass_corr);
$log_Mass_corr=log10($Mass_corr);
print "# Mass_corr = $log_Mass_corr vs. $log_Mass\n";
#$pdl_Sigma_Mass=$pdl_log_Mass_corr+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2
$pdl_Sigma_Mass=$pdl_log_Mass_corr+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2


$outfile=$name.".p_e.Sigma_Mass.fits";
$pdl_Sigma_Mass->wfits($outfile);

$pdl_Sigma_Mass_uncorr=$pdl_log_Mass+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2

#$dev="b4".$DEV;
#plot_XY_color_scale($pdl_Sigma_Mass_uncorr,$pdl_OH_O3N2,$pdl_log_EW_Ha,$pdl_ctable,$dev,-0.5,4.2,8.11,8.79,0,2.49,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","12+log(O/H)","log10 |EW(H\\ga)|","$name");

$dev="4".$DEV;
($pdl_MAP_SMass_OH,$pdl_fMAP_SMass_OH)=plot_XY_color_scale($pdl_Sigma_Mass,$pdl_OH_O3N2,$pdl_log_EW_Ha,$pdl_ctable,$dev,,-1.1,3.3,8.11,8.79,0,2.49,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","12+log(O/H)","log10 |EW(H\\ga)|","$name");

$outfile=$name.".p_e.MAP_SMass_OH.fits";
$pdl_MAP_SMass_OH->wfits($outfile);
$outfile=$name.".p_e.fMAP_SMass_OH_EW.fits";
$pdl_fMAP_SMass_OH->wfits($outfile);

#$pdl_MAP_SMass_OH->wfits("pdl0.fits");
#$pdl_fMAP_SMass_OH->wfits("pdl1.fits");

#$pdl_MAP_diag->wfits("bpt0.fits");
#$pdl_fMAP_diag->wfits("bpt1.fits");

$dev="5".$DEV;
plot_MAP_color_scale($pdl_MAP_SMass_OH,$pdl_fMAP_SMass_OH,$pdl_ctable,$dev,-1.1,3.3,8.11,8.79,0,2.49,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","12+log(O/H)","log10 |EW(H\\ga)|","$name",1,1);

#plot_MAP_color_scale($pdl_MAP_SMass_OH,$pdl_fMAP_SMass_OH,$pdl_ctable,"junk.ps/CPS",0.5,4.2,8.11,8.79,0,2.49,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","12+log(O/H)","log10 |EW(H\\ga)|","$name",1,1);

$dev="6".$DEV;
plot_diag_MAP_color_scale($pdl_MAP_diag,$pdl_fMAP_diag,$pdl_ctable,$dev,-1.59,0.74,-1.2,1.2,0,2.49,"log([NII]\\gl6583/H\\ga)","log([OIII]\\gl5007/H\\gb)","log10 |EW(H\\ga)|","$name",1,1);

#$dev="7".$DEV;
#plot_diag_MAP_color_scale($pdl_MAP_diag,$pdl_fMAP_diag,$pdl_ctable,"junk.ps/CPS",-1.2,0.6,-1.2,1.2,0,2.49,"log([NII]\\gl6583/H\\ga)","log([OIII]\\gl5007/H\\gb)","log10 |EW(H\\ga)|","$name",1,1);



$dev="8".$DEV;
#
($pdl_MAP_y200_EW_Ha,$pdl_fMAP_y200_EW_Ha)=plot_XY_color_scale($pdl_map_y_200,$pdl_log_EW_Ha,$pdl_Sigma_Mass,$pdl_ctable,$dev,0,99,-1.9,3.49,0,2.49,"% Lum. Age<200 Myr","log10 |EW(H\\ga)|","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","$name");
#($pdl_MAP_SMass_OH,$pdl_fMAP_SMass_OH)=plot_XY_color_scale($pdl_map_y_200,$pdl_log_EW_Ha,$pdl_map_y_200,$pdl_ctable,$dev,0,99,-1.9,3.49,0,99,"% Lum. Age<200 Myr","log10 |EW(H\\ga)|","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","$name");

$outfile=$name.".p_e.MAP_y200_EW_Ha.fits";
$pdl_MAP_y200_EW_Ha->wfits($outfile);
$outfile=$name.".p_e.fMAP_y200_EW_Ha_SMass.fits";
$pdl_fMAP_y200_EW_Ha->wfits($outfile);






$pdl_rad_sfh_lum_age=radial_sum_ring($pdl_sfh_lum_age,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_sfh_lum_age->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH";
$pdl_rad_sfh_lum_age->sethdr($h);


$pdl_rad_sfh_lum_age_V=radial_sum_ring($pdl_sfh_lum_age_V,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_sfh_lum_age->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH";
$pdl_rad_sfh_lum_age_V->sethdr($h);

$pdl_rad_sfh_lum_age_Mass=radial_sum_ring($pdl_sfh_lum_age_Mass,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_sfh_lum_age->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH";
$pdl_rad_sfh_lum_age_Mass->sethdr($h);


$pdl_rad_sfh_met=radial_sum_ring($pdl_sfh_met,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_sfh_met->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH met";
$pdl_rad_sfh_met->sethdr($h);

$pdl_rad_sfh_met_C=radial_sum_ring($pdl_sfh_met_C,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_sfh_met_C->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH met";
$pdl_rad_sfh_met_C->sethdr($h);




$dev="10".$DEV;
my $r_max=0.1*$n_r;
$pdl_label_Y=0.1*pdl([0..$n_r]);
$pdl_label_X=pdl(@age);
#print "pdl_age = $pdl_label_X\n";
$pdl_label_X->inplace->log10;
$pdl_label_X=rint(10*(9+$pdl_label_X))/10;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_sfh_lum_age, $kernel, {Boundary => Reflect});
$pdl_smooth->inplace->clip(0,1);
$pdl_smooth->inplace->setnantobad;
$pdl_smooth->inplace->setbadtoval(0);
plot_MAP_color_scale_table($pdl_smooth,$pdl_rad_sfh_lum_age,$pdl_ctable,$dev,-0.05,0.35,"log(age/yr)","R/Re","Lum. Weight","$name",0,1,$pdl_label_X,$pdl_label_Y);

$outfile=$name.".p_e.rad_SFH_lum.fits";
$pdl_rad_sfh_lum_age->wfits($outfile);

$dev="11".$DEV;
($pdl_MAP_SMass_SFR,$pdl_fMAP_SMass_SFR)=plot_XY_color_scale($pdl_Sigma_Mass,$pdl_Sigma_SFR,$pdl_map_y_200,$pdl_ctable,$dev,,-1.1,3.3,-10.5,-5.5,0,49,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","log(\\gS\\dSFR\\u/M\\d\\m9\\u yr\\u-1\\d pc\\u-2\\d)","% Lum. Age<200 Myr","$name");

$outfile=$name.".p_e.MAP_SMass_SFR.fits";
$pdl_MAP_SMass_SFR->wfits($outfile);
$outfile=$name.".p_e.fMAP_SMass_SFR_y_200.fits";
$pdl_fMAP_SMass_SFR->wfits($outfile);

#################################################################
# Here it was the pdl_Mass_corr
# 



# Mass normalized to time=0;
$pdl_Mass_Norm=$pdl_sfh_lum_age_Mass_C/$pdl_sfh_lum_age_Mass_C->slice(":,:,(0)");
$pdl_rad_Mass_Norm=radial_sum_ring($pdl_Mass_Norm,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_Mass_Norm->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH Mass Norm";
$pdl_rad_Mass_Norm->sethdr($h);

# Normalized metals at time=0
#$pdl_sfh_met_C->setnantobad();

$pdl_met_Norm=$pdl_sfh_met_C/$pdl_sfh_met_C->slice(":,:,(0)");
$pdl_rad_met_Norm=radial_sum_ring($pdl_met_Norm,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_met_Norm->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH Met Norm";
$pdl_rad_met_Norm->sethdr($h);


$pdl_met_Norm_Mass=$pdl_sfh_met_Mass_C/$pdl_sfh_met_Mass_C->slice(":,:,(0)");
$pdl_rad_met_Norm_Mass=radial_sum_ring($pdl_met_Norm_Mass,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_met_Norm_Mass->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH Met Norm Mass";
$pdl_rad_met_Norm->sethdr($h);

#exit;

#
#for ($I=1;$I<38;$I++) {
# Look for a loop
#
my $I=12; # 39
#$pdl_Sigma_Mass=$pdl_log_Mass+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2
my $pdl_Mass_Time_I=$pdl_sfh_lum_age_Mass_C->slice(":,:,($I)")*$img_mask_new;;
my $pdl_Sigma_Mass_Time=$pdl_Mass_Time_I->log10-log10(($DA)**2)-6;#+$pdl_scale->log10; # Msun/pc^2
$pdl_Sigma_Mass_Time=$pdl_Sigma_Mass_Time+0.4*$pdl_Av_ssp;
my $pdl_SFR_Time_I=$pdl_SFR_Time->slice(":,:,($I)")*$img_mask_new;;
$pdl_log_SFR_I=$pdl_SFR_Time_I->log10;
$pdl_Sigma_SFR_I=$pdl_log_SFR_I-log10(($DA)**2)-6;#+$pdl_scale->log10;; # Msun/yr/pc^2
$pdl_Sigma_SFR_I=$pdl_Sigma_SFR_I+0.4*$pdl_Av_ssp;
$pdl_Sigma_SFR_I=$pdl_Sigma_SFR_I-$V_img->log10+$V_img_seg->log10;#+$pdl_scale->log10;
$pdl_Sigma_Mass_Time=$pdl_Sigma_Mass_Time-$V_img->log10+$V_img_seg->log10;#+$pdl_scale->log10;

#$dev="22/xs";
$dev="11_Time".$DEV;
($pdl_MAP_SMass_SFR_Time,$pdl_fMAP_SMass_SFR_Time)=plot_XY_color_scale($pdl_Sigma_Mass_Time,$pdl_Sigma_SFR_I,$pdl_r,$pdl_ctable,$dev,,-1.1,3.3,-10.5,-5.5,0,3,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d) \@age=$age[$I]","log(\\gS\\dSFR\\u/M\\d\\m9\\u yr\\u-1\\d pc\\u-2\\d)","R/Re","$name");

#<stdin>;
#$pdl_Sigma_SFR_I->wfits("pdl_1.fits");
#$pdl_Sigma_Mass_Time->wfits("pdl_2.fits");
#$pdl_scale->wfits("scale.fits");
#}
#exit;


$dev="12".$DEV;
plot_MAP_color_scale($pdl_MAP_SMass_SFR,$pdl_fMAP_SMass_SFR,$pdl_ctable,$dev,-1.1,3.3,-10.5,-5.5,0,49,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","log(\\gS\\dSFR\\u/M\\d\\m9\\u yr\\u-1\\d pc\\u-2\\d)","% Lum. Age<200 Myr","$name",1,1);


$dev="13".$DEV;
($pdl_MAP_SMass_SFR,$pdl_fMAP_SMass_SFR)=plot_XY_color_scale($pdl_Sigma_Mass,$pdl_Sigma_SFR,$pdl_log_EW_Ha,$pdl_ctable,$dev,,-1.1,3.3,-10.5,-5.5,0,2.5,"log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","log(\\gS\\dSFR\\u/M\\d\\m9\\u yr\\u-1\\d pc\\u-2\\d)","log |EW(H\\ga)|","$name");

#$outfile=$name.".p_e.MAP_SMass_SFR.fits";
#$pdl_MAP_SMass_SFR->wfits($outfile);
$outfile=$name.".p_e.fMAP_SMass_SFR_EW.fits";
$pdl_fMAP_SMass_SFR->wfits($outfile);

$dev="14".$DEV;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_sfh_lum_age, $kernel, {Boundary => Reflect});
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Lum. Weight","R/Re","$name",$pdl_label_X,$pdl_label_Y);



$dev="14_V".$DEV;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_sfh_lum_age_V, $kernel, {Boundary => Reflect});
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Light Distribution","R/Re","$name",$pdl_label_X,$pdl_label_Y);

#print "pdl_age = $pdl_label_X\n";
$dev="14_Mass_norm".$DEV;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_Mass_Norm, $kernel, {Boundary => Reflect});
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Normalized Mass","R/Re","$name",$pdl_label_X,$pdl_label_Y);


$dev="14_Met_norm".$DEV;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_met_Norm, $kernel, {Boundary => Reflect});
#
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Normalized [Z/H]","R/Re","$name",$pdl_label_X,$pdl_label_Y);


$dev="14_Met_norm_TEST".$DEV;
my $kernel=ones(3,3)/9;
$pdl_rad_met_Norm_TEST=10**$pdl_rad_met_Norm;
my $pdl_smooth=conv2d($pdl_rad_met_Norm_TEST, $kernel, {Boundary => Reflect});
#plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Normalized [Z/H]","R/Re","$name",$pdl_label_X,$pdl_label_Y);
plot_MAP_line_color_scale_table($pdl_rad_met_Norm_TEST,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Normalized [Z/H]","R/Re","$name",$pdl_label_X,$pdl_label_Y);


$dev="14_Met_norm_Mass".$DEV;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_met_Norm_Mass, $kernel, {Boundary => Reflect});
#
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Normalized [Z/H] Mass","R/Re","$name",$pdl_label_X,$pdl_label_Y);


$dev="14_Met_norm_Mass_TEST".$DEV;
my $kernel=ones(3,3)/9;
$pdl_rad_met_Norm_Mass_TEST=10**$pdl_rad_met_Norm_Mass;
my $pdl_smooth=conv2d($pdl_rad_met_Norm_Mass_TEST, $kernel, {Boundary => Reflect});
#
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Normalized [Z/H] Mass","R/Re","$name",$pdl_label_X,$pdl_label_Y);
#plot_MAP_line_color_scale_table($pdl_rad_met_Norm,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Normalized [Z/H]","R/Re","$name",$pdl_label_X,$pdl_label_Y);




$dev="15".$DEV;
my $pdl_map=$pdl_elines->slice(":,:,($n_F_Ha)");
my $h = $pdl_elines->gethdr();
$$h{CRVAL1}=-$xc;#*0.5;
$$h{CDELT1}=1;#0.5;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc;#*0.5;
$$h{CDELT2}=1.0;#0.5;
$$h{CRPIX2}=1;
$pdl_log_map=$pdl_map->log10();
#$pdl_log_map=$pdl_log_map
$pdl_map->sethdr($h);
$V_img->sethdr($h);
#$pdl_map->sethdr($h);


plot_MAP_color_scale($pdl_map,$pdl_log_map,$pdl_ctable,$dev,0,0,0,0,-1.5,1.5,"RA (arcsec)","DEC (arcsec)","H\\ga Flux (10\\u-16\\d Erg/s/cm\\u2\\d/arcsec\\u2\\d)","$name",1,1,1);
#$pdl_map->wfits("pdl.fits");

$dev="16".$DEV;
my $pdl_sec_elines=$pdl_elines->slice(":,:,0:50");
#$pdl_sec_elines=$pdl_sec_elines*4;
my $h = $pdl_elines->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sec_elines->sethdr($h);
$pdl_log_sec_elines=$pdl_sec_elines->log10();
plot_CUBE_color_scale($pdl_sec_elines,$pdl_log_sec_elines,$pdl_ctable,$dev,0,0,0,0,-1.99,2.49,"RA (arcsec)","DEC (arcsec)","log (Flux/10\\u-16\\d Erg/s/cm\\u2\\d/arcsec\\u2\\d) $name","",0,1,1,9,6);

$dev="17".$DEV;
my $pdl_sec_velines=$pdl_elines->slice(":,:,51:101");
my $h = $pdl_elines->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sec_velines->sethdr($h);
my $vel_min=-350;
my $vel_max=+350;
$pdl_sec_velines=$pdl_sec_velines-$a_vel_Ha_cen[0]-17;#-$redshift*300000;#$a_vel_Ha_cen[0];
print "# VEL_SYS_GAS = $a_vel_Ha_cen[0] km/s\n";
plot_CUBE_color_scale($pdl_sec_elines,$pdl_sec_velines,$pdl_ctable_vel,$dev,0,0,0,0,$vel_min,$vel_max,"RA (arcsec)","DEC (arcsec)","Gas Vel. (km/s) $name","",0,1,1,9,6);


#
#
#
$test_t=1;
if ($test_t==1) {
$t_vel_Ha=$pdl_elines->slice(":,:,(96)");
$t_e_vel_Ha=$pdl_elines->slice(":,:,(300)");
$t_flux_Ha=$pdl_elines->slice(":,:,(45)");
$t_e_flux_Ha=$pdl_elines->slice(":,:,(249)");
$t_disp_Ha=$pdl_elines->slice(":,:,(147)");
$t_e_disp_Ha=$pdl_elines->slice(":,:,(351)");
$t_SN_flux_Ha=$t_flux_Ha/$t_e_flux_Ha;
$t_SN_vel_Ha=$t_flux_Ha/$t_e_vel_Ha;
$t_SN_disp_Ha=$t_flux_Ha/$t_e_disp_Ha;

$t_flux_Ha->wfits("t_flux_Ha.fits.gz");
$t_e_flux_Ha->wfits("t_e_flux_Ha.fits.gz");
$t_vel_Ha->wfits("t_vel_Ha.fits.gz");
$t_e_vel_Ha->wfits("t_e_vel_Ha.fits.gz");
$t_disp_Ha->wfits("t_disp_Ha.fits.gz");
$t_e_disp_Ha->wfits("t_e_disp_Ha.fits.gz");
$t_SN_flux_Ha->wfits("t_SN_flux_Ha.fits.gz");
$t_SN_vel_Ha->wfits("t_SN_vel_Ha.fits.gz");
$t_SN_disp_Ha->wfits("t_SN_disp_Ha.fits.gz");
}

$dev="b17".$DEV;
my $pdl_map=$pdl_sec_velines->slice(":,:,(45)");
my $h = $pdl_sec_elines->gethdr();
$$h{CRVAL1}=-$xc;
$$h{CDELT1}=1;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc;
$$h{CDELT2}=1;
$$h{CRPIX2}=1;
$pdl_map->sethdr($h);
plot_MAP_color_scale($pdl_map,$pdl_map,$pdl_ctable_vel,$dev,0,0,0,0,$vel_min,$vel_max,"RA (arcsec)","DEC (arcsec)","H\\ga Vel. (km/s)","$name",0,1,1);





#$pdl_sfh=rfits($sfh_file);
#print "# Reading sfh\n";
#my $h=$pdl_sfh->gethdr();
#$pdl_sfh_lum_age=$pdl_sfh->slice(":,:,156:193");

my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sfh_lum_age->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_sfh_lum_age, $kernel, {Boundary => Reflect});
$pdl_smooth->sethdr($h);

$dev="18".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth,$pdl_ctable,$dev,0,0,0,0,0.001,0.699,"RA (arcsec)","DEC (arcsec)","Lum. Fraction for each age  ($name)","",0,1,1,9,6);


my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sfh_lum_age_V->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_sfh_lum_age_V, $kernel, {Boundary => Reflect});
#$pdl_sfh_lum_age_V->wfits("pdl.fits");
$pdl_smooth->inplace->log10();
$pdl_smooth=$pdl_smooth-16;
$pdl_smooth->sethdr($h);
$dev="19".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth,$pdl_ctable,$dev,0,0,0,0,-19.99,-14.99,"RA (arcsec)","DEC (arcsec)","log(Flux/Erg/s/cm\\u2\\d) \@age,V-band  ($name)","",0,1,1,8,5);

$pdl_rad_sfh_lum_age_V=radial_sum_ring($pdl_sfh_lum_age_V,$pdl_r,0.1,3.5);
($n_ages,$n_r)=$pdl_rad_sfh_lum_age_V->dims;
my $h = {NAXIS=>2, NAXIS1=>$n_ages, NAXIS2=>$n_r, COMMENT=>"FITs header"};
my $head="NAME";
$$h{$head}="Radial SFH";
$pdl_rad_sfh_lum_age_V->sethdr($h);

$dev="10_V".$DEV;
my $r_max=0.1*$n_r;
$pdl_label_Y=0.1*pdl([0..$n_r]);
$pdl_label_X=pdl(@age);
$pdl_label_X->inplace->log10;
$pdl_label_X=rint(10*(9+$pdl_label_X))/10;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_sfh_lum_age_V, $kernel, {Boundary => Reflect});

#$pdl_smooth->inplace->log10;

plot_MAP_color_scale_table($pdl_smooth,$pdl_smooth,$pdl_ctable,$dev,-0.001,0.035,"log(age/yr)","R/Re","Light Distribution","$name",0,1,$pdl_label_X,$pdl_label_Y);

$outfile=$name.".p_e.rad_SFH_lum_V.fits";
$pdl_rad_sfh_lum_age_V->wfits($outfile);


my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sfh_lum_age_Mass->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth2=conv2d($pdl_sfh_lum_age_Mass, $kernel, {Boundary => Reflect});
$pdl_smooth2->inplace->log10();
$pdl_smooth2=$pdl_smooth2+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2
#$pdl_smooth2=$pdl_smooth2+log10(4); # MaNGA spaxels to arcsec!
$pdl_smooth2->sethdr($h);
$dev="20".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth2,$pdl_ctable,$dev,0,0,0,0,-1.99,3.99,"RA (arcsec)","DEC (arcsec)","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d) age ($name)","",0,1,1,8,5);


my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sfh_lum_age_Mass_Time->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth2=conv2d($pdl_sfh_lum_age_Mass_Time, $kernel, {Boundary => Reflect});
$pdl_smooth2->inplace->log10();
$pdl_smooth2=$pdl_smooth2+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2
#$pdl_smooth2=$pdl_smooth2+log10(4); # MaNGA spaxels to arcsec!
$pdl_smooth2->sethdr($h);
$dev="20_Time".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth2,$pdl_ctable,$dev,0,0,0,0,-1.99,3.99,"RA (arcsec)","DEC (arcsec)","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d) \@age ($name)","",0,1,1,8,5);


$dev="10_Mass".$DEV;
my $r_max=0.1*$n_r;
$pdl_label_Y=0.1*pdl([0..$n_r]);
$pdl_label_X=pdl(@age);
$pdl_label_X->inplace->log10;
$pdl_label_X=rint(10*(9+$pdl_label_X))/10;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_sfh_lum_age_Mass, $kernel, {Boundary => Reflect});
$pdl_smooth->inplace->log10;
$pdl_smooth=$pdl_smooth-log10(($DA)**2)-6; # Msun/pc^2
plot_MAP_color_scale_table($pdl_smooth,$pdl_smooth,$pdl_ctable,$dev,-2.5,3.5,"log(age/yr)","R/Re","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","$name",0,1,$pdl_label_X,$pdl_label_Y);

$outfile=$name.".p_e.rad_SFH_lum_Mass.fits";
$pdl_smooth->wfits($outfile);

#exit;

$dev="14_Mass".$DEV;
#my $kernel=ones(3,3)/9;
#my $pdl_smooth=conv2d($pdl_rad_sfh_lum_age_Mass, $kernel, {Boundary => Reflect});
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","R/Re","$name",$pdl_label_X,$pdl_label_Y,-2.5,3.5);
#plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,-0.1,2.5,"log(age/yr)","Light Distribution","R/Re","$name",$pdl_label_X,$pdl_label_Y);





my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sfh_lum_age_Mass_C->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth3=conv2d($pdl_sfh_lum_age_Mass_C, $kernel, {Boundary => Reflect});
$pdl_smooth3->inplace->log10();
$pdl_smooth3=$pdl_smooth3+$pdl_scale->log10-log10(($DA)**2)-6; # Msun/pc^2
#$pdl_smooth3=$pdl_smooth3+log10(4); # MaNGA spaxels to arcsec!
$pdl_smooth3->sethdr($h);
$dev="21".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth3,$pdl_ctable,$dev,0,0,0,0,-1.99,3.99,"RA (arcsec)","DEC (arcsec)","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d) \\gSage/Gyr ($name)","",0,1,1,8,5);
#$pdl_smooth3->wfits("pdl_C2.fits");

my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*1.0;
$$h{CDELT1}=1.0;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*1.0;
$$h{CDELT2}=1.0;
$$h{CRPIX2}=1;
$pdl_sfh_lum_age_V_C->sethdr($h);
my $kernel=ones(2,2)/9;
my $pdl_smooth3=conv2d($pdl_sfh_lum_age_V_C, $kernel, {Boundary => Reflect});
$pdl_smooth3->inplace->log10();
$pdl_smooth3=$pdl_smooth3+$pdl_scale-16;;#+$pdl_scale->log10;#-log10(($DA)**2)-6; # Msun/pc^2
#$pdl_smooth3=$pdl_smooth3;#+log10(4); # MaNGA spaxels to arcsec!
$pdl_smooth3->sethdr($h);
$dev="22".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth3,$pdl_ctable,$dev,0,0,0,0,-19.5,-15.5,,"RA (arcsec)","DEC (arcsec)","log(Flux/Erg/s/cm\\u2\\d) \\Sage/Gyr ,V-band  ($name)","",0,1,1,8,5);
#$pdl_smooth3->wfits("pdl_C2.fits");



$dev="23".$DEV;
my $pdl_map=$pdl_ssp->slice(":,:,(0)");
my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_map->sethdr($h);
plot_MAP_color_scale($pdl_map,$pdl_map,$pdl_ctable,$dev,0,0,0,0,-0.05,0.75,"RA (arcsec)","DEC (arcsec)","V-band Flux (10\\u-16\\d Erg/s/cm\\u2\\d/arcsec\\u2\\d)","$name",0,1,1);


$dev="24".$DEV;
my $pdl_map=$pdl_ssp->slice(":,:,(1)");
my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_map->sethdr($h);
plot_MAP_color_scale($V_img,$pdl_map,$pdl_ctable,$dev,0,0,0,0,-100,250,"RA (arcsec)","DEC (arcsec)","Segmentation Index","$name",1,1,1);

$dev="25".$DEV;
my $pdl_map=$pdl_ssp->slice(":,:,(5)");
$pdl_map=$pdl_map*$img_mask_new;
my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_map->sethdr($h);
$V_img->sethdr($h);
plot_MAP_color_scale($V_img,$pdl_map,$pdl_ctable,$dev,0,0,0,0,7.1,10.2,"RA (arcsec)","DEC (arcsec)","log(age/yr)","$name",1,1,1);


$dev="26".$DEV;
my $pdl_map=$pdl_ssp->slice(":,:,(8)");
$pdl_map=$pdl_map/$img_mask_new;
my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_map->sethdr($h);
plot_MAP_color_scale($V_img,$pdl_map,$pdl_ctable,$dev,0,0,0,0,-0.75,0.22,"RA (arcsec)","DEC (arcsec)","[Z/H]","$name",1,1,1);

$dev="27".$DEV;
my $pdl_map=$pdl_ssp->slice(":,:,(11)");
$pdl_map=$pdl_map/$img_mask_new;
my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_map->sethdr($h);
plot_MAP_color_scale($V_img,$pdl_map,$pdl_ctable,$dev,0,0,0,0,0,1.6,"RA (arcsec)","DEC (arcsec)","Av (mag)","$name",1,1,1);


$dev="28".$DEV;
my $pdl_map=$pdl_ssp->slice(":,:,(15)");
$pdl_map=$pdl_map/$img_mask_new;
my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_map->sethdr($h);
plot_MAP_color_scale($V_img,$pdl_map,$pdl_ctable,$dev,0,0,0,0,0,300,"RA (arcsec)","DEC (arcsec)","\gs (km/s)","$name",1,1,1);


$dev="29".$DEV;
my $pdl_map=$pdl_ssp->slice(":,:,(13)");
$pdl_map=$pdl_map/$img_mask_new;
my $h = $pdl_ssp->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_map=$pdl_map-$redshift*300000;#$a_vel_Ha_cen[0];
$pdl_map->sethdr($h);

#$pdl_sec_velines=$pdl_sec_velines-$redshift*300000;#$a_vel_Ha_cen[0];
print "# VEL_SYS_SSP = $a_vel_ssp_cen[0] km/s\n";
plot_MAP_color_scale($V_img,$pdl_map,$pdl_ctable_vel,$dev,0,0,0,0,$vel_min,$vel_max,"RA (arcsec)","DEC (arcsec)","Stellar Vel. (km/s) $name","",1,1,1,9,6);

#$pdl_map->wfits("pdl.fits");


$dev="30".$DEV;
$pdl_rat_OIII=$pdl_F_OIII/$pdl_F_OIIIw;
$pdl_log_F_OIIIw=$pdl_F_OIIIw->log10();
$pdl_SN_OIIIw=$pdl_F_OIIIw/$pdl_eF_OIIIw;
$pdl_l_SN_OIIIw=$pdl_SN_OIIIw->log10();
plot_XY_color_scale($pdl_log_F_OIIIw,$pdl_rat_OIII,$pdl_l_SN_OIIIw,$pdl_ctable,$dev,,-2.5,2,0,6,-2,1.5,"[OIII] 4959 log(Flux/10\\u-16\\d Erg/s/cm\\u2\\d/spaxel\\u2\\d)","[OIII]5007/[OIII]4959 Flux ratio","log(S/N) [OIII] 4959","$name");


$dev="31".$DEV;
$pdl_rat_NII=$pdl_F_NII/$pdl_F_NIIw;
$pdl_log_F_NIIw=$pdl_F_NIIw->log10();
$pdl_SN_NIIw=$pdl_F_NIIw/$pdl_eF_NIIw;
$pdl_l_SN_NIIw=$pdl_SN_NIIw->log10();
plot_XY_color_scale($pdl_log_F_NIIw,$pdl_rat_NII,$pdl_l_SN_NIIw,$pdl_ctable,$dev,,-2.5,2,0,6,-2,1.5,"[NII] 6548 log(Flux/10\\u-16\\d Erg/s/cm\\u2\\d/spaxel\\u2\\d)","[NII]6583/[NII]6458 Flux ratio","log(S/N) [NII] 6548","$name");

#
# More diagnostic diagrams!
#

my @cut_x,@cut_y;
for ($i=0;$i<200;$i++) {
    $cut_x[$i]=-2+0.02*$i;
    if ($cut_x[$i]<1.5) {
        $cut_y[$i]=0.15/(($cut_x[$i]-1.5))+1.10;
    } else {
        $cut_y[$i]=-100;
    }
}
my $pdl_x_line=pdl(@cut_x);
my $pdl_y_line=pdl(@cut_y);

$dev="32".$DEV;
($pdl_MAP_OIII_OII,$pdl_fMAP_OIII_OII)=plot_XY_color_scale($pdl_OII_Hb_clean,$pdl_OIII_Hb_clean,$pdl_r,$pdl_ctable,$dev,-1,2.25,-1.3,1.3,0,3.1,"log([OII]/H\\gb)","log([OIII]/H\\gb)","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);

$outfile=$name.".p_e.MAP_OIII_OII.fits";
$pdl_MAP_OIII_OII->wfits($outfile);
$outfile=$name.".p_e.fMAP_OIII_OII_r.fits";
$pdl_fMAP_OIII_OII->wfits($outfile);


my @cut_x,@cut_y;
for ($i=0;$i<200;$i++) {
    $cut_x[$i]=-2+0.02*$i;
    if ($cut_x[$i]<0.15) {
        #$cut_y[$i]=0.15/(($cut_x[$i]-1.5))+1.10;
	$cut_y[$i]=0.61/($cut_x[$i]-0.3)+1.3;
    } else {
        $cut_y[$i]=-100;
    }
}
my $pdl_x_line=pdl(@cut_x);
my $pdl_y_line=pdl(@cut_y);

$dev="33".$DEV;
($pdl_MAP_OIII_SII,$pdl_fMAP_OIII_SII)=plot_XY_color_scale($pdl_SII_Ha_clean,$pdl_OIII_Hb_clean,$pdl_r,$pdl_ctable,$dev,-1,0.5,-1.3,1.3,0,3.1,"log([SII]/H\\ga)","log([OIII]/H\\gb)","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);

$outfile=$name.".p_e.MAP_OIII_SII.fits";
$pdl_MAP_OIII_SII->wfits($outfile);
$outfile=$name.".p_e.fMAP_OIII_SII_r.fits";
$pdl_fMAP_OIII_SII->wfits($outfile);


my @cut_x,@cut_y;
for ($i=0;$i<200;$i++) {
    $cut_x[$i]=-2+0.02*$i;

    $cut_y3[$i]=0.61/($cut_x[$i]-0.05)+1.3;


    if ($cut_x[$i]<0) {
        $cut_y3_d[$i]=0.15/(($cut_x[$i]-1.5))+1.10;
	$cut_y[$i]=$cut_y3[$i]/$cut_y3_d[$i];
    } else {
        $cut_y3_d[$i]=-100;
	$cut_y[$i]=-100;
    }




}
my $pdl_x_line=pdl(@cut_x);
my $pdl_y_line=pdl(@cut_y);

$dev="34".$DEV;
($pdl_MAP_OIII_OII_NII,$pdl_fMAP_OIII_OII_NII)=plot_XY_color_scale($pdl_NII_Ha_clean,$pdl_OIII_OII_clean,$pdl_r,$pdl_ctable,$dev,-1,0.5,-1.3,1.3,0,3.1,"log([NII]/H\\ga)","log([OIII]/[OII])","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);

$outfile=$name.".p_e.MAP_OIII_OII_NII.fits";
$pdl_MAP_OIII_OII_NII->wfits($outfile);
$outfile=$name.".p_e.fMAP_OIII_OII_NII_r.fits";
$pdl_fMAP_OIII_OII_NII->wfits($outfile);


my @cut_x,@cut_y;
for ($i=0;$i<200;$i++) {
    $cut_x[$i]=-2+0.02*$i;

    $cut_y3[$i]=0.61/($cut_x[$i]-0.3)+1.3;

    if ($cut_x[$i]<0.15) {
        $cut_y3_d[$i]=0.15/(($cut_x[$i]-1.5))+1.10;
	$cut_y[$i]=$cut_y3[$i]/$cut_y3_d[$i];
    } else {
	$cut_y3_d[$i]=-100;
	$cut_y[$i]=-100;
    }


}
my $pdl_x_line=pdl(@cut_x);
my $pdl_y_line=pdl(@cut_y);

$dev="35".$DEV;
($pdl_MAP_OIII_OII_SII,$pdl_fMAP_OIII_OII_SII)=plot_XY_color_scale($pdl_SII_Ha_clean,$pdl_OIII_OII_clean,$pdl_r,$pdl_ctable,$dev,-1,0.5,-1.3,1.3,0,3.1,"log([SII]/H\\ga)","log([OIII]/[OII])","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);

$outfile=$name.".p_e.MAP_OIII_OII_SII.fits";
$pdl_MAP_OIII_OII_SII->wfits($outfile);
$outfile=$name.".p_e.fMAP_OIII_OII_SII_r.fits";
$pdl_fMAP_OIII_OII_SII->wfits($outfile);


my @cut_x,@cut_y;
for ($i=0;$i<100;$i++) {
    $cut_x[$i]=7.5+0.02*$i;
    $cut_y[$i]=$cut_x[$i];
}
my $pdl_x_line=pdl(@cut_x);
my $pdl_y_line=pdl(@cut_y);
$dev="36".$DEV;
($pdl_MAP_OH_O3N2_N2,$pdl_fMAP_OH_O3N2_N2)=plot_XY_color_scale($pdl_OH_O3N2,$pdl_OH_N2,$pdl_r,$pdl_ctable,$dev,7.9,8.9,7.9,8.9,0,3.1,"12+log(O/H) O3N2","12+log(O/H) N2","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);
$outfile=$name.".p_e.MAP_OH_O3N2_N2.fits";
$pdl_MAP_OH_O3N2_N2->wfits($outfile);
$outfile=$name.".p_e.fMAP_OH_O3N2_N2_r.fits";
$pdl_fMAP_OH_O3N2_N2->wfits($outfile);

$dev="37".$DEV;
($pdl_MAP_OH_O3N2_ONS,$pdl_fMAP_OH_O3N2_ONS)=plot_XY_color_scale($pdl_OH_O3N2,$pdl_OH_ONS,$pdl_r,$pdl_ctable,$dev,7.9,8.9,7.9,8.9,0,3.1,"12+log(O/H) O3N2","12+log(O/H) ONS","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);
$outfile=$name.".p_e.MAP_OH_O3N2_ONS.fits";
$pdl_MAP_OH_O3N2_ONS->wfits($outfile);
$outfile=$name.".p_e.fMAP_OH_O3N2_ONS_r.fits";
$pdl_fMAP_OH_O3N2_ONS->wfits($outfile);

$dev="38".$DEV;
($pdl_MAP_OH_O3N2_R23,$pdl_fMAP_OH_O3N2_R23)=plot_XY_color_scale($pdl_OH_O3N2,$pdl_OH_R23,$pdl_r,$pdl_ctable,$dev,7.9,8.9,7.9,8.9,0,3.1,"12+log(O/H) O3N2","12+log(O/H) R23","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);
$outfile=$name.".p_e.MAP_OH_O3N2_R23.fits";
$pdl_MAP_OH_O3N2_ONS->wfits($outfile);
$outfile=$name.".p_e.fMAP_OH_O3N2_R23_r.fits";
$pdl_fMAP_OH_O3N2_ONS->wfits($outfile);



$pdl_rat_SII=($pdl_F_SII6717/$pdl_F_SII6731)*($pdl_mask_e_Ha)*($pdl_mask_e_Ha)*($pdl_mask_e_SII);
#$pdl_rat_SII->inplace->log10;
$pdl_rat_SII->inplace->setnantobad;
$pdl_nEL=abs(10000*($pdl_rat_SII/1.49-1)/(3.77-(12.8/1.49)*$pdl_rat_SII));
# $nEL{$id}=abs(10000*($rat_SII{$id}/1.49-1)/(3.77-(12.8/1.49)*$rat_SII{$id}));


$dev="39".$DEV;
my $h = $pdl_elines->gethdr();
$$h{CRVAL1}=-$xc;#*0.5;
$$h{CDELT1}=1;#0.5;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc;#*0.5;
$$h{CDELT2}=1.0;#0.5;
$$h{CRPIX2}=1;
$pdl_nEL->sethdr($h);
$pdl_F_Ha->sethdr($h);
#$pdl_F_Ha->inplace->hclip(0,1e12);
$pdl_l_nEL=$pdl_nEL->log10;
plot_MAP_color_scale($pdl_F_Ha,$pdl_l_nEL,$pdl_ctable,$dev,0,0,0,0,-1,3.5,"RA (arcsec)","DEC (arcsec)","Electron Density (cm\\u-3\\d)","$name",1,1,1);
#$pdl_nEL->wfits("pdl.fits");


$dev="40".$DEV;
my $h = $pdl_elines->gethdr();
$$h{CRVAL1}=-$xc*$spax_scale;
$$h{CDELT1}=$spax_scale;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*$spax_scale;
$$h{CDELT2}=$spax_scale;
$$h{CRPIX2}=1;
$pdl_Av=$pdl_Av*($pdl_mask_e_Ha)*$pdl_mask_e_Hb;
$pdl_Av->sethdr($h);
plot_MAP_color_scale($pdl_F_Ha,$pdl_Av,$pdl_ctable,$dev,0,0,0,0,0,2.6,"RA (arcsec)","DEC (arcsec)","Av gas (mag)","$name",1,1,1);


my @cut_x,@cut_y;
for ($i=0;$i<100;$i++) {
    $cut_x[$i]=0+0.025*$i;
    $cut_y[$i]=0.45*$cut_x[$i];
}
my $pdl_x_line=pdl(@cut_x);
my $pdl_y_line=pdl(@cut_y);

($pdl_e_Av_ssp_clean,$pdl_mask_e_Av,$stats_e_Av)=create_mask($pdl_e_Av_ssp/$pdl_Av_ssp,-1,3);
$pdl_Av_ssp_clean=$pdl_Av_ssp*$pdl_mask_e_Av*$img_mask_new;
$pdl_Av_ssp_clean->sethdr($h);


$dev="41".$DEV;
($pdl_MAP_Av_gas_ssp,$pdl_fMAP_Av_gas_ssp)=plot_XY_color_scale($pdl_Av,$pdl_Av_ssp_clean,$pdl_r,$pdl_ctable,$dev,,0,2.6,0,2.6,0,3.1,"Av gas (mag)","Av stellar (mag)","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);
$outfile=$name.".p_e.MAP_Av_gas_ssp.fits";
$pdl_MAP_OH_O3N2_ONS->wfits($outfile);
$outfile=$name.".p_e.fMAP_Av_gas_ssp_r.fits";
$pdl_fMAP_OH_O3N2_ONS->wfits($outfile);

my @cut_x,@cut_y;
for ($i=0;$i<100;$i++) {
    $cut_x[$i]=-4+0.04*$i;
    $cut_y[$i]=$cut_x[$i];
}
my $pdl_x_line=pdl(@cut_x);
my $pdl_y_line=pdl(@cut_y);
$dev="42".$DEV;
#$pdl_log_U_P
$pdl_log_U_iter->inplace->setvaltobad(0);
$pdl_log_U_P->inplace->setvaltobad(0);

plot_XY_color_scale($pdl_log_U_P,$pdl_log_U_iter,$pdl_r,$pdl_ctable,$dev,,-5,-0.5,-5,-0.5,0,3.1,"log(u) [OIII]/[OII]","log(u) Izotov","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);


$dev="43".$DEV;
($pdl_MAP_SMass_Av_gas,$pdl_fMAP_SMass_Av_gas)=plot_XY_color_scale($pdl_Av_org,$pdl_Sigma_Mass,$pdl_r,$pdl_ctable,$dev,,0,2.6,-1.1,3.3,0,3.1,"Av gas (mag)","log(\\gS\\dMass\\u/M\\d\\m9\\u pc\\u-2\\d)","R/Re","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);
$outfile=$name.".p_e.MAP_SMass_Av_gas.fits";
$pdl_MAP_SMass_Av_gas->wfits($outfile);
$outfile=$name.".p_e.fMAP_SMass_Av_gas_r.fits";
$pdl_fMAP_SMass_Av_gas->wfits($outfile);

$dev="44".$DEV;
($pdl_MAP_S_SFR_Av_gas,$pdl_fMAP_S_SFR_Av_gas)=plot_XY_color_scale($pdl_Sigma_SFR,$pdl_Av_org,$pdl_log_EW_Ha,$pdl_ctable,$dev,,-9.5,-7.5,0,2.6,0,2.5,"log(\\gS\\dSFR\\u/M\\d\\m9\\u yr\\u-1\\d pc\\u-2\\d)","Av gas (mag)","log10 |EW(H\\ga)|","$name",0,0,1,1,$pdl_x_line,$pdl_y_line,1);
$outfile=$name.".p_e.MAP_S_SFR_Av_gas.fits";
$pdl_MAP_S_SFR_Av_gas->wfits($outfile);
$outfile=$name.".p_e.fMAP_S_SFR_Av_gas_EW_Ha.fits";
$pdl_fMAP_S_SFR_Av_gas->wfits($outfile);

#$pdl_Av->wfits("pdl.fits");

$dev="45".$DEV;
my $r_max=0.1*$n_r;
$pdl_label_Y=0.1*pdl([0..$n_r]);
$pdl_label_X=pdl(@age);
$pdl_label_X->inplace->log10;
$pdl_label_X=rint(10*(9+$pdl_label_X))/10;
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_rad_sfh_met_C, $kernel, {Boundary => Reflect});
#$pdl_smooth->wfits("pdl.fits");
$pdl_smooth->inplace->log10;
$pdl_smooth=$pdl_smooth-log10(0.02);
#$pdl_smooth->wfits("pdl.fits");
#$pdl_smooth=$pdl_smooth-log10(($DA)**2)-6; # Msun/pc^2
plot_MAP_color_scale_table($pdl_smooth,$pdl_smooth,$pdl_ctable,$dev,-0.75,0.22,"log(age/yr)","R/Re","[Z/H]","$name",0,1,$pdl_label_X,$pdl_label_Y);

$outfile=$name.".p_e.rad_SFH_Met.fits";
$pdl_smooth->wfits($outfile);

$dev="45_rad".$DEV;
plot_MAP_line_color_scale_table($pdl_smooth,$pdl_label_Y,$pdl_ctable,$dev,0,2.5,"log(age/yr)","[Z/H]","R/Re","$name",$pdl_label_X,$pdl_label_Y,-0.75,0.22);


my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*0.5;
$$h{CDELT1}=0.5;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*0.5;
$$h{CDELT2}=0.5;
$$h{CRPIX2}=1;
$pdl_sfh_met_C->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_sfh_met_C, $kernel, {Boundary => Reflect});
$pdl_smooth->inplace->log10();
$pdl_smooth=$pdl_smooth-log10(0.02);
$pdl_smooth->sethdr($h);
$dev="46".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth,$pdl_ctable,$dev,0,0,0,0,-0.75,0.22,"RA (arcsec)","DEC (arcsec)","[Z/H] LW ($name)","",0,1,1,8,5);


my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*0.5;
$$h{CDELT1}=0.5;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*0.5;
$$h{CDELT2}=0.5;
$$h{CRPIX2}=1;
$pdl_sfh_met_Mass_C->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_sfh_met_Mass_C, $kernel, {Boundary => Reflect});
$pdl_smooth->inplace->log10();
$pdl_smooth=$pdl_smooth-log10(0.02);
$pdl_smooth->sethdr($h);
$dev="46Mass".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth,$pdl_ctable,$dev,0,0,0,0,-0.75,0.22,"RA (arcsec)","DEC (arcsec)","[Z/H] MW ($name)","",0,1,1,8,5);

my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*0.5;
$$h{CDELT1}=0.5;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*0.5;
$$h{CDELT2}=0.5;
$$h{CRPIX2}=1;
$pdl_sfh_met->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_sfh_met, $kernel, {Boundary => Reflect});
$pdl_smooth->inplace->log10();
$pdl_smooth=$pdl_smooth-log10(0.02);
$pdl_smooth->sethdr($h);
$dev="46_age".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth,$pdl_ctable,$dev,0,0,0,0,-0.75,0.22,"RA (arcsec)","DEC (arcsec)","[Z/H] ($name) \@age","",0,1,1,8,5);


my $h = $pdl_sfh_lum_age->gethdr();
$$h{CRVAL1}=-$xc*0.5;
$$h{CDELT1}=0.5;
$$h{CRPIX1}=1;
$$h{CRVAL2}=-$yc*0.5;
$$h{CDELT2}=0.5;
$$h{CRPIX2}=1;
$pdl_sfh_met_Mass->sethdr($h);
my $kernel=ones(3,3)/9;
my $pdl_smooth=conv2d($pdl_sfh_met_Mass, $kernel, {Boundary => Reflect});
$pdl_smooth->inplace->log10();
$pdl_smooth=$pdl_smooth-log10(0.02);
$pdl_smooth->sethdr($h);
$dev="46_age_Mass".$DEV;
plot_CUBE_color_scale($pdl_sfh_lum_age,$pdl_smooth,$pdl_ctable,$dev,0,0,0,0,-0.75,0.22,"RA (arcsec)","DEC (arcsec)","[Z/H] ($name) \@age","",0,1,1,8,5);




#	$log_U_iter=$log_q-log10(3e10);
#print "$pdl_r\n";
#print "$pdl_OH_O3N2\n";

my @tmp_X;
my @tmp_Y;
my $n_X_Y=0;
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	my $val_X=$pdl_r->at($i,$j);
	my $val_Y=$pdl_OH_O3N2->at($i,$j);
	if (($val_X>0.3)&&($val_X<2.1)&&($val_Y>0)) {
	    $tmp_X[$n_X_Y]=$val_X;
	    $tmp_Y[$n_X_Y]=$val_Y;
	    $n_X_Y++;
	}
    }
}

if ($n_X_Y>0) {
    ($pdl_r_OH_O3N2_mod,$coeff_r_OH_O3N2) = fitpoly1d(pdl(@tmp_X),pdl(@tmp_Y),2);
    $OH_O3N2_Re=$coeff_r_OH_O3N2->at(0)+$coeff_r_OH_O3N2->at(1);
    $alpha_OH_O3N2_Re=$coeff_r_OH_O3N2->at(1);
} else {
    $OH_O3N2_Re="nan";
    $alpha_OH_O3N2_Re="nan";
}
#print "# weighted Av = $Av_w  (mag)\n";	    
$z_gas=$a_vel_Ha_cen[0]/$speed_of_light;
$z_stars=$a_vel_ssp_cen[0]/$speed_of_light;

$pdl_Av_ssp_clean->setvaltobad(0);
@Av_ssp_stats=stats($pdl_Av_ssp);

print "# Av_ssp = $Av_ssp_stats[0]\n";
open(OUT,">proc_elines.csv");
print OUT "$name,$log_Mass_corr,$lSFR,$OH_O3N2,$e_OH_O3N2,$a_ion_cen[3],$frac_area_GAS,$frac_area_SFR_pure,$Av_w,$OH_O3N2_Re,$alpha_OH_O3N2_Re,$a_vel_Ha_cen[0],$a_vel_ssp_cen[0],$z_gas,$z_stars,$FoV,$nx,$ny,$Re_kpc,$mag_g,$log_Mass,$Av_ssp_stats[2]\n";
close(OUT);

#$frac_area_SFR_pure=$subm_SFR_pure/(6*($nx/2)*sqrt(($nx/2)**2-($nx/4)##########################################
# IONIZATION 
# -1 -> No class **2));
#$frac_area_GAS=$sum_GAS/(6*($nx/2)*sqrt(($nx/2)**2-($nx/4)**2));

#print "# CENTRAL IONIZATION = $a_ion_cen[3]\n";

exit;

sub create_mask {
    my $pdl_input=$_[0];
    my $min_val=$_[1];
    my $max_val=$_[2];
    my $pdl_output=$pdl_input->clip($min_val,$max_val);
    $pdl_output->inplace->setvaltobad($min_val);
    $pdl_output->inplace->setvaltobad($max_val);
    my $pdl_mask=$pdl_output/$pdl_output;
    $pdl_output->inplace->setnantobad;
    my @stats=stats($pdl_output);
    $pdl_mask->inplace->setbadtoval(0);
    my $pdl_stats=pdl(@stats);
    return ($pdl_output,$pdl_mask,$pdl_stats);
}


sub plot_diag {
    my $pdl_X=$_[0];
    my $pdl_Y=$_[1];    
    my $color=$_[2];
    my $dev=$_[3];    
    my $pdl_color;
    my $i,$j,$k;
    my @stats_color;
    my ($nx,$ny)=$pdl_X->dims;
    
    
    pgbegin(0,$dev,1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsvp(0.1,0.75,0.15,0.98);
    pgswin(-1.2,0.6,-1.2,1.2);


    pglabel("log([NII]\\gl6583/H\\ga)","log([OIII]\\gl5007/H\\gb)","");
    pgsci(1);
    my $sym=17;
    pgsci(1);
    pgsch(2.1);
    my $n=0;
    my $n0=0;
    my $nALL=0;
    my $nGOOD=0;
    for ($i=0;$i<$nx;$i++) {
	for ($j=0;$j<$ny;$j++) {
	my $X=$pdl_X->at($i,$j);
	my $Y=$pdl_Y->at($i,$j);
	my $c=$color;
	if (($X ne "BAD")&&($Y ne "BAD")) {
	    pgsci($c);
	    pgsch(1.1);
	    pgpoint(1,[$X],[$Y],17);
	}
	}
    }
    pgsci(1);
    pgsch(1.2);
    my @cut_x,@cut_y,@cut_y2,@cut_y3,@cut_y4;
    for ($i=0;$i<200;$i++) {
	$cut_x[$i]=-2+0.02*$i;
	$cut_y[$i]=-0.7+0.2-3.67*$cut_x[$i];
	$cut_y2[$i]=-1.7+0.5-3.67*$cut_x[$i];
	
	$cut_y3[$i]=0.61/($cut_x[$i]-0.05)+1.3;
	$cut_y4[$i]=0.61/($cut_x[$i]-0.47)+1.19;
	if ($cut_y3[$i]>1.1) {
	    $cut_y3[$i]=-10;
	}
	if ($cut_y4[$i]>1.1) {
	    $cut_y4[$i]=-10;
	}
	
    }
    pgsls(1);
    pgslw(3);
    pgsci(2);
    pgline($#cut_x,\@cut_x,\@cut_y3);
    pgsls(2);
    pgslw(3);
    pgsci(4);
    pgline($#cut_x,\@cut_x,\@cut_y4);
    
    pgsch(1.7);
    pgsls(1);
    pgsci(1);
    pgptxt(-1.1,1,0,0,"$name");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);
    pgclos();
    pgend();
}

sub plot_diag_color {
    my $pdl_X=$_[0];
    my $pdl_Y=$_[1];    
    my $pdl_color=$_[2];
    my $pdl_ctab=$_[3];
    my $dev=$_[4];    
    my $label=$_[5];
    my $label2=$_[6];
    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    my ($nx,$ny)=$pdl_X->dims;
    
    
    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1.0; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);


    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsvp(0.1,0.75,0.15,0.98);
    pgswin(-1.2,0.6,-1.2,1.2);
    pglabel("log([NII]\\gl6583/H\\ga)","log([OIII]\\gl5007/H\\gb)","");
    pgsci(1);
    my $sym=17;
    pgsci(1);
    pgsch(2.1);
    my $n=0;
    my $n0=0;
    my $nALL=0;
    my $nGOOD=0;
    my $pdl_c=50+50*$pdl_color/$stats_color[0];
    for ($i=0;$i<$nx;$i++) {
	for ($j=0;$j<$ny;$j++) {
	my $X=$pdl_X->at($i,$j);
	my $Y=$pdl_Y->at($i,$j);
	my $c;
	if (($X ne "BAD")&&($Y ne "BAD")) {
	    $c=int($pdl_c->at($i,$j));#50+int(40*$c_val/$stats_color[0]);
	    if ($c<51) {
		$c=51;
	    }
	    if ($c>99) {
		$c=99;
	    }
#	    if (($j>20)&&($j<30)) {
#		print "$c\n"; <stdin>;
#	    }
	    pgsci($c);
	    pgsch(0.9);
	    pgpoint(1,[$X],[$Y],17);
	}
	}
    }
    pgsci(1);
    pgsch(1.2);
    my @cut_x,@cut_y,@cut_y2,@cut_y3,@cut_y4;
    for ($i=0;$i<200;$i++) {
	$cut_x[$i]=-2+0.02*$i;
	$cut_y[$i]=-0.7+0.2-3.67*$cut_x[$i];
	$cut_y2[$i]=-1.7+0.5-3.67*$cut_x[$i];
	
	$cut_y3[$i]=0.61/($cut_x[$i]-0.05)+1.3;
	$cut_y4[$i]=0.61/($cut_x[$i]-0.47)+1.19;
	if ($cut_y3[$i]>1.1) {
	    $cut_y3[$i]=-10;
	}
	if ($cut_y4[$i]>1.1) {
	    $cut_y4[$i]=-10;
	}
	
    }
    pgsls(1);
    pgslw(3);
    pgsci(2);
    pgline($#cut_x,\@cut_x,\@cut_y3);
    pgsls(2);
    pgslw(3);
    pgsci(4);
    pgline($#cut_x,\@cut_x,\@cut_y4);
    
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgptxt(-1.1,1,0,0,"$label");
    pgptxt(0.4,1,0,0.5,"$label2");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);
#
# PGWEB
#

    pgsci(1);
    for ($i=60;$i<100;$i=$i+10) {
#    my $pdl_c=50+40*$pdl_color/$stats_color[0];
	my $val=apr(($i-50)/(40/$stats_color[0]));
	my $x=-1.2+($i-50)*(0.6+1.2)/50;
	my $y=1.1;
	pgptxt($x,$y,0,0.5,$val);
    }

    pgwedg("TI",0.0,1.5,50,100,"");
#pgwedg("RI",0,1.0,1,500,"");
    pgbox("SBC",0,0,"SBC",0,0);



    pgclos();
    pgend();
}

sub plot_diag_color_scale {
    my $pdl_X=$_[0];
    my $pdl_Y=$_[1];    
    my $pdl_color=$_[2];
    my $pdl_ctab=$_[3];
    my $dev=$_[4];    
    my $label=$_[5];
    my $label2=$_[6];
    my $c_min=$_[7];
    my $c_max=$_[8];
    my $NX=$_[9];
    my $NY=$_[10];
    if ($NX==0) {
	$NX=100;
    }
    if ($NY==0) {
	$NY=100;
    }


    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    my ($nx,$ny)=$pdl_X->dims;
    
    
    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1.0; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);


    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsvp(0.1,0.75,0.15,0.98);
    my $x_max=0.75;
    my $x_min=-1.59;
    my $y_max=1.2;
    my $y_min=-1.2;
    pgswin($x_min,$x_max,$y_min,$y_max);
    my $dx=($x_max-$x_min)/$NX;
    my $dy=($y_max-$y_min)/$NY;
    my $pdl_MAP=zeroes($NX,$NY);
    my $pdl_fMAP=zeroes($NX,$NY);

    pglabel("log([NII]\\gl6583/H\\ga)","log([OIII]\\gl5007/H\\gb)","");
    pgsci(1);
    my $sym=17;
    pgsci(1);
    pgsch(2.1);
    my $n=0;
    my $n0=0;
    my $nALL=0;
    my $nGOOD=0;
#    my $pdl_c=50+40*$pdl_color/$stats_color[0];

    my $pdl_c=50+50*($pdl_color-$c_min)/($c_max-$c_min);

    for ($i=0;$i<$nx;$i++) {
	for ($j=0;$j<$ny;$j++) {
	my $X=$pdl_X->at($i,$j);
	my $Y=$pdl_Y->at($i,$j);
	my $c;
	if (($X ne "BAD")&&($Y ne "BAD")) {
	    $c=int($pdl_c->at($i,$j));#50+int(40*$c_val/$stats_color[0]);
	    $val_c=$pdl_color->at($i,$j);
	    #print "$c $val_c\n"; <stdin>;
	    if ($c<51) {
		$c=51;
	    }
	    if ($c>99) {
		$c=99;
	    }
#	    if (($j>20)&&($j<30)) {
#		print "$c\n"; <stdin>;
#	    }
	    pgsci($c);
	    pgsch(0.9);
	    pgpoint(1,[$X],[$Y],17);

#
# Density maps
#
	    my $IX=int(($X-$x_min+0.5*$dx)/$dx);
	    my $IY=int(($Y-$y_min+0.5*$dy)/$dy);
	    if ($IX<1) {
		$IX=1;
	    }
	    if ($IY<1) {
		$IY=1;
	    }
	    if ($IX>$NX-2) {
		$IX=$NX-2;
	    }
	    if ($IY>$NY-2) {
		$IY=$NY-2;
	    }
	    my $ii,$jj;
	    for ($ii=$IX-1;$ii<$IX+2;$ii++) {
		for ($jj=$IY-1;$jj<$IY+2;$jj++) {
		    if (($ii==$IX)&&($jj==$IY)) {
			my $val=$pdl_MAP->at($ii,$jj);
			set($pdl_MAP,$ii,$jj,$val+0.2);
		    } else {
			my $val=$pdl_MAP->at($ii,$jj);
			set($pdl_MAP,$ii,$jj,$val+0.1);
		    }
		}
	    }
	    my $W=$pdl_color->at($i,$j);
	    if (($W ne "nan")&&($W ne "BAD")) {
		for ($ii=$IX-1;$ii<$IX+2;$ii++) {
		    for ($jj=$IY-1;$jj<$IY+2;$jj++) {
			if (($ii==$IX)&&($jj==$IY)) {
			    my $val=$pdl_fMAP->at($ii,$jj);
			    set($pdl_fMAP,$ii,$jj,$val+0.2*$W);
			} else {
			    my $val=$pdl_fMAP->at($ii,$jj);
			    set($pdl_fMAP,$ii,$jj,$val+$W*(0.1));
			}
		    }    
		}
	    }

	}
	}
    }
    pgsci(1);
    pgsch(1.2);
 #   $pdl_fMAP=$pdl_fMAP/$pdl_MAP;
    my @cut_x,@cut_y,@cut_y2,@cut_y3,@cut_y4;
    for ($i=0;$i<200;$i++) {
	$cut_x[$i]=-2+0.02*$i;
	$cut_y[$i]=-0.7+0.2-3.67*$cut_x[$i];
	$cut_y2[$i]=-1.7+0.5-3.67*$cut_x[$i];
	
	$cut_y3[$i]=0.61/($cut_x[$i]-0.05)+1.3;
	$cut_y4[$i]=0.61/($cut_x[$i]-0.47)+1.19;
	if ($cut_y3[$i]>1.1) {
	    $cut_y3[$i]=-10;
	}
	if ($cut_y4[$i]>1.1) {
	    $cut_y4[$i]=-10;
	}
	
    }
    pgsls(1);
    pgslw(3);
    pgsci(3);
    pgline($#cut_x,\@cut_x,\@cut_y3);
    pgsls(2);
    pgslw(3);
    pgsci(1);
    pgline($#cut_x,\@cut_x,\@cut_y4);
    
    pgsls(1);
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgptxt(-1.1,1,0,0,"$label");
#    pgptxt(0.4,1,0,0.7,"$label2");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);
#
# PGWEB
#


    pgsci(1);
    pgsch(1.8);
    pgwedg("RI",0.0,3,$c_min,$c_max,"$label2");
    pgsci(1);
#    pgsci(1);
#    for ($i=60;$i<100;$i=$i+10) {
#    my $pdl_c=50+40*$pdl_color/$stats_color[0];
#	my $val=apr(($i-50)/(40/$stats_color[0]));
#	my $val=$c_min+apr(($i-50)/(40/($c_max-$c_min)));
#    my $pdl_c=50+40*($pdl_color-$c_min)/($c_max-$c_min);
#	my $x=-1.2+($i-50)*(0.6+1.2)/50;
#	my $y=1.1;
#	pgptxt($x,$y,0,0.5,$val);
 #   }

 #   pgwedg("TI",0.0,1.5,50,100,"");
#pgwedg("RI",0,1.0,1,500,"");
 #   pgbox("SBC",0,0,"SBC",0,0);



    pgclos();
    pgend();
    
    $pdl_MAP->inplace->setnantobad;
    $pdl_MAP->inplace->setbadtoval(0);
    $pdl_fMAP->inplace->setnantobad;
    $pdl_fMAP->inplace->setbadtoval(0);


    my $h = {NAXIS=>2, NAXIS1=>$NX, NAXIS2=>$NY, COMMENT=>"Den MAP $x_label,$yl_abel"};
    $$h{CRPIX1}=1;
    $$h{CRVAL1}=$x_min;
    $$h{CDELT1}=$dx;
    $$h{CRPIX2}=1;
    $$h{CRVAL2}=$y_min;
    $$h{CDELT2}=$dy;
    $pdl_MAP->sethdr($h);

    $pdl_fMAP=$pdl_fMAP/$pdl_MAP;
    my $h = {NAXIS=>2, NAXIS1=>$NX, NAXIS2=>$NY, COMMENT=>"Den fMAP $x_label,$yl_abel,$z_label"};
    $$h{CRPIX1}=1;
    $$h{CRVAL1}=$x_min;
    $$h{CDELT1}=$dx;
    $$h{CRPIX2}=1;
    $$h{CRVAL2}=$y_min;
    $$h{CDELT2}=$dy;
    $pdl_fMAP->sethdr($h);

    return $pdl_MAP,$pdl_fMAP;

}

sub plot_XY_color_scale {
    my $pdl_X=$_[0];
    my $pdl_Y=$_[1];    
    my $pdl_color=$_[2];
    my $pdl_ctab=$_[3];
    my $dev=$_[4];    
    my $x_min=$_[5];
    my $x_max=$_[6];
    my $y_min=$_[7];
    my $y_max=$_[8];
    my $c_min=$_[9];
    my $c_max=$_[10];
    my $x_label=$_[11];
    my $y_label=$_[12];
    my $c_label=$_[13];
    my $label=$_[14];
    my $NX=$_[15];
    my $NY=$_[16];
    my $do_cont=$_[17];
    my $do_line=$_[18];
    my $pdl_x_line=$_[19];
    my $pdl_y_line=$_[20];
    my $do_sq=$_[21];
    if ($NX==0) {
	$NX=100;
    }
    if ($NY==0) {
	$NY=100;
    }

    #DENSITY MAPS
    my $dx=($x_max-$x_min)/$NX;
    my $dy=($y_max-$y_min)/$NY;
    my $pdl_MAP=zeroes($NX,$NY);
    my $pdl_fMAP=zeroes($NX,$NY);


    $pdl_color->inplace->setnantobad;
    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    if (($c_min==0)&&($c_max==0)) {
	$c_min=$stats_color[3]->at(0);
	$c_max=$stats_color[4]->at(0);
    }
    my ($nx,$ny)=$pdl_X->dims;

    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1.0; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);


    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
#    pgsvp(0.1,0.95,0.15,0.98);
    if ($do_sq==1) {
	pgsvp(0.1,0.8,0.15,0.9);
    } else {
	pgsvp(0.1,0.95,0.15,0.9);
    }
    pgswin($x_min,$x_max,$y_min,$y_max);
    pglabel("$x_label","$y_label","");
    pgsci(1);
    my $sym=17;
    pgsci(1);
    pgsch(2.1);
    my $n=0;
    my $n0=0;
    my $nALL=0;
    my $nGOOD=0;
    my $pdl_c=50+50*($pdl_color-$c_min)/($c_max-$c_min);

    for ($i=0;$i<$nx;$i++) {
	for ($j=0;$j<$ny;$j++) {
	my $X=$pdl_X->at($i,$j);
	my $Y=$pdl_Y->at($i,$j);
	my $c;
	if (($X ne "BAD")&&($Y ne "BAD")) {
	    $c=int($pdl_c->at($i,$j));
	    if ($c<51) {
		$c=51;
	    }
	    if ($c>99) {
		$c=99;
	    }
	    pgsci($c);
	    pgsch(0.9);
	    pgpoint(1,[$X],[$Y],17);
	    
#
# Density maps
#
	    my $IX=int(($X-$x_min+0.5*$dx)/$dx);
	    my $IY=int(($Y-$y_min+0.5*$dy)/$dy);
	    if (($IX>0)&&($IX<$NX-2)&&($IY>0)&&($IY<$NY-2)) {
		if ($IX<1) {
		    $IX=1;
		}
		if ($IY<1) {
		    $IY=1;
		}
		if ($IX>$NX-2) {
		    $IX=$NX-2;
		}
		if ($IY>$NY-2) {
		    $IY=$NY-2;
		}
		my $ii,$jj;
		for ($ii=$IX-1;$ii<$IX+2;$ii++) {
		    for ($jj=$IY-1;$jj<$IY+2;$jj++) {
			if (($ii==$IX)&&($jj==$IY)) {
			    my $val=$pdl_MAP->at($ii,$jj);
			    set($pdl_MAP,$ii,$jj,$val+0.2);
			} else {
			    my $val=$pdl_MAP->at($ii,$jj);
			set($pdl_MAP,$ii,$jj,$val+0.1);
			}
		    }
		}
		my $W=$pdl_color->at($i,$j);
		if (($W ne "nan")&&($W ne "BAD")) {
		    for ($ii=$IX-1;$ii<$IX+2;$ii++) {
			for ($jj=$IY-1;$jj<$IY+2;$jj++) {
			    if (($ii==$IX)&&($jj==$IY)) {
				my $val=$pdl_fMAP->at($ii,$jj);
				set($pdl_fMAP,$ii,$jj,$val+0.2*$W);
			    } else {
				my $val=$pdl_fMAP->at($ii,$jj);
				set($pdl_fMAP,$ii,$jj,$val+$W*(0.1));
			    }
			}    
		    }
		}
	    }
	}
	}
    }
    pgsci(1);
    pgsch(1.2);
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgptxt($x_min+0.05*($x_max-$x_min),$y_max-0.09*($y_max-$y_min),0,0,"$label");
#    pgptxt($x_max-0.05*($x_max-$x_min),$y_max-0.085*($y_max-$y_min),0,1,"$c_label");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);

#
#
#
    if ($do_cont==1) {
	my @stats_MAP=stats($pdl_MAP);
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @tr=($x_min-1*$dx,$dx,0,$y_min-1*$dy,0,$dy);
	my @levels;
	my $nlevels=5;
	my $ii;
	my @map=list($pdl_MAP);
	my $min=$stats_MAP[3]->at(0);
	my $max=$stats_MAP[4]->at(0);
	for ($ii=0;$ii<$nlevels;$ii++) {
	    $levels[$ii]=0.05*($max-$min)+($max-$min)/5*$ii;
	}

	pgcont(\@map,$NX,$NY,1,$NX,1,$NY,\@levels,$nlevels,\@tr);
	pgsci(1);
	pgslw(1);

    }


#
# Plot a line!
#
    if ($do_line==1) {
	my @X=list($pdl_x_line);
	my @Y=list($pdl_y_line);
	my $n=$#X;
	pgslw(3);
	pgline($n,\@X,\@Y);
	pgslw(1);
    }


#
# PGWEB
#

#    pgsci(1);
 #   for ($i=60;$i<90;$i=$i+5) {
#	my $val=$c_min+apr(($i-50)/(40/($c_max-$c_min)));
#	my $x=$x_min+($i-50)*($x_max-$x_min)/40;
#	my $y=$y_max-0.05*($y_max-$y_min);
#	pgptxt($x,$y,0,0.5,$val);
 #   }

  #  pgwedg("TI",0.0,1.5,50,100,"");
   # pgbox("SBC",0,0,"SBC",0,0);

    pgsci(1);
    pgsch(1.8);
    pgwedg("RI",0.0,3,$c_min,$c_max,"$c_label");
    pgsci(1);

    pgclos();
    pgend();

    $pdl_MAP->inplace->setnantobad;
    $pdl_MAP->inplace->setbadtoval(0);
    $pdl_fMAP->inplace->setnantobad;
    $pdl_fMAP->inplace->setbadtoval(0);
    
    my $h = {NAXIS=>2, NAXIS1=>$NX, NAXIS2=>$NY, COMMENT=>"Den MAP $x_label,$yl_abel"};
    $$h{CRPIX1}=1;
    $$h{CRVAL1}=$x_min;
    $$h{CDELT1}=$dx;
    $$h{CRPIX2}=1;
    $$h{CRVAL2}=$y_min;
    $$h{CDELT2}=$dy;
    $pdl_MAP->sethdr($h);

    $pdl_fMAP=$pdl_fMAP/$pdl_MAP;
    my $h = {NAXIS=>2, NAXIS1=>$NX, NAXIS2=>$NY, COMMENT=>"Den fMAP $x_label,$yl_abel,$z_label"};
    $$h{CRPIX1}=1;
    $$h{CRVAL1}=$x_min;
    $$h{CDELT1}=$dx;
    $$h{CRPIX2}=1;
    $$h{CRVAL2}=$y_min;
    $$h{CDELT2}=$dy;
    $pdl_fMAP->sethdr($h);



    return $pdl_MAP,$pdl_fMAP;
    # RETURN
}



sub plot_CUBE_color_scale {
    my $pdl_CUBE=$_[0];
    my $pdl_CUBE_color=$_[1];
    my $pdl_ctab=$_[2];
    my $dev=$_[3];    
    my $x_min=$_[4];
    my $x_max=$_[5];
    my $y_min=$_[6];
    my $y_max=$_[7];
    my $c_min=$_[8];
    my $c_max=$_[9];
    my $x_label=$_[10];
    my $y_label=$_[11];
    my $c_label=$_[12];
    my $label=$_[13];
    my $do_cont=$_[14];
    my $do_imag=$_[15];
    my $do_sq=$_[16];
    my $i_out=$_[17];
    my $j_out=$_[18];
    my $h_CUBE=$pdl_CUBE->gethdr();
    my ($NX,$NY,$NZ)=$pdl_CUBE->dims;
    my $tmp=sqrt($NZ/1.2);
    if (($i_out==0)&&($j_out==0)) {
	$j_out=int(sqrt($NZ/1.1));
	$i_out=int(1.4*$j_out);
    }
    
    


#    my $i_out=$i_out_pdl->at(0);
#    my $j_out=$j_out_pdl->at(0);

#    print "$i_out,$j_out | $tmp ($NX,$NY,$NZ)\n";

    my $nx=$i_out*$NX;
    my $ny=$j_out*$NY;
    my $pdl_MAP=zeroes($nx,$ny);
    my $pdl_color=zeroes($nx,$ny)-1e14;
    my $k=0;
    my @xl,@yl;
    for ($j=0;$j<$j_out;$j++) {
	for ($i=0;$i<$i_out;$i++) {
	    if ($k<$NZ) {
		my $nx0=$i*$NX;
		my $nx1=($i+1)*$NX-1;
		my $ny0=$j*$NY;
		my $ny1=($j+1)*$NY-1;
		my $sec_MAP=$pdl_CUBE->slice(":,:,($k)");
		my $t=$pdl_MAP->slice("$nx0:$nx1,$ny0:$ny1");
		$t .= $sec_MAP;
		
		my $sec_color=$pdl_CUBE_color->slice(":,:,($k)");
		my $t=$pdl_color->slice("$nx0:$nx1,$ny0:$ny1");
		$t .= $sec_color;

		$k++;
	    }

	}
    }
    
    my $h_MAP=$h_CUBE;
    my $crval1=$$h_MAP{CRVAL1};
    my $crpix1=$$h_MAP{CRPIX1};
    my $cdelt1=$$h_MAP{CDELT1};
    my $crval2=$$h_MAP{CRVAL2};
    my $crpix2=$$h_MAP{CRPIX2};
    my $cdelt2=$$h_MAP{CDELT2};

    $k=0;
    for ($j=0;$j<$j_out;$j++) {
	for ($i=0;$i<$i_out;$i++) {
	    if ($k<$NZ) {
		my $nx0=$i*$NX;
		my $nx1=($i+1)*$NX-1;
		my $ny0=$j*$NY;
		my $ny1=($j+1)*$NY-1;
		$xl[$k]=$crval1+$nx0*$cdelt1+0.2*$NX*$cdelt1-($nx1-$nx0)*$cdelt1*0.05;
		$yl[$k]=$crval2+$ny1*$cdelt2-($ny1-$ny0)*$cdelt2*0.05;
		$k++;
	    }

	}
    }

    my @tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
    my ($NX,$NY)=$pdl_MAP->dims();
    $pdl_color->inplace->setnantobad;
    $pdl_color->inplace->lclip($c_min);
    $pdl_color->inplace->setvaltobad($c_min);
    $pdl_color->inplace->setbadtoval(-1e12);
    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    if (($c_min==0)&&($c_max==0)) {
        $pdl_color->inplace->setvaltobad(-1e12);
        @stats_color=stats($pdl_color);
	$c_min=$stats_color[3]->at(0);
	$c_max=$stats_color[4]->at(0);
    }


    if (($x_min==0)&&($x_max==0)) {
	$x_min=$crval1;
	$x_max=$crval1+$NX*$cdelt1;
    }
    if (($y_min==0)&&($y_max==0)) {
	$y_min=$crval2;
	$y_max=$crval2+$NY*$cdelt2;
    }
    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    if ($do_sq==1) {
	pgsvp(0.1,0.9,0.15,0.95);
    } else {
	pgsvp(0.1,0.95,0.15,0.9);
    }
    pgswin($x_min,$x_max,$y_min,$y_max);
#    print "($x_min,$x_max,$y_min,$y_max);\n";
    pglabel("$x_label","$y_label","");
    pgsci(1);
    pgsci(1);
    pgsch(2.1);
    my $pdl_c=50+50*($pdl_color-$c_min)/($c_max-$c_min);

    pgsci(1);
    pgsch(1.2);
    pgsls(1);
    pgsci(1);

#
#
#
    if ($do_imag==1) {
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @map=list($pdl_color);
#	print "$c_min,$c_max\n";
	pgimag(\@map,$NX,$NY,1,$NX,1,$NY,$c_min,$c_max,\@tr);
	pgsci(1);
	pgslw(1);

    }

    
#    print "DO = $do_imag,$do_cont\n";
    if ($do_cont==1) {
	my @stats_MAP=stats($pdl_MAP);
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @levels;
	my $nlevels=5;
	my $ii;
	my @map=list($pdl_MAP);
	my $min=$stats_MAP[3]->at(0);
	my $max=$stats_MAP[4]->at(0);
	for ($ii=0;$ii<$nlevels;$ii++) {
#	    $levels[$ii]=0.5*$min+($max-$min)/5*$ii;
	    $levels[$ii]=0.05*($max-$min)+($max-$min)/5*$ii;		
	}

	pgcont(\@map,$NX,$NY,1,$NX,1,$NY,\@levels,$nlevels,\@tr);
	pgsci(1);
	pgslw(1);

    }


#    pgsci(3);
    pgptxt($x_max-0.4*($x_max-$x_min),$y_max-0.05*($y_max-$y_min),0,0,"$label");
#    pgsci(1);
    my $k=0;
    my $size=1.2-$NZ/100;
    if ($size<0.7) {
	$size=0.7;
    }
    pgsch($size);
    pgslw(3);
    for ($j=0;$j<$j_out;$j++) {
	for ($i=0;$i<$i_out;$i++) {
	    if ($k<$NZ) {
		my $NAME="NAME".$k;
		my $label=$$h_MAP{$NAME};
		$label =~ s/flux//g;
		$label =~ s/vel//g;
		$label =~ s/disp//g;
		pgsci(1);
		pgptxt($xl[$k],$yl[$k],0,0,"$label");
		pgsci(1);
		$k++;
	    }
	}
    }
    pgslw(1);
#    print "K=$k | $NZ\n";

    pgsch(1.2);
    pgsls(1);
    pgsci(1);
#    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);




#
# PGWEB
#

    pgsci(1);
#    for ($i=60;$i<90;$i=$i+5) {
#	my $val=$c_min+apr(($i-50)/(40/($c_max-$c_min)));
#	my $x=$x_min+($i-50)*($x_max-$x_min)/40;
#	my $y=$y_max-0.05*($y_max-$y_min);
#	pgptxt($x,$y,0,0.5,$val);
#    }

    pgsci(1);
#    pgwedg("TI",0.0,1.5,50,100,"");    
    pgsch(1.8);
    pgwedg("RI",0.0,3,$c_min,$c_max,"$c_label");
    pgsci(1);
#    pgbox("SBC",0,0,"SBC",0,0);



    pgclos();
    pgend();

    # RETURN
}


sub plot_MAP_color_scale {
    my $pdl_MAP=$_[0];
    my $pdl_color=$_[1];
    my $pdl_ctab=$_[2];
    my $dev=$_[3];    
    my $x_min=$_[4];
    my $x_max=$_[5];
    my $y_min=$_[6];
    my $y_max=$_[7];
    my $c_min=$_[8];
    my $c_max=$_[9];
    my $x_label=$_[10];
    my $y_label=$_[11];
    my $c_label=$_[12];
    my $label=$_[13];
    my $do_cont=$_[14];
    my $do_imag=$_[15];
    my $do_sq=$_[16];
    my $h_MAP=$pdl_MAP->gethdr();
    my $crval1=$$h_MAP{CRVAL1};
    my $crpix1=$$h_MAP{CRPIX1};
    my $cdelt1=$$h_MAP{CDELT1};
    my $crval2=$$h_MAP{CRVAL2};
    my $crpix2=$$h_MAP{CRPIX2};
    my $cdelt2=$$h_MAP{CDELT2};
    my @tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
    my ($NX,$NY)=$pdl_MAP->dims();
    $pdl_color->inplace->setnantobad;
    $pdl_color->inplace->lclip($c_min);
    $pdl_color->inplace->setvaltobad($c_min);
    $pdl_color->inplace->setbadtoval(-1e12);
#    $pdl_color->wfits("pdl2.fits");
    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    #if (($c_min==0)&&($c_max==0)) {
#	$c_min=$stats_color[3]->at(0);
#	$c_max=$stats_color[4]->at(0);
 #   }

    if (($c_min==0)&&($c_max==0)) {
	$pdl_color->inplace->setvaltobad(-1e12);
        @stats_color=stats($pdl_color);
	$c_min=$stats_color[3]->at(0);
	$c_max=$stats_color[4]->at(0);
    }
    
    
    if (($x_min==0)&&($x_max==0)) {
	$x_min=$crval1;
	$x_max=$crval1+$NX*$cdelt1;
    }
    if (($y_min==0)&&($y_max==0)) {
	$y_min=$crval2;
	$y_max=$crval2+$NX*$cdelt2;
    }


    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    if ($do_sq==1) {
	pgsvp(0.1,0.8,0.15,0.95);
    } else {
	pgsvp(0.1,0.95,0.15,0.9);
    }
    pgswin($x_min,$x_max,$y_min,$y_max);
    pglabel("$x_label","$y_label","");
    pgsci(1);
    pgsci(1);
    pgsch(2.1);
    my $pdl_c=50+50*($pdl_color-$c_min)/($c_max-$c_min);

    pgsci(1);
    pgsch(1.2);
    pgsls(1);
    pgsci(1);

#
#
#
    if ($do_imag==1) {
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @map=list($pdl_color);
#	print "$c_min,$c_max\n";
	pgimag(\@map,$NX,$NY,1,$NX,1,$NY,$c_min,$c_max,\@tr);
	pgsci(1);
	pgslw(1);

    }

    
#    print "DO = $do_imag,$do_cont\n";
    if ($do_cont==1) {
	$pdl_MAP=$pdl_MAP->inplace->setvaltobad(0);
	$pdl_MAP=$pdl_MAP->inplace->setvaltobad(-1000000000000);
	$pdl_MAP=$pdl_MAP->inplace->setnantobad;
	my @stats_MAP=stats($pdl_MAP);
#	print "@stats_MAP\n";
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @levels;
	my $nlevels=5;
	my $ii;
	my @map=list($pdl_MAP);
	my $min=$stats_MAP[3]->at(0);
	my $max=$stats_MAP[4]->at(0);
	for ($ii=0;$ii<$nlevels;$ii++) {
#	    $levels[$ii]=0.5*$min+($max-$min)/5*$ii;
	    $levels[$ii]=0.05*($max-$min)+($max-$min)/5*$ii;
#	    print "$ii $levels[$ii] $min $max\n";

	}

	pgcont(\@map,$NX,$NY,1,$NX,1,$NY,\@levels,$nlevels,\@tr);
	pgsci(1);
	pgslw(1);

    }



    pgptxt($x_min+0.05*($x_max-$x_min),$y_max-0.09*($y_max-$y_min),0,0,"$label");
#    pgptxt($x_max-0.05*($x_max-$x_min),$y_max-0.085*($y_max-$y_min),0,1,"$c_label");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);




#
# PGWEB
#

    pgsci(1);
#    for ($i=60;$i<90;$i=$i+5) {
#	my $val=$c_min+apr(($i-50)/(40/($c_max-$c_min)));
#	my $x=$x_min+($i-50)*($x_max-$x_min)/40;
#	my $y=$y_max-0.05*($y_max-$y_min);
#	pgptxt($x,$y,0,0.5,$val);
#    }

    pgsci(1);
#    pgwedg("TI",0.0,1.5,50,100,"");    
    pgsch(1.8);
    pgwedg("RI",0.0,3,$c_min,$c_max,"$c_label");
    pgsci(1);
#    pgbox("SBC",0,0,"SBC",0,0);



    pgclos();
    pgend();

    # RETURN
}

sub plot_MAP_color_scale_table {
    my $pdl_MAP=$_[0];
    my $pdl_color=$_[1];
    my $pdl_ctab=$_[2];
    my $dev=$_[3];    
    my $c_min=$_[4];
    my $c_max=$_[5];
    my $x_label=$_[6];
    my $y_label=$_[7];
    my $c_label=$_[8];
    my $label=$_[9];
    my $do_cont=$_[10];
    my $do_imag=$_[11];
    my $pdl_table_X=$_[12];
    my $pdl_table_Y=$_[13];
    my $h_MAP=$pdl_MAP->gethdr();
    my $crval1=0;
    my $crpix1=1;
    my $cdelt1=1;
    my $crval2=0;
    my $crpix2=1;
    my $cdelt2=1;
    my @tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
    my ($NX,$NY)=$pdl_MAP->dims();

    my $x_min=0;
    my $x_max=$NX;
    my $y_min=0;
    my $y_max=$NY;


    $pdl_color->inplace->setnantobad;
    $pdl_color->inplace->lclip($c_min);
    $pdl_color->inplace->setvaltobad($c_min);
    $pdl_color->inplace->setbadtoval(-1e12);
#    $pdl_color->wfits("pdl2.fits");
    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    if (($c_min==0)&&($c_max==0)) {
	$c_min=$stats_color[3]->at(0);
	$c_max=$stats_color[4]->at(0);
    }

#    print "C = $c_min,$c_max\n";
    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1.0; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsvp(0.1,0.95,0.15,0.9);
    pgswin(0.5,$NX+0.5,0.5,$NY+0.5);
    pglabel("$x_label","$y_label","");
    pgsci(1);
    pgsci(1);
    pgsch(2.1);
    my $pdl_c=50+50*($pdl_color-$c_min)/($c_max-$c_min);

    pgsci(1);
    pgsch(1.2);
    pgsls(1);
    pgsci(1);

#
#
#
    pgsitf(0);
    if ($do_imag==1) {
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @map=list($pdl_color);
#	print "$c_min,$c_max\n";
	pgimag(\@map,$NX,$NY,1,$NX,1,$NY,$c_min,$c_max,\@tr);
	pgsci(1);
	pgslw(1);

    }

    
#    print "DO = $do_imag,$do_cont\n";
    if ($do_cont==1) {
	my @stats_MAP=stats($pdl_MAP);
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @levels;
	my $nlevels=5;
	my $ii;
	my @map=list($pdl_MAP);
	my $min=$stats_MAP[3]->at(0);
	my $max=$stats_MAP[4]->at(0);
	for ($ii=0;$ii<$nlevels;$ii++) {
#	    $levels[$ii]=0.5*$min+($max-$min)/5*$ii;
	    $levels[$ii]=0.05*($max-$min)+($max-$min)/5*$ii;
	}

	pgcont(\@map,$NX,$NY,1,$NX,1,$NY,\@levels,$nlevels,\@tr);
	pgsci(1);
	pgslw(1);

    }



    pgptxt($x_min+0.05*($x_max-$x_min),$y_max-0.09*($y_max-$y_min),0,0,"$label");
#    pgptxt($x_max-0.05*($x_max-$x_min),$y_max-0.085*($y_max-$y_min),0,1,"$c_label");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
#    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);
    pgbox("BC",0,0,"BC",0,0);

    for ($i=1;$i<$NX;$i=$i+2) {
	my $age=$pdl_table_X->at($i);
	pgptxt($i+1,-0.025*$NY,0,0.5,"$age");
	pgline(2,[$i+1,$i+1],[0,0.05*$NY]);
    }



#

    for ($i=1;$i<$NY;$i=$i+5) {
	$R=apr($i/$Re);
	my $R=$pdl_table_Y->at($i);
	pgptxt(-0.01*$NX,$i-0.5,0,0.5,"$R");
	pgline(2,[0,0.035*$NY],[$i,$i]);
    }
#    print "PASO*****\n";




#
# PGWEB
#

    pgsci(1);
#    for ($i=60;$i<90;$i=$i+5) {
#	my $val=$c_min+apr(($i-50)/(40/($c_max-$c_min)));
#	my $x=$x_min+($i-50)*($x_max-$x_min)/40;
#	my $y=$y_max-0.05*($y_max-$y_min);
#	pgptxt($x,$y,0,0.5,$val);
#    }

    pgsci(1);
#    pgwedg("TI",0.0,1.5,50,100,"");    
    pgsch(1.8);
    pgwedg("RI",0.0,3,$c_min,$c_max,"$c_label");
    pgsci(1);
#    pgbox("SBC",0,0,"SBC",0,0);



    pgclos();
    pgend();

    # RETURN
}


#
# Radial plots!
#
sub plot_MAP_line_color_scale_table {
    my $pdl_MAP=$_[0];
    my $pdl_color=$_[1];
    my $pdl_ctab=$_[2];
    my $dev=$_[3];    
    my $c_min=$_[4];
    my $c_max=$_[5];
    my $x_label=$_[6];
    my $y_label=$_[7];
    my $c_label=$_[8];
    my $label=$_[9];
    my $pdl_table_X=$_[10];
    my $pdl_table_Y=$_[11];
 my $Y_MIN=$_[12];
    my $Y_MAX=$_[13];

    my $h_MAP=$pdl_MAP->gethdr();
    my $crval1=0;
    my $crpix1=1;
    my $cdelt1=1;
    my $crval2=0;
    my $crpix2=1;
    my $cdelt2=1;
    my @tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
    my ($NX,$NY)=$pdl_MAP->dims();
    my @stats_MAP=stats($pdl_MAP);
    my @stats_X=stats($pdl_table_X);
    my $x_min=$stats_X[3]->at(0)-0.15;
    my $x_max=$stats_X[4]->at(0)+0.15;
    my $y_min,$y_max;
if (($Y_MIN==0)&&($Y_MAX==0)) {
	$y_min=0;#$stats_MAP[3]->at(0);
	$y_max=$stats_MAP[4]->at(0);
    } else {
	$y_min=$Y_MIN;
	$y_max=$Y_MAX;
    }

    $pdl_color->inplace->setnantobad;
    $pdl_color->inplace->lclip($c_min);
    $pdl_color->inplace->setvaltobad($c_min);
    $pdl_color->inplace->setbadtoval(-1e12);
    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    if (($c_min==0)&&($c_max==0)) {
	$c_min=$stats_color[3]->at(0);
	$c_max=$stats_color[4]->at(0);
    }

    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1.0; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsvp(0.1,0.95,0.15,0.9);
    pgswin($x_min,$x_max,$y_min,$y_max);
    pglabel("$x_label","$y_label","");
    pgsci(1);
    pgsci(1);
    pgsch(2.1);
    my $pdl_c=50+50*($pdl_color-$c_min)/($c_max-$c_min);
    pgsci(1);
    pgsch(1.2);
    pgsls(1);
    pgsci(1);

#
#
#
    pgsitf(0);
    pgslw(2);
    for ($j=0;$j<$NY;$j++) {
	my $color=$pdl_c->at($j+1);
	pgsci($color);
	my $pdl_Y=$pdl_MAP->slice(":,($j)");
	#print "$j | $pdl_Y\n";
	my @a_Y=list($pdl_Y);
	my @a_X=list($pdl_table_X);
	pgline($#a_X+1,\@a_X,\@a_Y);
    }
    pgsci(1);
    pgslw(1);

    pgptxt($x_min+0.05*($x_max-$x_min),$y_max-0.09*($y_max-$y_min),0,0,"$label");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);
    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);
#    pgbox("BC",0,0,"BC",0,0);


#

    pgsci(1);
#    for ($i=60;$i<90;$i=$i+5) {
#	my $val=$c_min+apr(($i-50)/(40/($c_max-$c_min)));
#	my $x=$x_min+($i-50)*($x_max-$x_min)/40;
#	my $y=$y_max-0.05*($y_max-$y_min);
#	pgptxt($x,$y,0,0.5,$val);
#    }

    pgsci(1);
#    pgwedg("TI",0.0,1.5,50,100,"");    
    pgsch(1.8);
    pgwedg("RI",0.0,3,$c_min,$c_max,"$c_label");
    pgsci(1);
#    pgbox("SBC",0,0,"SBC",0,0);



    pgclos();
    pgend();

    # RETURN
}



sub plot_diag_MAP_color_scale {
    my $pdl_MAP=$_[0];
    my $pdl_color=$_[1];
    my $pdl_ctab=$_[2];
    my $dev=$_[3];    
    my $x_min=$_[4];
    my $x_max=$_[5];
    my $y_min=$_[6];
    my $y_max=$_[7];
    my $c_min=$_[8];
    my $c_max=$_[9];
    my $x_label=$_[10];
    my $y_label=$_[11];
    my $c_label=$_[12];
    my $label=$_[13];
    my $do_cont=$_[14];
    my $do_imag=$_[15];
    my $h_MAP=$pdl_MAP->gethdr();
    my $crval1=$$h_MAP{CRVAL1};
    my $crpix1=$$h_MAP{CRPIX1};
    my $cdelt1=$$h_MAP{CDELT1};
    my $crval2=$$h_MAP{CRVAL2};
    my $crpix2=$$h_MAP{CRPIX2};
    my $cdelt2=$$h_MAP{CDELT2};
    my @tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
    my ($NX,$NY)=$pdl_MAP->dims();
    $pdl_color->inplace->setnantobad;
    $pdl_color->inplace->lclip($c_min);
    $pdl_color->inplace->setvaltobad($c_min);
    $pdl_color->inplace->setbadtoval(-1e12);
    my $i,$j,$k;
    my @stats_color=stats($pdl_color);
    if (($c_min==0)&&($c_max==0)) {
	$c_min=$stats_color[3]->at(0);
	$c_max=$stats_color[4]->at(0);
    }
    pgbegin(0,$dev,1,1);
    my $nc=254;
    my @r=list($pdl_ctab->slice(":,(0)"));
    my @g=list($pdl_ctab->slice(":,(1)"));
    my @b=list($pdl_ctab->slice(":,(2)"));
    my @l=list($pdl_ctab->slice(":,(3)"));
    my $bright=1; 
    my $contrast=0.5;
    pgscir(50,100);
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgscf(2.0);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsvp(0.1,0.95,0.15,0.9);
    pgswin($x_min,$x_max,$y_min,$y_max);
    pglabel("$x_label","$y_label","");
    pgsci(1);
    pgsci(1);
    pgsch(2.1);
    my $pdl_c=50+50*($pdl_color-$c_min)/($c_max-$c_min);

    pgsci(1);
    pgsch(1.2);
    pgsls(1);
    pgsci(1);

#
#
#
    if ($do_imag==1) {
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @map=list($pdl_color);
	pgimag(\@map,$NX,$NY,1,$NX,1,$NY,$c_min,$c_max,\@tr);
	pgsci(1);
	pgslw(1);

    }

    
#    print "DO = $do_imag,$do_cont\n";
    if ($do_cont==1) {
	my @stats_MAP=stats($pdl_MAP);
	pgslw(4);
	pgsci(1);
	pgsls(1);
	my @levels;
	my $nlevels=5;
	my $ii;
	my @map=list($pdl_MAP);
	my $min=$stats_MAP[3]->at(0);
	my $max=$stats_MAP[4]->at(0);
	for ($ii=0;$ii<$nlevels;$ii++) {
	    $levels[$ii]=0.05*($max-$min)+($max-$min)/5*$ii;
	}

	pgcont(\@map,$NX,$NY,1,$NX,1,$NY,\@levels,$nlevels,\@tr);
	pgsci(1);
	pgslw(1);

    }


    pgsls(1);
    pgptxt($x_min+0.05*($x_max-$x_min),$y_max-0.09*($y_max-$y_min),0,0,"$label");
    pgsch(1.2);
    pgsls(1);
    pgsci(1);


    my @cut_x,@cut_y,@cut_y2,@cut_y3,@cut_y4;
    for ($i=0;$i<200;$i++) {
	$cut_x[$i]=-2+0.02*$i;
	$cut_y[$i]=-0.7+0.2-3.67*$cut_x[$i];
	$cut_y2[$i]=-1.7+0.5-3.67*$cut_x[$i];
	
	$cut_y3[$i]=0.61/($cut_x[$i]-0.05)+1.3;
	$cut_y4[$i]=0.61/($cut_x[$i]-0.47)+1.19;
	if ($cut_y3[$i]>1.1) {
	    $cut_y3[$i]=-10;
	}
	if ($cut_y4[$i]>1.1) {
	    $cut_y4[$i]=-10;
	}
	
    }
    pgsls(1);
    pgslw(3);
    pgsci(3);
    pgline($#cut_x,\@cut_x,\@cut_y3);
    pgsls(2);
    pgslw(3);
    pgsci(1);
    pgline($#cut_x,\@cut_x,\@cut_y4);

    pgslw(1);
    pgsls(1);
    pgbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);




#
# PGWEB
#

    pgsci(1);
#    for ($i=60;$i<90;$i=$i+5) {
#	my $val=$c_min+apr(($i-50)/(40/($c_max-$c_min)));
#	my $x=$x_min+($i-50)*($x_max-$x_min)/40;
#	my $y=$y_max-0.05*($y_max-$y_min);
#	pgptxt($x,$y,0,0.5,$val);
#    }

    pgsci(1);
#    pgwedg("TI",0.0,1.5,50,100,"");    
    pgsch(1.8);
    pgwedg("RI",0.0,3,$c_min,$c_max,"$c_label");
    pgsci(1);
#    pgbox("SBC",0,0,"SBC",0,0);



    pgclos();
    pgend();

    # RETURN
}





sub CALIFA_cmap {
    my $nc=254;
    my @r,@l,@g,@b;
    while($cmap=<DATA>) {
	chop($cmap);
	@data=split(" ",$cmap);
	$nc=$data[0];
	$r[$nc-1]=$data[1]/255;
	$g[$nc-1]=$data[2]/255;
	$b[$nc-1]=$data[3]/255;
	$l[$nc]=$nc/255;
    }
    my $bright=1.0; 
    my $contrast=0.5;
    $r[0]=1.0;
    $g[0]=1.0;
    $b[0]=1.0;
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);
    return 0;
}

sub radial_sum_ring {
    my $pdl_cube=$_[0];
    my $pdl_r=$_[1];
    my $Dr=$_[2];
    my $r_max=$_[3];
    my ($nx,$ny,$nz)=$pdl_cube->dims();   
    my @stats_r=stats($pdl_r);
    if ($r_max==0) {
	$r_max=$stats_r[4]->at(0);
    }
    my $nr=int($r_max/$Dr)+1;
    my $pdl_out=zeroes($nz,$nr);
    my $pdl_r_sum=zeroes($nz,$nr);
    my $i,$j;
    for ($i=0;$i<$nx;$i++) {
	for ($j=0;$j<$ny;$j++) {
	    my $r_now=$pdl_r->at($i,$j);
	    my $i_r=int($r_now/$Dr);
	    if ($i_r<$nr) {
		my $t=$pdl_out->slice(":,($i_r)");
		my $pdl_slice=$pdl_cube->slice("($i),($j),:");
		my $pdl_norm=$pdl_slice/$pdl_slice;
		$pdl_norm->inplace->setnantobad;
		$pdl_norm->inplace->setbadtoval(0);
		$t .= $t+$pdl_slice*$pdl_norm;
		my $t=$pdl_r_sum->slice(":,($i_r)");
		$t .= $t+$pdl_norm;
	    }
	}
    }
#    my $pdl_xchg=$pdl_out->xchg(0,1);
 #   $pdl_xchg=$pdl_xchg/$pdl_r_sum;
  #  $pdl_out=$pdl_xchg->xchg(0,1);
    $pdl_out=$pdl_out/$pdl_r_sum;
    return $pdl_out;
}

sub radial_sum_ring_20150526 {
    my $pdl_cube=$_[0];
    my $pdl_r=$_[1];
    my $Dr=$_[2];
    my $r_max=$_[3];
    my ($nx,$ny,$nz)=$pdl_cube->dims();   
    my @stats_r=stats($pdl_r);
    if ($r_max==0) {
	$r_max=$stats_r[4]->at(0);
    }
    my $nr=int($r_max/$Dr)+1;
    my $pdl_out=zeroes($nz,$nr);
    my $pdl_r_sum=zeroes($nr);
    my $i,$j;
    for ($i=0;$i<$nx;$i++) {
	for ($j=0;$j<$ny;$j++) {
	    my $r_now=$pdl_r->at($i,$j);
	    my $i_r=int($r_now/$Dr);
	    if ($i_r<$nr) {
		my $t=$pdl_out->slice(":,($i_r)");
		my $pdl_slice=$pdl_cube->slice("($i),($j),:");
		$t .= $t+$pdl_slice;
		my @stats_slice=stats($pdl_slice);
		if (($stats_slice[0]->at(0))>0) {
		    my $val=$pdl_r_sum->at($i_r);
		    $val=$val+1;
		    set($pdl_r_sum,$i_r,$val);
		}	    
	    }
	}
    }
    my $pdl_xchg=$pdl_out->xchg(0,1);
    $pdl_xchg=$pdl_xchg/$pdl_r_sum;
    $pdl_out=$pdl_xchg->xchg(0,1);
    return $pdl_out;
}




__DATA__
      1.00000      0.00000      5.64000      130.080
      2.00000      0.00000      7.64000      133.080
      3.00000      0.00000      11.4600      135.620
      4.00000      0.00000      15.2800      138.160
      5.00000      0.00000      19.1000      140.700
      6.00000      0.00000      22.9200      143.240
      7.00000      0.00000      26.7400      145.780
      8.00000      0.00000      30.5600      148.320
      9.00000      0.00000      34.3800      150.860
      10.0000      0.00000      38.2000      153.400
      11.0000      0.00000      42.0200      155.940
      12.0000      0.00000      45.8400      158.480
      13.0000      0.00000      49.6600      161.020
      14.0000      0.00000      53.4800      163.560
      15.0000      0.00000      57.3000      166.100
      16.0000      0.00000      61.1200      168.640
      17.0000      0.00000      64.9400      171.180
      18.0000      0.00000      68.7600      173.720
      19.0000      0.00000      72.5800      176.260
      20.0000      0.00000      76.4000      178.800
      21.0000      0.00000      80.2200      181.340
      22.0000      0.00000      84.0400      183.880
      23.0000      0.00000      87.8600      186.420
      24.0000      0.00000      91.6800      188.960
      25.0000      0.00000      95.5000      191.500
      26.0000      0.00000      99.3200      194.040
      27.0000      0.00000      103.140      196.580
      28.0000      0.00000      106.960      199.120
      29.0000      0.00000      110.780      201.660
      30.0000      0.00000      114.600      204.200
      31.0000      0.00000      118.420      206.740
      32.0000      0.00000      122.240      209.280
      33.0000      0.00000      126.060      211.820
      34.0000      0.00000      129.880      214.360
      35.0000      0.00000      133.700      216.900
      36.0000      0.00000      137.520      219.440
      37.0000      0.00000      141.340      221.980
      38.0000      0.00000      145.160      224.520
      39.0000      0.00000      148.980      227.060
      40.0000      0.00000      152.800      229.600
      41.0000      0.00000      156.620      232.140
      42.0000      0.00000      160.440      234.680
      43.0000      0.00000      164.260      237.220
      44.0000      0.00000      168.080      239.760
      45.0000      0.00000      171.900      242.300
      46.0000      0.00000      175.720      244.840
      47.0000      0.00000      179.540      247.380
      48.0000      0.00000      183.360      249.920
      49.0000      0.00000      187.180      252.460
      50.0000      0.00000      191.000      255.000
      51.0000      5.10000      187.180      249.900
      52.0000      10.2000      183.360      244.800
      53.0000      15.3000      179.540      239.700
      54.0000      20.4000      175.720      234.600
      55.0000      25.5000      171.900      229.500
      56.0000      30.6000      168.080      224.400
      57.0000      35.7000      164.260      219.300
      58.0000      40.8000      160.440      214.200
      59.0000      45.9000      156.620      209.100
      60.0000      51.0000      152.800      204.000
      61.0000      56.1000      148.980      198.900
      62.0000      61.2000      145.160      193.800
      63.0000      66.3000      141.340      188.700
      64.0000      71.4000      137.520      183.600
      65.0000      76.5000      133.700      178.500
      66.0000      81.6000      129.880      173.400
      67.0000      86.7000      126.060      168.300
      68.0000      91.8000      122.240      163.200
      69.0000      96.9000      118.420      158.100
      70.0000      102.000      114.600      153.000
      71.0000      107.100      110.780      147.900
      72.0000      112.200      106.960      142.800
      73.0000      117.300      103.140      137.700
      74.0000      122.400      99.3200      132.600
      75.0000      127.500      95.5000      127.500
      76.0000      132.600      91.6800      122.400
      77.0000      137.700      87.8600      117.300
      78.0000      142.800      84.0400      112.200
      79.0000      147.900      80.2200      107.100
      80.0000      153.000      76.4000      102.000
      81.0000      158.100      72.5800      96.9000
      82.0000      163.200      68.7600      91.8000
      83.0000      168.300      64.9400      86.7000
      84.0000      173.400      61.1200      81.6000
      85.0000      178.500      57.3000      76.5000
      86.0000      183.600      53.4800      71.4000
      87.0000      188.700      49.6600      66.3000
      88.0000      193.800      45.8400      61.2000
      89.0000      198.900      42.0200      56.1000
      90.0000      204.000      38.2000      51.0000
      91.0000      209.100      34.3800      45.9000
      92.0000      214.200      30.5600      40.8000
      93.0000      219.300      26.7400      35.7000
      94.0000      224.400      22.9200      30.6000
      95.0000      229.500      19.1000      25.5000
      96.0000      234.600      15.2800      20.4000
      97.0000      239.700      11.4600      15.3000
      98.0000      244.800      7.64001      10.2000
      99.0000      249.900      3.82000      5.10000
      100.000      255.000      0.00000      0.00000
      101.000      255.000      3.30000      0.00000
      102.000      255.000      6.60000      0.00000
      103.000      255.000      9.90000      0.00000
      104.000      255.000      13.2000      0.00000
      105.000      255.000      16.5000      0.00000
      106.000      255.000      19.8000      0.00000
      107.000      255.000      23.1000      0.00000
      108.000      255.000      26.4000      0.00000
      109.000      255.000      29.7000      0.00000
      110.000      255.000      33.0000      0.00000
      111.000      255.000      36.3000      0.00000
      112.000      255.000      39.6000      0.00000
      113.000      255.000      42.9000      0.00000
      114.000      255.000      46.2000      0.00000
      115.000      255.000      49.5000      0.00000
      116.000      255.000      52.8000      0.00000
      117.000      255.000      56.1000      0.00000
      118.000      255.000      59.4000      0.00000
      119.000      255.000      62.7000      0.00000
      120.000      255.000      66.0000      0.00000
      121.000      255.000      69.3000      0.00000
      122.000      255.000      72.6000      0.00000
      123.000      255.000      75.9000      0.00000
      124.000      255.000      79.2000      0.00000
      125.000      255.000      82.5000      0.00000
      126.000      255.000      85.8000      0.00000
      127.000      255.000      89.1000      0.00000
      128.000      255.000      92.4000      0.00000
      129.000      255.000      95.7000      0.00000
      130.000      255.000      99.0000      0.00000
      131.000      255.000      102.300      0.00000
      132.000      255.000      105.600      0.00000
      133.000      255.000      108.900      0.00000
      134.000      255.000      112.200      0.00000
      135.000      255.000      115.500      0.00000
      136.000      255.000      118.800      0.00000
      137.000      255.000      122.100      0.00000
      138.000      255.000      125.400      0.00000
      139.000      255.000      128.700      0.00000
      140.000      255.000      132.000      0.00000
      141.000      255.000      135.300      0.00000
      142.000      255.000      138.600      0.00000
      143.000      255.000      141.900      0.00000
      144.000      255.000      145.200      0.00000
      145.000      255.000      148.500      0.00000
      146.000      255.000      151.800      0.00000
      147.000      255.000      155.100      0.00000
      148.000      255.000      158.400      0.00000
      149.000      255.000      161.700      0.00000
      150.000      255.000      165.000      0.00000
      151.000      251.000      162.800      1.10000
      152.000      247.000      160.600      2.20000
      153.000      243.000      158.400      3.30000
      154.000      239.000      156.200      4.40000
      155.000      235.000      154.000      5.50000
      156.000      231.000      151.800      6.60000
      157.000      227.000      149.600      7.70000
      158.000      223.000      147.400      8.80000
      159.000      219.000      145.200      9.90000
      160.000      215.000      143.000      11.0000
      161.000      211.000      140.800      12.1000
      162.000      207.000      138.600      13.2000
      163.000      203.000      136.400      14.3000
      164.000      199.000      134.200      15.4000
      165.000      195.000      132.000      16.5000
      166.000      191.000      129.800      17.6000
      167.000      187.000      127.600      18.7000
      168.000      183.000      125.400      19.8000
      169.000      179.000      123.200      20.9000
      170.000      175.000      121.000      22.0000
      171.000      171.000      118.800      23.1000
      172.000      167.000      116.600      24.2000
      173.000      163.000      114.400      25.3000
      174.000      159.000      112.200      26.4000
      175.000      155.000      110.000      27.5000
      176.000      151.000      107.800      28.6000
      177.000      147.000      105.600      29.7000
      178.000      143.000      103.400      30.8000
      179.000      139.000      101.200      31.9000
      180.000      135.000      99.0000      33.0000
      181.000      131.000      96.8000      34.1000
      182.000      127.000      94.6000      35.2000
      183.000      123.000      92.4000      36.3000
      184.000      119.000      90.2000      37.4000
      185.000      115.000      88.0000      38.5000
      186.000      111.000      85.8000      39.6000
      187.000      107.000      83.6000      40.7000
      188.000      103.000      81.4000      41.8000
      189.000      99.0000      79.2000      42.9000
      190.000      95.0000      77.0000      44.0000
      191.000      91.0000      74.8000      45.1000
      192.000      87.0000      72.6000      46.2000
      193.000      83.0000      70.4000      47.3000
      194.000      79.0000      68.2000      48.4000
      195.000      75.0000      66.0000      49.5000
      196.000      71.0000      63.8000      50.6000
      197.000      67.0000      61.6000      51.7000
      198.000      63.0000      59.4000      52.8000
      199.000      59.0000      57.2000      53.9000
      200.000      55.0000      55.0000      55.0000
      201.000      58.0182      56.9091      58.0182
      202.000      61.0364      58.8182      61.0364
      203.000      64.0545      60.7273      64.0545
      204.000      67.0727      62.6364      67.0727
      205.000      70.0909      64.5455      70.0909
      206.000      73.1091      66.4545      73.1091
      207.000      76.1273      68.3636      76.1273
      208.000      79.1455      70.2727      79.1455
      209.000      82.1636      72.1818      82.1636
      210.000      85.1818      74.0909      85.1818
      211.000      88.2000      76.0000      88.2000
      212.000      91.2182      77.9091      91.2182
      213.000      94.2364      79.8182      94.2364
      214.000      97.2545      81.7273      97.2545
      215.000      100.273      83.6364      100.273
      216.000      103.291      85.5455      103.291
      217.000      106.309      87.4546      106.309
      218.000      109.327      89.3636      109.327
      219.000      112.345      91.2727      112.345
      220.000      115.364      93.1818      115.364
      221.000      118.382      95.0909      118.382
      222.000      121.400      97.0000      121.400
      223.000      124.418      98.9091      124.418
      224.000      127.436      100.818      127.436
      225.000      130.455      102.727      130.455
      226.000      133.473      104.636      133.473
      227.000      136.491      106.545      136.491
      228.000      139.509      108.455      139.509
      229.000      142.527      110.364      142.527
      230.000      145.545      112.273      145.545
      231.000      148.564      114.182      148.564
      232.000      151.582      116.091      151.582
      233.000      154.600      118.000      154.600
      234.000      157.618      119.909      157.618
      235.000      160.636      121.818      160.636
      236.000      163.655      123.727      163.655
      237.000      166.673      125.636      166.673
      238.000      169.691      127.545      169.691
      239.000      172.709      129.455      172.709
      240.000      175.727      131.364      175.727
      241.000      178.745      133.273      178.745
      242.000      181.764      135.182      181.764
      243.000      184.782      137.091      184.782
      244.000      187.800      139.000      187.800
      245.000      190.818      140.909      190.818
      246.000      193.836      142.818      193.836
      247.000      196.855      144.727      196.855
      248.000      199.873      146.636      199.873
      249.000      202.891      148.545      202.891
      250.000      205.909      150.455      205.909
      251.000      208.927      152.364      208.927
      252.000      211.945      154.273      211.945
      253.000      214.964      156.182      214.964
      254.000      217.982      158.091      217.982
      255.000      221.000      160.000      221.000
      256.000      255.000      165.000      0.00000
