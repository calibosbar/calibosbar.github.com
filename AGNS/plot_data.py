# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 16:38:32 2015

@author: jeffrey
"""

import pyfits
import numpy as np
import matplotlib.pyplot as plt

hdulist = pyfits.open('/home/jeffrey/Documentos/laptop/dataproducts/AGNS/UGC11680NED01.p_e.rad_SFH_lum_V.fits.gz')
len(hdulist)
hdulist.info()
hdulist[0].header.keys()
data = hdulist[0].data
plt.imshow(data)
plt.imshow(data[20,:])

