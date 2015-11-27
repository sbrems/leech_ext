import numpy as np
from leech import dewarp
import jFits
import datetime
from skimage import transform as tf

image=np.zeros((1024,1024))
image[::30,:]=1
image[:,::30]=1

a=[-2.14789,1.01091,0.00581381,-2.11645e-05,-2.38265e-05,-7.63967e-06,1.29840e-08,5.31154e-08,2.84586e-08,2.53596e-09,0,-4.12540e-11,-6.63156e-11,-9.32065e-12,0,0,0,5.16370e-14,2.28884e-14,0,0,0,0,0,-1.59884e-17,0,0,0] 

b=[9.27179,-0.0136178,0.987767,1.13132e-05,-3.95261e-05,4.35146e-06,1.62831e-09,6.71273e-08,8.12042e-08,9.34507e-09,0,-2.77232e-11,-1.65320e-10,-5.20488e-11,0,0,0,8.21180e-14,1.06565e-13,0,0,0,0,0,-5.36950e-17,0,0,0]

params=np.array([a,b])
#Kx = [[-2.1478925,    0.0058138110,  -7.6396687e-06,   2.5359596e-09],
#      [1.0109149,  -2.3826537e-05,   2.8458629e-08,  -9.3206482e-12],
#      [-2.1164521e-05,   5.3115381e-08,  -6.6315643e-11,   2.2888432e-14],
#      [1.2983972e-08,  -4.1253977e-11,   5.1637044e-14,  -1.5988376e-17]]
#Ky = [[9.2717864,      0.98776733,   4.3514612e-06,   9.3450739e-09],
#      [-0.013617797,  -3.9526096e-05,   8.1204222e-08,  -5.2048768e-11],
#      [1.1313247e-05,   6.7127301e-08,  -1.6531988e-10,   1.0656544e-13],
#      [1.6283111e-09,  -2.7723216e-11,   8.2118035e-14,  -5.3695050e-17]]
mapping=tf.PolynomialTransform(params=params[::-1,:])
tt=datetime.datetime.now()
dim=tf.warp(image, inverse_map=mapping, order=3)
print datetime.datetime.now()-tt

disp=jFits.jInteractive_Display(dim,vmin=0,vmax=1)
disp.a.figure.savefig('python_dewarped.pdf',dpi=1000)


jFits.pyfits.writeto('grid.fits',image,clobber=True)
jFits.pyfits.writeto('python_dewarped2_grid.fits',dim,clobber=True)
