import numpy as np
from leech import dewarp
import jFits
import datetime

image=np.zeros((1024,1024))
image[::30,:]=1
image[:,::30]=1

Kx = [[-2.1478925,    0.0058138110,  -7.6396687e-06,   2.5359596e-09],
      [1.0109149,  -2.3826537e-05,   2.8458629e-08,  -9.3206482e-12],
      [-2.1164521e-05,   5.3115381e-08,  -6.6315643e-11,   2.2888432e-14],
      [1.2983972e-08,  -4.1253977e-11,   5.1637044e-14,  -1.5988376e-17]]
#Kx = [[-2.1478925,    0.0058138110,  -7.6396687e-06,   2.5359596e-09],
#      [1.0109149,  4e-3,   2.8458629e-08,  -9.3206482e-12],
#      [-2.1164521e-05,   5.3115381e-08,  -6.6315643e-11,   2.2888432e-14],
#      [1.2983972e-08,  -4.1253977e-11,   5.1637044e-14,  -1.5988376e-17]]

Ky = [[9.2717864,      0.98776733,   4.3514612e-06,   9.3450739e-09],
      [-0.013617797,  -3.9526096e-05,   8.1204222e-08,  -5.2048768e-11],
      [1.1313247e-05,   6.7127301e-08,  -1.6531988e-10,   1.0656544e-13],
      [1.6283111e-09,  -2.7723216e-11,   8.2118035e-14,  -5.3695050e-17]]
print 'start dewarp...'
#dim=dewarp(image,np.array(Kx).T,np.array(Ky).T,order=0)
tt=datetime.datetime.now()
dim=dewarp(image,np.array(Kx).T,np.array(Ky).T,order=3)
print datetime.datetime.now()-tt
#dim2=dewarp(image,np.array(Kx).T,np.array(Ky).T,order=2)
#dim1=dewarp(image,np.array(Kx).T,np.array(Ky).T,order=1)
print 'end dewarp...'

disp=jFits.jInteractive_Display(dim,vmin=0,vmax=1)
disp.a.figure.savefig('python_dewarped.pdf',dpi=1000)


jFits.pyfits.writeto('grid.fits',image,clobber=True)
jFits.pyfits.writeto('python_dewarped_grid.fits',dim,clobber=True)
