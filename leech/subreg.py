import numpy as np
from scipy.signal import fftconvolve

from gaussfitter import gaussfit

from leech.fix_pix import bfixpix
from leech import find_max_star

def subreg(reference, images):
    bfixpix(reference,np.isnan(reference),n=8)
    kernel=reference[::-1,::-1]
    shifts=[]
    for im in images:
        bfixpix(im,np.isnan(im),n=8)
        cor = fftconvolve(im,kernel,mode='same')
        y,x=find_max_star(cor)
        g=gaussfit(cor[max(0, y-40):min(y+40, cor.shape[0]),
                       max(0, x-40):min(x+40, cor.shape[1])])
        shiftx=np.rint(cor.shape[1]/2.) - max(0,x-40)-g[2]
        shifty=np.rint(cor.shape[0]/2.) - max(0,y-40)-g[3]
        shifts.append((shifty,shiftx))
    return shifts


                       

