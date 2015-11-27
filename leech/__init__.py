import os
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import rotate, map_coordinates
from astropy.io import fits

__all__=['os','np','griddata','rotate','fits','rot','dewarp',
         'bin_median','find_max_star']

path=os.path.dirname(__file__)
bpm=fits.open(path+'/bpm.fits')[0].data
background=fits.open(path+'/first_order_background.fits')[0].data

def rot(img, angle, pivot):
    '''from 
    http://stackoverflow.com/questions/25458442/rotate-a-2d-image-around-specified-origin-in-python
    rotate an input image by an angle about a pivot.
    INPUTS:
    img: 2d array (image)
    angle: float [Degrees]
    pivot: 2-tuple, only integers are supported'''
    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
    imgP = np.pad(img, [padY, padX], 'constant')
    imgR = rotate(imgP, angle, reshape=False)
    return imgR[padY[0] : -padY[1], padX[0] : -padX[1]]

#def dewarp(image,P,Q,order=3):
#    '''2d polynomial dewarping. Output x,y are mapped to input xt, yt
#    according to:
#    xt=F(x,y|P)
#    yt=F(x,y|Q)
#    see the documentation for numpy.polynomial.polynomial.polygrid2d 
#    for information on how to organize polynomial coefficients in P, and Q
#    INPUTS:
#    image: 2D array
#    P:   2D array of polynomial coefficients for generating xt from x
#    Q:   2D array of polynomial coefficients for generating yt from y'''
#    sh=image.shape
#    xt=np.polynomial.polynomial.polygrid2d(range(sh[0]),range(sh[1]),P)
#    yt=np.polynomial.polynomial.polygrid2d(range(sh[0]),range(sh[1]),Q)
#    return map_coordinates(image,[yt,xt],order=3)

def bin_median(arr,smaller_by_factor=1,returnStd=False):
    '''bin an array arr by creating super pixels of size
    smaller_by_factor*smaller_by_factor and taking the median. 
    Can optionally also return the standard deviation within
    the super-pixels.
    INPUTS:
    arr: 2d array
    smaller_by_factor: integer
    returnStd: bool, default False
    RETURNS binned_array OR binned_array, bin_std'''
    sub_arrs0=[]
    for i in xrange(smaller_by_factor):
        for j in xrange(smaller_by_factor):
            sub_arrs0.append(arr[i::smaller_by_factor,j::smaller_by_factor])
    sub_arrs=[s[:sub_arrs0[-1].shape[0]-1,:sub_arrs0[-1].shape[1]-1] for s in sub_arrs0]
    if returnStd:
        #zip truncates each sub-arr at the length of the minimum length subarr.
        #this ensures every bin has the same number of datapoints, but throws
        #away data if the last bin doesn't have a full share of datapoints.
        return np.median(sub_arrs,axis=0),np.std(sub_arrs,axis=0)
    else:
        return np.median(sub_arrs,axis=0)

def find_max_star(image):
    '''Median smooth an image and find the max pixel.
    The median smoothing helps filter hot pixels and 
    cosmic rays. The median is taken by using bin_median
    with a smaller_by_factor=16'''
    image[np.isnan(image)]=np.median(image[~np.isnan(image)])
    binned=bin_median(image,smaller_by_factor=16)
    y,x=np.transpose((binned==binned.max()).nonzero())[0]
    y*=16
    x*=16
    while True:
        x0=max(x-15,0)
        y0=max(y-15,0)
        patch = image[y0:min(y+15,image.shape[0]),x0:min(x+15,image.shape[1])]
        dy,dx=np.transpose(np.where(patch==patch.max()))[0]
        dy-=(y-y0)
        dx-=(x-x0)
        y=min( max(y+dy,0), image.shape[0]-1)
        x=min( max(x+dx,0), image.shape[1]-1)
        if (dx==0) and (dy==0):
            break
    return y,x




