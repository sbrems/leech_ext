import os
import fnmatch
import pickle

import numpy as np
import matplotlib.pyplot as mpl
from scipy.interpolate import griddata

from astropy.io import fits
#import jFits

from leech import bpm
from leech import dewarp
from leech import fix_pix

def nod_subtract(side, directory='../../data/sat/', 
                      pickle_dir='../../processed_data/sat/',
                      outdirectory='../../processed_data/sat/bcponlm/'):
    #check if folders exist and create
    for director in [directory,pickle_dir,outdirectory,outdirectory+'SX/',outdirectory+'DX/']:
        if not os.path.exists(director):
            os.makedirs(director)
    #check folders are empty and delete old files
    for director in [outdirectory+'SX/',outdirectory+'DX/']:
        for the_file in os.listdir(director):
            file_path = os.path.join(director, the_file)
            if os.path.isfile(file_path):
                os.remove(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path) #del also subfolders
    X_DIM=1024
    Y_DIM=1024

    #/////////////////////////////////////////////////////////
    #get the filenames and nod positions
    #read saved star positions, sort files into nods, save that info
    #/////////////////////////////////////////////////////////

    filenames=sorted(os.listdir(directory))#this sort is important!
    n_files=len(filenames)
    if side == 0: 
        fo=open(pickle_dir+'/SX_sat_star_positions.pkl','rb')
        dat=pickle.load(fo)
        fo.close()
        filenumbers = np.array(dat['filenumbers'])
        star_cen_x =  np.array(dat['star_cen_x'])
        star_cen_y =  np.array(dat['star_cen_y'])
    if side == 1: 
        fo=open(pickle_dir+'/DX_sat_star_positions.pkl','rb')
        dat=pickle.load(fo)
        fo.close()
        filenumbers = np.array(dat['filenumbers'])
        star_cen_x =  np.array(dat['star_cen_x'])
        star_cen_y =  np.array(dat['star_cen_y'])

    nod_beams=[0]
    n_nods=1
    nod_beginnings=[0]
    nod_ends=[]

    for i in xrange(1,n_files):
        no_nod=0
        if (star_cen_y[i] > 600 and star_cen_y[i-1] <= 600) or (star_cen_y[i] <= 600 and star_cen_y[i-1] > 600):
            n_nods+=1
            nod_beginnings.append(i)
            nod_ends.append(i-1)
            nod_beams.append(1-nod_beams[i-1])
        else: 
            nod_beams.append(nod_beams[i-1])

    nod_ends.append(n_files-1)
    nods=[]

    #save nod-beams
    if side == 0: 
        fo=open(pickle_dir+'/SX_nod_beams.pkl','wb')
        pickle.dump({'nod_beams':nod_beams, 'nod_beginnings':nod_beginnings, 'nod_ends':nod_ends},fo)
    if side == 1: 
        fo=open(pickle_dir+'/DX_nod_beams.pkl','wb')
        pickle.dump({'nod_beams':nod_beams, 'nod_beginnings':nod_beginnings, 'nod_ends':nod_ends},fo)


    mpl.plot(range(n_files), star_cen_y)
    mpl.plot(nod_beginnings, star_cen_y[nod_beginnings],marker='*')
    mpl.draw()

    #/////////////////////////////////////////////////////////
    #make nod images (median is time consuming...)
    #/////////////////////////////////////////////////////////
    print 'Make Nod Images'
    for h in xrange(n_nods):
        print '\r Nod %i of %i'%(h,n_nods-1)
        nod_scratch=[]
        for k in xrange(nod_beginnings[h], nod_ends[h]+1):
            nod_scratch.append(fits.open(directory+filenames[k])[0].data[0,:,:])
        #if nod_ends[h]-nod_beginnings[h] == 0: 
        #    nods.append(np.array(nod_scratch)) #not needed and gives wrong dim
        nods.append(np.median(nod_scratch, axis=0))
    print '\n'
    #/////////////////////////////////////////////////////////
    #do nod subtraction, overscan subtraction, pad, crop, and bin in sequence
    #/////////////////////////////////////////////////////////

    #prepare things that only need to be computed once.
    good_neighbors=fix_pix.find_neighbors(bpm,n=8)
    #distortion (from Anne-Lise email)  
    Kx = [[-2.1478925,    0.0058138110,  -7.6396687e-06,   2.5359596e-09], 
          [1.0109149,  -2.3826537e-05,   2.8458629e-08,  -9.3206482e-12], 
          [-2.1164521e-05,   5.3115381e-08,  -6.6315643e-11,   2.2888432e-14], 
          [1.2983972e-08,  -4.1253977e-11,   5.1637044e-14,  -1.5988376e-17]]
    Ky = [[9.2717864,      0.98776733,   4.3514612e-06,   9.3450739e-09], 
          [-0.013617797,  -3.9526096e-05,   8.1204222e-08,  -5.2048768e-11], 
          [1.1313247e-05,   6.7127301e-08,  -1.6531988e-10,   1.0656544e-13], 
          [1.6283111e-09,  -2.7723216e-11,   8.2118035e-14,  -5.3695050e-17]]
    #transpose is necessary for IDL vs. Python definition of polynomial coefficients
    dewarp_coords=dewarp.make_dewarp_coordinates(bpm.shape,
                                                 np.array(Kx).T,np.array(Ky).T)
    print 'Subtract Nods'
    for h in xrange(n_files):
        print 'Subtracting image %i of %i '%(h,n_files-1)
        #get the image
        print 'get image'
        image=(fits.open(directory+filenames[h])[0].data[0,:,:]).astype(float)
        #change if cubes are used

        #determine which nod to subtract
        print 'determine which nod to subtract'
        this_nod=-1
        for k in xrange(n_nods):
            if (h >= nod_beginnings[k]) and (h <= nod_ends[k]): 
                this_nod=k
        if this_nod != -1:
            if this_nod % 2 == 1: 
                nod_to_use=this_nod-1
            elif this_nod % 2 == 0: 
                nod_to_use=this_nod+1
            elif (n_nods % 2 == 1) and (this_nod == n_nods-1): 
                nod_to_use=this_nod-1 #for odd numbers of nods,last nod

            #subtract the nod

#            if h== 290:
#                print 'Debugging...'
#                ipdb.set_trace()
            print 'subtract the nod'
            image-=nods[nod_to_use]

            #column subtract
            #not exactly an overscan correction...
            for i in xrange(X_DIM):
                column=image[:,i]
                stdev_column=np.std(column)
                median_column=np.median(column)
                stdev_column=np.std(column[(np.abs(column-median_column) < (3*stdev_column))])
                col_bias=np.median(column[(np.abs(column-median_column) < (3*stdev_column))])
                image[:,i]-=col_bias

            print 'bfixpix'
            fix_pix.correct_with_precomputed_neighbors(image,good_neighbors)#inplace

            print 'dewarp'
            image=dewarp.dewarp_with_precomputed_coords(image,dewarp_coords)

            #pad
            large_image=np.zeros((1024+512,1024+512))
            large_image[256:1280,256:1280]=image

            #crop
            if (star_cen_x[h]+256-300 > 0) and (star_cen_x[h]+256+299 < 1535) and (star_cen_y[h]+256-300 > 0) and (star_cen_y[h]+256+299 < 1535):
                cropped_image=large_image[star_cen_y[h]+256-300:star_cen_y[h]+256+299, star_cen_x[h]+256-300:star_cen_x[h]+256+299]

            #bin
            bin_image=np.zeros((300,300))
            for i in xrange(300):
                for j in xrange(300):
                    pixel=cropped_image[0+i*2:1+i*2, 0+j*2:1+j*2]
                    bin_image[i,j]=np.median(pixel)

            #save
            save_dir=outdirectory
            if side == 0: 
                save_filename=save_dir+'SX/bcpon_SX_'+filenames[h]
            if side == 1: 
                save_filename=save_dir+'DX/bcpon_DX_'+filenames[h]

            fits.writeto(save_filename, bin_image)
