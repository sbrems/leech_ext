import os
import time
import gc
import pickle

import numpy as np
from astropy.io import fits
def combine_frames(side, nod_number,
                   directory='../../processed_data/sat/xsxbcponlm/SX_0/',
                   raw_dir='../../data/sat/',
                   save_dir='../../processed_data/sat/cxsxbcponlm/SX_0/',
                   date=None,n_combine=10):
    #make directories
    for director in [directory,raw_dir,save_dir]:
        if not os.path.exists(director):
            os.makedirs(director)
    print 'combine_frames'
    X_DIM=300
    Y_DIM=300

    #get date

    #////////////////////////////////////////////////////////
    #Get filenames and rotation info
    #////////////////////////////////////////////////////////
    filenames=sorted(os.listdir(directory))

    if (nod_number == 0) and (side == 0):
        save_PA_filename=save_dir+'/SX_0_PA.pkl'
    elif (nod_number == 1) and (side == 0):
        save_PA_filename=save_dir+'/SX_1_PA.pkl'
    elif (nod_number == 0) and (side == 1):
        save_PA_filename=save_dir+'/DX_0_PA.pkl'
    elif (nod_number == 1) and (side == 1):
        save_PA_filename=save_dir+'/DX_1_PA.pkl'

    filenumbers=[f[22:27] for f in filenames]
#    print filenumbers
    n_images=len(filenumbers)

    #open raw data to get PAs
    par_angle=[]
    for h in xrange(n_images):
        print h, filenumbers[h]
        par_filename=raw_dir+'lm_'+date+'_'+filenumbers[h]+'.fits'
        header_image=fits.open(par_filename)[0]
        par_angle.append(header_image.header['LBT_PARA'])
    print par_angle
    print 'done_PA'
    #time.sleep(5)#necessary?? JMS

    #////////////////////////////////////////////////////////
    #Get the files, n_combine at a time
    #////////////////////////////////////////////////////////

    combined_n_images=n_images/n_combine
    combined_images_PA=[]

    for h in xrange(combined_n_images):
        individual_images=[]
        individual_images_PA=[]

        which_images=np.arange(n_combine)+(h*n_combine)
        print 'which'
        print which_images
        print len(which_images)
        for k in xrange(n_combine):
            print filenumbers[which_images[k]]
            if side == 0: 
                filename=directory+'/xsxbcpon_SX_lm_'+date+'_'+filenumbers[which_images[k]]+'.fits'
            if side == 1: 
                filename=directory+'/xsxbcpon_DX_lm_'+date+'_'+filenumbers[which_images[k]]+'.fits'
            hdul=fits.open(filename)
            individual_images.append(hdul[0].data.copy())
            individual_images_PA.append(par_angle[which_images[k]])
            del hdul
            gc.collect()
        print len(individual_images_PA)
        print len(individual_images)
          
        median_PA=np.median(individual_images_PA)
        good_PA=np.where(np.abs(individual_images_PA-median_PA) < 2)[0]#I added the abs JMS
        combined_image=np.median(np.take(individual_images,good_PA, axis=0), axis=0)

        save_filename=save_dir+('%05i'%h)+'.fits'
        fits.writeto(save_filename, combined_image)

        combined_images_PA.append(np.median(np.take(individual_images_PA,good_PA)))
        print h, len(good_PA), combined_images_PA[-1], np.min(individual_images_PA[0]), np.max(individual_images_PA)

    combined_images_PA=-1*np.array(combined_images_PA)
    fo=open(save_PA_filename,'wb')
    pickle.dump(combined_images_PA,fo)
    fo.close()
