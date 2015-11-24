import gc
import os
import fnmatch
import pickle
import numpy as np
from astropy.io import fits
from leech import rot
import datetime

def annular_PCA(nod_beam, x_cen, y_cen, n_PCA, side, 
                            directory='processed_data/cxsxbcponlm/SX_0/',
                            outdir='processed_data/PCA/SX_0/'):

    side={0:'S',1:'D'}[side]

    #////////////////////////////////////////////////////////
    #Make polar and cartisian coordinate systems
    #////////////////////////////////////////////////////////
    X_DIM=300
    Y_DIM=300
    y,x=np.indices((Y_DIM,X_DIM))
    #polar
    radius=np.abs((x-x_cen)+1j*(y-y_cen))
    radius_1d=radius.flatten()

    #////////////////////////////////////////////////////////
    #Get the data
    #////////////////////////////////////////////////////////
    filenames=sorted(fnmatch.filter(os.listdir(directory),'*.fits'))

    n_images=len(filenames)
    #data=fltarr(X_DIM,Y_DIM,n_images)

    #////////////////////////////////////////////////////////
    #start a loop to do region-by-region
    #////////////////////////////////////////////////////////
    n_r=140

    r_subt_in=np.arange(n_r)*1.+3.
    r_subt_out=np.arange(n_r)*1.+3.+1.

    r_opt_in=r_subt_in-7.
    r_opt_out=r_subt_in+9.

    r_opt_in[r_opt_in < 6]=5

    #////////////////////////////////////////////////////////
    #Do the whole thing in loops for memory allocation
    #////////////////////////////////////////////////////////

    combined_median=np.zeros((Y_DIM,X_DIM))
    combined_data=np.zeros((X_DIM,Y_DIM,n_images))

    before=datetime.datetime.now()
    for l in xrange(n_r):
        print l, n_r
        data=[]
        for fn in filenames: 
            hdul=fits.open(directory+fn)
            data.append(hdul[0].data.copy())
            del hdul
            gc.collect()
        data=np.array(data)
        data-=np.median(data,axis=0)[None,:,:]

        #////////////////////////////////////////////////////////
        #construct a masked region
        #////////////////////////////////////////////////////////

        #opt_mask=np.zeros(X_DIM*Y_DIM)
        opt_mask=np.logical_and(radius_1d>r_opt_in[l],radius_1d<=r_opt_out[l])#bool
        #opt_mask[where(radius_1d lt r_opt_in[l])]=1
        #opt_mask[where(radius_1d ge r_opt_out[l])]=1

        n_pix_opt_masked=np.sum(opt_mask.astype(int))

        #sub_mask=fltarr(X_DIM*Y_DIM)
        sub_mask=np.logical_and(radius_1d>r_subt_in[l],radius_1d<=r_subt_out[l])#bool
        #sub_mask[where(radius_1d lt r_subt_in[l])]=1
        #sub_mask[where(radius_1d ge r_subt_out[l])]=1

        n_pix_sub_masked=np.sum(sub_mask.astype(int))

        #////////////////////////////////////////////////////////
        #make matrix of non-opt_masked pixels for PCA
        #////////////////////////////////////////////////////////

        A=np.zeros((opt_mask.sum(),n_images))
        for h in xrange(n_images):
            data_1d=data[h,:,:].flatten()
            A[:,h]=data_1d[opt_mask]

        #////////////////////////////////////////////////////////
        #Find the principal components (eigenvectors)
        #////////////////////////////////////////////////////////

        mean_image=np.mean(A,axis=0)
        A-=mean_image[None,:] #subtract off the mean image(the median is already gone JMS)

        #linear algebra described on the wiki page for eigenfaces
        #https://en.wikipedia.org/wiki/Eigenface#Computing_the_eigenvectors
        AT_A=np.dot(A.T,A)
        eigenvalues, eigenvectors = np.linalg.eig(AT_A)
        idx=np.argsort(eigenvalues)[::-1]#python doesn't sort like idl...
        eigenvectors=eigenvectors[:,idx]
        eigenvectors=np.dot(A,eigenvectors)
        #normalize the eigenvectors
        for h in xrange(n_images): 
            eigenvectors[:,h]/=(np.sum(eigenvectors[:,h]**2.))**0.5#vector in rows for idl, columns in python...

        #////////////////////////////////////////////////////////
        #Fit the individual images and reform to 2d
        #////////////////////////////////////////////////////////
        #loop through image by image
        for h in xrange(n_images):
            #fit the image
            #project each image onto each (requested) eigenvector
            PCA_image_1d=np.zeros(n_pix_opt_masked)
            for k in xrange(n_PCA):
                coefficient=np.dot(A[:,h],eigenvectors[:,k])   
                PCA_image_1d+=(eigenvectors[:,k]*coefficient)

            #subtract from original
            PCA_subtracted_image=A[:,h]-PCA_image_1d

            #insert the opt_masked pixels back in
            PCA_subtracted_with_opt_mask_1d=np.zeros(X_DIM*Y_DIM,dtype=float)
            PCA_subtracted_with_opt_mask_1d[opt_mask]=PCA_subtracted_image

            #only save the subtraction region (not the opt region)
            #PCA_subtracted_image*=sub_mask

            #reform to 2-d
            #overwriting data... 
            data[h,:,:]=PCA_subtracted_with_opt_mask_1d.reshape(Y_DIM,X_DIM)


        #////////////////////////////////////////////////////////
        #rotate and combine
        #////////////////////////////////////////////////////////

        fo=open(directory+side+'X_'+str(nod_beam)+'_PA.pkl','rb')
        combined_images_PA=pickle.load(fo)

        for h in xrange(n_images):
            # wrote in python, works for integer pivot coords JMS
            # rotating in the correct direction? JMS
            data[h,:,:]=rot(data[h,:,:],combined_images_PA[h], (x_cen, y_cen))

        #median combine and save
        sub_mask_2d=sub_mask.reshape(Y_DIM,X_DIM) 
        combined_median[sub_mask_2d]=np.median(data,axis=0)[sub_mask_2d] #only save the subtraction region, not the optimization region
        print datetime.datetime.now()-before
        before=datetime.datetime.now()

    save_filename=outdir+'/median_annular_'+str(n_PCA)+'.fits'
    fits.writeto(save_filename,combined_median)
