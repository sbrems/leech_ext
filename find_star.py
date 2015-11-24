import gc
import fnmatch
import pickle
import matplotlib.pyplot as mpl
from leech import *
from leech import bpm as bpm0 
from leech import background as background0
import jFits

def find_star(side, directory='../../data/sat/',outdirectory='../../processed_data/sat/'):
    #/////////////////////////////////////////////////////////
    #get the filenames
    #/////////////////////////////////////////////////////////
    #this sort is important for all the codes to work!
    filenames=sorted(os.listdir(directory))
    n_files=len(filenames)

    #/////////////////////////////////////////////////////////
    #get the x and y centers
    #/////////////////////////////////////////////////////////

    filenumbers=[]
    star_cen_x=[]
    star_cen_y=[]

    for h in xrange(n_files):
        filenumbers.append(filenames[h][10:16])#check JMS
        filename=directory+filenames[h]
        hdul=fits.open(filename)
        image=hdul[0].data[0,:,:]-background0
        del hdul[0].data
        
        if side == 0: 
            #overscan, but this is the whole chip... not left
            image=image[4:1020,4:1020]
            bpm=bpm0[4:1020,4:1020]
        if side == 1: 
            image=image[4:1020,512:1020] # right
            bpm=bpm0[4:1020,512:1020]

        #find bright pixels
        threshold=0.90*np.max(image)
        test=np.logical_and(image > threshold, bpm == 0)
        y_high_counts, x_high_counts = test.nonzero()
        if side == 0:  
            star_cen_x.append(np.median(x_high_counts) + 4)
        if side == 1:  
            star_cen_x.append(np.median(x_high_counts) + 512)
        star_cen_y.append(np.median(y_high_counts) + 4)

        print filenumbers[-1], star_cen_x[-1], star_cen_y[-1]
        gc.collect()#astropy.io.fits is inefficient. This is a fix.

    #if side eq 0 then save, filenumbers, star_cen_x, star_cen_y, filename='../../processed_data/sat/SX_sat_star_positions.sav'
    if side == 0:
        fo=open(outdirectory+'SX_sat_star_positions.pkl','wb')
        pickle.dump({'filenumbers':filenumbers, 'star_cen_x':star_cen_x, 'star_cen_y':star_cen_y},fo)
        fo.close()
    #if side eq 1 then save, filenumbers, star_cen_x, star_cen_y, filename='../../processed_data/sat/DX_sat_star_positions.sav'
    if side == 1:
        fo=open(outdirectory+'DX_sat_star_positions.pkl','wb')
        pickle.dump({'filenumbers':filenumbers, 'star_cen_x':star_cen_x, 'star_cen_y':star_cen_y},fo)
        fo.close()
    mpl.plot(range(n_files), star_cen_y)
