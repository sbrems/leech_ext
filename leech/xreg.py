import os
import pickle
import gc
import contextlib
import sys
import cStringIO

import numpy as np
from scipy.ndimage.interpolation import shift
from astropy.io import fits

from subreg import subreg

def xreg(side, nod_number,directory='../../processed_data/sat/bcponlm/SX', pickle_dir=None,
         out_dir='../../processed_data/sat/xsxbcponlm/'):
    #check if folders exist and create
    for director in [directory,pickle_dir,out_dir]:
        if not os.path.exists(director):
            os.makedirs(director)
    #check folders are empty and delete old files
#   for director in [outdirectory+'SX/',outdirectory+'DX/']:
#     for the_file in os.listdir(director):
#         file_path = os.path.join(director, the_file)
#         if os.path.isfile(file_path):
#             os.remove(file_path)
#         elif os.path.isdir(file_path): shutil.rmtree(file_path) #del also subfolders
    X_DIM=300
    Y_DIM=300

    #/////////////////////////////////////////////////////////
    #/////////////////////////////////////////////////////////
    #DO NODS SEPARATELY
    #/////////////////////////////////////////////////////////
    #/////////////////////////////////////////////////////////
    if side == 0:
        fo=open(pickle_dir+'/SX_nod_beams.pkl','rb')
        data=pickle.load(fo)
        nod_beams=data['nod_beams']
        nod_beginnings=data['nod_beams']
        nod_ends=data['nod_beams']
    if side == 1: 
        fo=open(pickle_dir+'/DX_nod_beams.pkl','rb')
        data=pickle.load(fo)
        nod_beams=data['nod_beams']
        nod_beginnings=data['nod_beams']
        nod_ends=data['nod_beams']

    #nod 0 or 1
    where_nod=(np.array(nod_beams) == nod_number).nonzero()[0]

    #/////////////////////////////////////////////////////////
    #get the files
    #/////////////////////////////////////////////////////////

    filenames=sorted(os.listdir(directory))

    #this is dangerous. Relies on filenames being sorted correctly in find_star.py...
    filenames0=[]
    for ii in where_nod:
        filenames0.append(filenames[ii])
    filenames=filenames0
    n_files=len(filenames)
    print 'Getting the %i images for nod %i...'%(n_files,nod_number)
    images=[]
    for h in xrange(n_files): 
        hdul=fits.open(directory+filenames[h])
        image=hdul[0].data.copy()
        images.append(image)
        del hdul
        gc.collect()

    #/////////////////////////////////////////////////////////
    #median combine and first xreg
    #/////////////////////////////////////////////////////////
#    import ipdb;ipdb.set_trace()
    first_median=np.median(images, axis=0)
    print 'first subreg...'
    with nostdout(): #suppress messages of badpix
        shifts=subreg(first_median,images)
    print 'shifting the images'
    for h in xrange(n_files): 
        images[h]=shift(images[h], shifts[h])#make sure shift is right direction...

    #/////////////////////////////////////////////////////////
    #keep the best 70% of images
    #/////////////////////////////////////////////////////////
    cross_reg=[]
    for h in xrange(n_files):
        cross_reg.append(np.sum((images[h]-first_median)**2.))
    sorted_cross_reg=np.argsort(cross_reg)
    selected_cross_reg=sorted_cross_reg[0:int(0.7*n_files)]
    n_selected=len(selected_cross_reg)

    #/////////////////////////////////////////////////////////
    #median combine and second xreg
    #/////////////////////////////////////////////////////////

    images=np.array(images)[selected_cross_reg,:,:]
    second_median=np.median(images,axis=0)

    print 'second subreg'
    shifts=subreg(second_median,images)
    for h in xrange(n_selected):
        images[h,:,:]=shift(images[h,:,:], shifts[h])#make sure shift is right direction...

    #/////////////////////////////////////////////////////////
    #save
    #/////////////////////////////////////////////////////////

    if side == 0:
        odirectory=out_dir+'SX_'+str(nod_number)+'/'
    if side == 1: 
        odirectory=out_dir+'DX_'+str(nod_number)+'/'
    #make directories
    if not os.path.exists(odirectory):
        os.makedirs(odirectory)
    for h in xrange(n_selected):
        filename=odirectory+'xsx'+filenames[selected_cross_reg[h]]
        print h
        print filename
        try:
            fits.writeto(filename,images[h,:,:])
        except:
            print 'problem'
            raw_input()


#for silencing functions efficiently:
@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = cStringIO.StringIO()
    yield
    sys.stdout = save_stdout
