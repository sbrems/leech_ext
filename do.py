from leech.find_star import find_star 
from leech.nod_subtract import nod_subtract
from leech.xreg import xreg
from leech.combine_frames import combine_frames
from leech.annular_PCA import annular_PCA

find_star(0,directory='/disk1/brems/131024/KapAnd/sat/',
            outdirectory='/disk1/brems/131024/KapAnd/processed_data/sat/')

#nod_subtract(0,directory='/home/jstone/Research/LEECH/data/test/sat/',
#               pickle_dir='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/',
#               outdirectory='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/bcponlm/')
#

#xreg(0,0,directory='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/bcponlm/SX/',
#         pickle_dir='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/',
#         out_dir='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/xsxbcponlm/')

#combine_frames(0,0,directory='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/xsxbcponlm/SX_0/',
#                   save_dir='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/cxsxbcponlm/SX_0/',
#                   raw_dir='/home/jstone/Research/LEECH/data/test/sat/',
#                   date='150105')

#annular_PCA(0,150,150,1,0,directory='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/cxsxbcponlm/SX_0/',
#                          outdir='/home/jstone/Research/LEECH/crunched/test/processed_data/sat/PCA/SX_0/')
