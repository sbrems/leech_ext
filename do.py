from leech.find_star import find_star 
from leech.nod_subtract import nod_subtract
from leech.xreg import xreg
from leech.combine_frames import combine_frames
from leech.annular_PCA import annular_PCA
import ipdb
#from embed_shell import ipsh

# find_star(0,directory='/disk1/brems/131024/KapAnd/sat/',
#             outdirectory='/disk1/brems/131024/KapAnd/processed_data/sat/')

# nod_subtract(0,directory='/disk1/brems/131024/KapAnd/sat/',
#                pickle_dir='/disk1/brems/131024/KapAnd/processed_data/sat/',
#                outdirectory='/disk1/brems/131024/KapAnd/processed_data/sat/bcponlm/')
#
#
#xreg(0,0,directory='/disk1/brems/131024/KapAnd/processed_data/sat/bcponlm/SX/',
#         pickle_dir='/disk1/brems/131024/KapAnd/processed_data/sat/',
#         out_dir='/disk1/brems/131024/KapAnd/processed_data/sat/xsxbcponlm/')
#
#combine_frames(0,0,directory='/disk1/brems/131024/KapAnd/processed_data/sat/xsxbcponlm/SX_0/',
#                   save_dir='/disk1/brems/131024/KapAnd/processed_data/sat/cxsxbcponlm/SX_0/',
#                   raw_dir='/disk1/brems/131024/KapAnd/sat/',
#                   date='131024',n_combine=20)

annular_PCA(0,150,150,1,0,directory='/disk1/brems/131024/KapAnd/processed_data/sat/cxsxbcponlm/SX_0/',
                          outdir='/disk1/brems/131024/KapAnd/processed_data/sat/PCA/SX_0/')
