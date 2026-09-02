import stereo as st
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', help='input_gef_file', required=True)
parser.add_argument('-b', '--binsize', help='square bin size', required=True)
parser.add_argument('-o', '--outdir', help='output dir', required=True)
parser.add_argument('-s', '--sample_name', help='sample name', required=True)
parser.add_argument('-im', '--image', help='input regist.tif', required=True)
args = parser.parse_args()

data_path = args.infile
data = st.io.read_gef(file_path=data_path, bin_size=int(args.binsize), gene_name_index=True)
data.tl.cal_qc()
data.tl.raw_checkpoint()

os.makedirs(args.outdir + '/' + args.sample_name)

adata = st.io.stereo_to_anndata(data, 
                                flavor='seurat', 
                                output=(args.outdir + '/' + args.sample_name + '/temp.h5ad'), 
                                sample_id=args.sample_name, 
                                image=args.image,
                                im_library_id='slice1')
