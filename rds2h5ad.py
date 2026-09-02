import scanpy as sc
import pandas as pd
import scipy.io
import json
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--indir', help='input dir', required=True)
parser.add_argument('-o', '--outdir', help='output dir', required=True)
parser.add_argument('-s', '--sample_name', help='sample name', required=True)
parser.add_argument('-im', '--image_dir', help='image dir', required=True)
parser.add_argument('-r', '--res', help='hires or lowres', default='lowres', required=True)
args = parser.parse_args()

# 1. 读取矩阵及基因/细胞信息
X = scipy.io.mmread(args.indir + "/counts.mtx").T.tocsr()
genes = pd.read_csv((args.indir + "/features.tsv"), header=None)[0].values
barcodes = pd.read_csv((args.indir + "/barcodes.tsv"), header=None)[0].values
###注意由于前面没有设置reindex，所以这里的细胞名称为全数字构成，因此很容易被识别为数字类型，所以在后续处理中全部要进行字符串化。
barcode = barcodes.astype(str)

# 2. 构建 AnnData 对象
adata = sc.AnnData(X=X)
adata.var_names = genes
adata.obs_names = barcode

# 3. 导入 Metadata
metadata = pd.read_csv((args.indir + "/metadata.csv"), index_col=0)
###这里metadata的index也为全数字，因此也要特别处理。
metadata.index = metadata.index.astype(str)
adata.obs = metadata.loc[adata.obs_names] # 确保顺序一致

# 4. 导入空间坐标 (关键步：Scanpy 默认将空间坐标存在 obsm['spatial'] 中)
spatial_coords = pd.read_csv((args.indir + "/spatial_coordinates_" + args.res + ".csv"), index_col=0)
###spatial_coords的index同上。
spatial_coords.index = spatial_coords.index.astype(str)
# Seurat 的坐标列名通常是 imagecol 和 imagerow，或者 x 和 y
# 确保顺序与 adata.obs 一致
spatial_coords = spatial_coords.loc[adata.obs_names]
adata.obsm['spatial'] = spatial_coords[['imagecol', 'imagerow']].values

# 5. 加入底图
img = mpimg.imread(args.image_dir + "/tissue_" + args.res + "_image.png")
###非RGB图像这里要手动伪造三通道图像
img_rgb = np.stack([img, img, img], axis=-1).astype(np.float32)

# 构建 Scanpy 要求的 uns 空间字典结构
# 'spatial' -> 'library_id' -> 'images' & 'scalefactors'

with open((args.image_dir + '/scalefactors_json.json'), 'r', encoding='utf-8') as f:
    data = json.load(f)

library_id = "slice1" 
scale_factor_from_r = data['tissue_' + args.res + '_scalef']  # 请替换为你在R里查到的真实 tissue_lowres_scalef

if (args.res == 'lowres'):
    adata.uns['spatial'] = {
        library_id: {
            'images': {
                'lowres': img_rgb
            },
            'scalefactors': {
                'tissue_lowres_scalef': 1, 
                'spot_diameter_fullres': data['spot_diameter_fullres']  # 替换为实际分辨率spot直径
            }
        }
    }
else:
    adata.uns['spatial'] = {
        library_id: {
            'images': {
                'hires': img_rgb
            },
            'scalefactors': {
                'tissue_hires_scalef': 1,
                'spot_diameter_fullres': data['spot_diameter_fullres'] 
            }
        }
    }

plt.rcParams["figure.figsize"] = [12, 12]
sc.pl.spatial(adata, color='nCount_Spatial', library_id='slice1', img_key=args.res, spot_size=1, alpha=0.8)
plt.savefig(args.outdir + '/spatial_expre.png')

# 6. 保存为 h5ad
adata.write_h5ad(args.outdir + '/' + args.sample_name + "_spatial_data_scanpy.h5ad")
