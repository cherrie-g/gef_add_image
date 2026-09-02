## This suite is code for add background tissue image to seurat object converted from Stereo-seq gef file. And further replace the lowres image to hires image in seurat object. Also, the rds file can convert back to h5ad file.


Environment requirement: stereopy 1.6.2; Seurat v4.

> Maybe not compatible with Seurat v5.

Input files include gef files and the corresponding tissue image that has aligned to expression matrix. 
```
input_gef=/path/to/stereoseq.tissue.gef
bin_size=50  ###bin-size
output=/outdir/
sample_name='sample1'   ###sample name
input_img=/path/to/stereoseq_regist.tif
res="hires"
```

First step is read gef using stereopy package. 
```
###Step1: read gef file
./read_gef_image.py -i ${input_gef} \
                    -b ${bin_size} \
                    -o ${output} \
                    -s ${sample_name} \
                    -im ${input_img}
```

Then use the script provided by stereopy teams [h5ad2rds.R](https://github.com/STOmics/Stereopy/blob/main/docs/source/_static/h5ad2rds.R) (with modification) to aquire rds file.
```
###Step2: convert h5ad to rds
./h5ad2rds.R --infile "${output}/${sample_name}/temp.h5ad" \
             --outfile "${output}/${sample_name}/${sample_name}_pre.rds" \
             --image_dir "${output}/${sample_name}/spatial/slice1/"
```

And replace high-resolution image in rds file.
```
###Step2-1: replace hires image in rds
./convert_to_hires.R -r "${output}/${sample_name}/${sample_name}_pre.rds" \
                     --im_dir "${output}/${sample_name}/spatial/slice1/" \
                     -o "${output}/${sample_name}/"
```

Finally, if you need h5ad file for analysis under python environment, the script below can achieve this purpose.
```
###Step3: convert rds with image back to h5ad file
./rds2h5ad_pre.R -i "${output}/${sample_name}/${sample_name}_pre.rds" \
                 -r ${res} \
                 -o "${output}/${sample_name}/"

./rds2h5ad.py -i "${output}/${sample_name}/" \
              -o "${output}/${sample_name}/" \
              -s ${sample_name} \
              -im "${output}/${sample_name}/spatial/slice1/" \
              -r ${res}
```
