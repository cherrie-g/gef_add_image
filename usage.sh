input_gef=/path/to/stereoseq.tissue.gef
bin_size=50
output=/outdir/
sample_name='sample1'
input_img=/path/to/stereoseq_regist.tif
res="hires"

###Step1: read gef file
./read_gef_image.py -i ${input_gef} \
                                      -b ${bin_size} \
                                      -o ${output} \
                                      -s ${sample_name} \
                                      -im ${input_img}

###Step2: convert h5ad to rds
./h5ad2rds.R --infile "${output}/${sample_name}/temp.h5ad" \
                        --outfile "${output}/${sample_name}/${sample_name}_pre.rds" \
                        --image_dir "${output}/${sample_name}/spatial/slice1/"

###Step2-1: replace hires image in rds
./convert_to_hires.R -r "${output}/${sample_name}/${sample_name}_pre.rds" \
                                     --im_dir "${output}/${sample_name}/spatial/slice1/" \
                                     -o "${output}/${sample_name}/"

###Step3: convert rds with image back to h5ad file
./rds2h5ad_pre.R -i "${output}/${sample_name}/${sample_name}_pre.rds" \
                                -r ${res} \
                                -o "${output}/${sample_name}/"

./rds2h5ad.py -i "${output}/${sample_name}/" \
                          -o "${output}/${sample_name}/" \
                          -s ${sample_name} \
                          -im "${output}/${sample_name}/spatial/slice1/" \
                          -r ${res}
