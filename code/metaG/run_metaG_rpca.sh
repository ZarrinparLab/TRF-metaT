#to run RPCA you need to be in environment qiime2-2021.4

path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG
mtb_file=$path/pfam_notnorm/pfam_clean_noNT.txt
md_file=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT_wfasting.txt
output_path=$path/pfam_notnorm/metaG_rpca_results
#comparison=condition
comparison=cond_fasted

#biom convert \
# -i $mtb_file \
# -o ${mtb_file/%.txt/.biom} \
# --table-type="OTU table" \
# --to-hdf5

#qiime tools import \
# --input-path ${mtb_file/%.txt/.biom} \
# --output-path ${mtb_file/%.txt/.qza} \
# --type FeatureTable[Frequency]

#qiime deicode rpca \
#  --i-table ${mtb_file/%.txt/.qza}\
#  --p-min-feature-count 10 \
#  --p-min-sample-count 500 \
#  --o-biplot $output_path/ordination.qza \
#  --o-distance-matrix $output_path/distance_matrix.qza

#qiime emperor biplot \
#  --i-biplot $output_path/ordination.qza \
#  --m-sample-metadata-file $md_file \
#  --p-ignore-missing-samples \
#  --o-visualization $output_path/biplot.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison}-significance.qzv
