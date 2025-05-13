#to run RPCA you need to be in environment qiime2-2021.4

path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/
mtb_file=$path/m16s/07_taxon_filtered_asv_table_dada2.qza
md_file=$path/metadata.TRF_combined_wLD.tab
output_path=$path/diversity_analysis/rpca_results_m16s/
#comparison=condition
comparison=cond_fasted
#comparison=condphase

qiime deicode rpca \
  --i-table $mtb_file \
  --p-min-feature-count 15 \
  --p-min-sample-count 500 \
  --o-biplot $output_path/ordination.qza \
  --o-distance-matrix $output_path/distance_matrix.qza

qiime emperor biplot \
  --i-biplot $output_path/ordination.qza \
  --m-sample-metadata-file $md_file \
  --p-ignore-missing-samples \
  --o-visualization $output_path/biplot.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison}-significance.qzv
