#to run RPCA you need to be in environment qiime2-2021.4
#need to clean the files using the R scripts first

#untargetted
#path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/
#mtb_file=$path/species_pfam_metaT/species_pfam_BSH_RPOB_clean_rmzero_noNT.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/species_pfam_metaT/rpca_results_BSH_RPOB_rmzero
#comparison=condition

#targetted
#path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/
#mtb_file=$path/BSH_proteindb_metaT/genome_noNT_rmdbton.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/BSH_proteindb_metaT/rpca_results_genome
#comparison=condition

biom convert \
 -i $mtb_file \
 -o ${mtb_file/%.tsv/.biom} \
 --table-type="OTU table" \
 --to-hdf5

qiime tools import \
 --input-path ${mtb_file/%.tsv/.biom} \
 --output-path ${mtb_file/%.tsv/.qza} \
 --type FeatureTable[Frequency]

qiime deicode rpca \
  --i-table ${mtb_file/%.tsv/.qza}\
  --p-min-feature-count 0 \
  --p-min-sample-count 0 \
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

qiime qurro loading-plot \
  --i-table ${mtb_file/%.tsv/.qza} \
  --i-ranks $output_path/ordination.qza \
  --m-sample-metadata-file $md_file \
  --o-visualization $output_path/qurro.qzv
