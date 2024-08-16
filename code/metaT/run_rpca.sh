#to run RPCA you need to be in environment qiime2-2021.4

path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
mtb_file=$path/pfam/pfam_clean_noNT_rmFA09b.tsv
md_file=$path/metaT_metadata_ztcat_noNT_wfasting.txt
output_path=$path/pfam/rpca_results_rmFA09b
comparison=condition

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/pfam/pfam_clean_noNT.tsv
#md_file=$path/metaT_metadata_ztcat_noNT_wfasting.txt
#output_path=$path/pfam/rpca_results
#comparison=condition
#comparison=cond_fasted

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/genus_pfam/rpca_results_BSH_rmzero
#comparison=condition

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/genus_pfam/genus_pfam_BSH_RPOB_clean_rmzero_noNT.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/genus_pfam/rpca_results_BSH_RPOB_rmzero
#comparison=condition

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/genus_pfam/genus_pfam_clean_rmdoubletons_noNT_dark.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/genus_pfam/rpca_results_dark
#comparison=condition

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/species_pfam/species_pfam_BSH_RPOB_clean_rmzero_noNT.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/species_pfam/rpca_results_BSH_RPOB_rmzero
#comparison=condition

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/species_pfam/species_pfam_BSHonly_clean_rmzero_noNT_dark.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/species_pfam/rpca_results_BSH_dark
#comparison=condition

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/BSH/genomeL_noNT_rmdbton.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/BSH/rpca_results_genomeL
#comparison=condition

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results
#mtb_file=$path/BSH/genome_noNT_rmdbton.tsv
#md_file=$path/metaT_metadata_ztcat_noNT.txt
#output_path=$path/BSH/rpca_results_genome
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
