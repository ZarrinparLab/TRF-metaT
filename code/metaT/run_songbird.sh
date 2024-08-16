#to run Songbird you need to be in environment songbird-qiime2-2020.6

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/genus_pfam
#feat_file=$path/genus_pfam_clean_rmdoubletons_noNT.qza
#md_file=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt
#output_path=$path/songbird_rmdoubletons_FA

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/genus_pfam
#feat_file=$path/genus_pfam_clean_rmdoubletons_noNT_light.qza
#md_file=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt
#output_path=$path/songbird_rmdl_light

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/genus_pfam
#feat_file=$path/genus_pfam_BSHonly_clean_rmzero_noNT.qza
#md_file=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt
#output_path=$path/songbird_BSH_rmzero_FA

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/genus_pfam
#feat_file=$path/genus_pfam_BSHonly_clean_rmzero_noNT_light.qza
#md_file=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt
#output_path=$path/songbird_BSH_rmzero_light

#path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/species_pfam
#feat_file=$path/species_pfam_BSHonly_clean_rmzero_noNT_light.qza
#md_file=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt
#output_path=$path/songbird_BSH_rmzero_light_FA

path=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/BSH
feat_file=$path/genomeL_noNT_rmdbton.qza
md_file=/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt
output_path=$path/songbird_genomeL_FA

qiime songbird multinomial \
  --i-table $feat_file \
  --m-metadata-file $md_file \
  --p-formula "C(condition, Treatment('FA'))" \
  --p-epochs 5000 \
  --p-differential-prior 1 \
  --p-summary-interval 1 \
  --p-min-feature-count 0 \
  --p-min-sample-count 0 \
  --o-differentials $output_path/differentials.qza \
  --o-regression-stats $output_path/regression-stats.qza \
  --o-regression-biplot $output_path/regression-biplot.qza


qiime songbird multinomial \
  --i-table $feat_file \
  --m-metadata-file $md_file \
  --p-formula "1" \
  --p-epochs 5000 \
  --p-differential-prior 1 \
  --p-summary-interval 1 \
  --p-min-feature-count 0 \
  --p-min-sample-count 0 \
  --o-differentials $output_path/null_differentials.qza \
  --o-regression-stats $output_path/null-stats.qza \
  --o-regression-biplot $output_path/null-biplot.qza

qiime songbird summarize-paired \
	--i-regression-stats $output_path/regression-stats.qza \
	--i-baseline-stats $output_path/null-stats.qza \
	--o-visualization $output_path/paired-summary.qzv
