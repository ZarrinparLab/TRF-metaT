#to run Songbird you need to be in environment songbird-qiime2-2020.6

path=/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/species_func
feat_file=$path/species_pfam_clean_noNT_dark.qza
md_file=/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/metaT_metadata_noNT_dark.txt
output_path=$path/songbird_species_byCond_dark_ftfa

#path=/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/genus_func/genusNfunc
#feat_file=$path/genus_pfam_clean_rmdoubletons_noNT_dark.qza
#md_file=/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/metaT_metadata_noNT_dark.txt
#output_path=$path/songbird_genuspfam_dark_byCond_ftfa

qiime songbird multinomial \
  --i-table $feat_file \
  --m-metadata-file $md_file \
  --p-formula "C(condition, Treatment('FA'))" \
  --p-epochs 5000 \
  --p-differential-prior 1 \
  --p-summary-interval 1 \
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
  --o-differentials $output_path/null_differentials.qza \
  --o-regression-stats $output_path/null-stats.qza \
  --o-regression-biplot $output_path/null-biplot.qza

qiime songbird summarize-paired \
	--i-regression-stats $output_path/regression-stats.qza \
	--i-baseline-stats $output_path/null-stats.qza \
	--o-visualization $output_path/paired-summary.qzv
