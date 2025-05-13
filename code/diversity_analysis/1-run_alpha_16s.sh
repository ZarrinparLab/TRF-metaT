#environment qiime2-2021.4

path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/
mtb_file=$path/m16s/07_taxon_filtered_asv_table_dada2.qza
md_file=$path/metadata.TRF_combined_wLD.tab
output_path=$path/diversity_analysis/diversity-core-metrics6k_16s/

qiime diversity core-metrics \
  --i-table $mtb_file \
  --p-sampling-depth 6000 \
  --m-metadata-file $md_file \
  --output-dir $output_path