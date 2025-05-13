#environment qiime2-2021.4

path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/
mtb_file=$path/pfam_metaT/pfam_clean_noNT.tsv
md_file=$path/metaT_metadata_ztcat_noNT.txt
output_path=$path/diversity_analysis/diversity-core-metrics2.3M/

biom convert \
 -i $mtb_file \
 -o ${mtb_file/%.tsv/.biom} \
 --table-type="OTU table" \
 --to-hdf5

qiime tools import \
 --input-path ${mtb_file/%.tsv/.biom} \
 --output-path ${mtb_file/%.tsv/.qza} \
 --type FeatureTable[Frequency]

qiime diversity core-metrics \
  --i-table ${mtb_file/%.tsv/.qza} \
  --p-sampling-depth 2300000 \
  --m-metadata-file $md_file \
  --output-dir $output_path