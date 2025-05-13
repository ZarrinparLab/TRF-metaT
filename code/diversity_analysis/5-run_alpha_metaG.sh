#environment qiime2-2021.4

path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/
mtb_file=$path/pfam_metaG/pfam_clean_noNT.txt
md_file=$path/metaG_metadata_noNT.txt
output_path=$path/diversity_analysis/diversity-core-metrics12k/

biom convert \
 -i $mtb_file \
 -o ${mtb_file/%.txt/.biom} \
 --table-type="OTU table" \
 --to-hdf5

qiime tools import \
 --input-path ${mtb_file/%.txt/.biom} \
 --output-path ${mtb_file/%.txt/.qza} \
 --type FeatureTable[Frequency]

qiime diversity core-metrics \
  --i-table ${mtb_file/%.txt/.qza} \
  --p-sampling-depth 12000 \
  --m-metadata-file $md_file \
  --output-dir $output_path