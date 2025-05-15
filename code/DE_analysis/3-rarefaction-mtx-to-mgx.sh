#environment qiime2-2021.4

#get depth mgx
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/
mtb_file=$path/genome_metaG/genome_clean_noNT.txt
md_file=$path/metaG_metadata_noNT.txt

biom convert \
 -i $mtb_file \
 -o ${mtb_file/%.txt/.biom} \
 --table-type="OTU table" \
 --to-hdf5

qiime tools import \
 --input-path ${mtb_file/%.txt/.biom} \
 --output-path ${mtb_file/%.txt/.qza} \
 --type FeatureTable[Frequency]
 
 qiime feature-table summarize \
  --i-table ${mtb_file/%.txt/.qza} \
  --o-visualization ${mtb_file/%.txt/.qzv} \
  --m-sample-metadata-file $md_file
#mgx genome depth is 113,935

#get depth mtx
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/
mtb_file=$path/genome_metaT/genome_clean_noNT.tsv
md_file=$path/metaT_metadata_ztcat_noNT.txt

biom convert \
 -i $mtb_file \
 -o ${mtb_file/%.tsv/.biom} \
 --table-type="OTU table" \
 --to-hdf5

qiime tools import \
 --input-path ${mtb_file/%.tsv/.biom} \
 --output-path ${mtb_file/%.tsv/.qza} \
 --type FeatureTable[Frequency]
 
 qiime feature-table summarize \
  --i-table ${mtb_file/%.tsv/.qza} \
  --o-visualization ${mtb_file/%.tsv/.qzv} \
  --m-sample-metadata-file $md_file
 #mtx genome depth is 60,610,350
 
#get rarefied table that matches depth of MGX
#output_path=$path/DE_analysis/g-diversity-core-metrics114K/

#qiime diversity core-metrics \
#  --i-table ${mtb_file/%.tsv/.qza} \
#  --p-sampling-depth 113935 \
#  --m-metadata-file $md_file \
#  --output-dir $output_path
