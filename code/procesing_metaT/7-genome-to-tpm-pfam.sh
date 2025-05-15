#input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/genome_metaT/genome.tsv
#input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/DE_analysis/g-diversity-core-metrics114K/rarefied_table/genome.tsv
#input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaT/genome.tsv
input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/species_pfam_metaT/species_pfam.tsv
dbdir=/mnt/zarrinpar/Pynchon/Databases/wol2

#get TPM tables
woltka normalize -i $input --scale 1k --digits 3 -o ${input/%.tsv/-rpk.tsv}
woltka normalize -i ${input/%.tsv/-rpk.tsv} --scale 1M -o ${input/%.tsv/-TPM.tsv}

# collapse to get pfam-TPM

fundir=$dbdir/function
input2=${input/%.tsv/-TPM.tsv}
md=$fundir/pfam

woltka tools collapse -i $input2 -m $md/orf-to-pfam.map.xz -n $md/pfam_description.txt -o ${input/%.tsv/pfam-TPM.tsv}
