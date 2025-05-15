#input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/genome_metaT/genome.tsv
#input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/DE_analysis/g-diversity-core-metrics114K/rarefied_table/genome.tsv
#input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaT/genome.tsv
input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/species_pfam_metaT/species_pfam.tsv
dbdir=/mnt/zarrinpar/Pynchon/Databases/wol2

# collapse to get pfam-TPM

fundir=$dbdir/function
md=$fundir/pfam

woltka tools collapse -i $input -m $md/orf-to-pfam.map.xz -n $md/pfam_description.txt -o ${input/%.tsv/pfam.tsv}
