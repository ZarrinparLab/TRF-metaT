input=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/genome_metaG/genome.tsv
dbdir=/mnt/zarrinpar/Pynchon/Databases/wol2

# collapse to get pfam-TPM

fundir=$dbdir/function
md=$fundir/pfam

woltka tools collapse -i $input -m $md/orf-to-pfam.map.xz -n $md/pfam_description.txt -o ${input/%.tsv/pfam.tsv}
