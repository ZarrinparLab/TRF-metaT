#input=/home/sfloresr/scratch/TRF_multiomics/metatranscript/woltka2_m_results/genome.tsv
#input=/home/sfloresr/scratch/TRF_multiomics/metatranscript/woltka2_m_results/BSH/genome.tsv
input=/home/sfloresr/scratch/TRF_multiomics/metatranscript/woltka2_m_results/species_pfam/species_pfam.tsv
#dbdir=/mnt/zarrinpar/Pynchon/Databases/wol2

#get rpk tables
#woltka normalize -i $input --sizes $dbdir/proteins/coords.txt.xz --scale 1k --digits 3 -o ${input/%.tsv/-rpk.tsv}
woltka normalize -i $input --scale 1k --digits 3 -o ${input/%.tsv/-rpk.tsv}

#get TPM tables
woltka normalize -i ${input/%.tsv/-rpk.tsv} --scale 1M -o ${input/%.tsv/-TPM.tsv}

# collapse using multiple functional classification systems
#fundir=$dbdir/function
#input2=${input/%.tsv/-TPM.tsv}

#output=/home/sfloresr/scratch/TRF_multiomics/metatranscript/woltka2_m_results 
#md=$fundir/uniref

#woltka tools collapse \
#        -i $input2 \
#        -m $md/orf-to-uniref.map.xz \
#        -n $md/uniref_name.txt.xz \
#        -o $output/uniref/uniref-TPM.tsv

#md=$fundir/go

#for x in component function process; do
#  woltka tools collapse -i $output/uniref/uniref-TPM.tsv -m $md/uniref/$x.map.xz -n $md/go_name.txt -o $output/go/$x-TPM.tsv
#  woltka tools collapse -i $output/go/$x-TPM.tsv -m $md/generic/$x.map -n $md/go_name.txt -o $output/go/$x.generic-TPM.tsv
#done

#md=$fundir/pfam
#woltka tools collapse -i $input2 -m $md/orf-to-pfam.map.xz -n $md/pfam_description.txt -o $output/pfam/pfam-TPM.tsv
#woltka tools collapse -i $output/pfam/pfam-TPM.tsv -m $md/pfam_type.txt -o $output/pfam/type-TPM.tsv
#woltka tools collapse -i $output/pfam/pfam-TPM.tsv -m $md/pfam-to-clan.map -n $md/clan_description.txt -o $output/pfam/clan-TPM.tsv

#md=$fundir/kegg
#woltka tools collapse -i $input2 -m $md/orf-to-ko.map.xz -n $md/ko_name.txt -o $output/kegg/ko-TPM.tsv
#woltka tools collapse -i $output/kegg/ko-TPM.tsv -m $md/ko-to-ec.map -o $output/kegg/ec-TPM.tsv
#woltka tools collapse -i $output/kegg/ko-TPM.tsv -m $md/ko-to-reaction.map -n $md/reaction_name.txt -o $output/kegg/reaction-TPM.tsv
#woltka tools collapse -i $output/kegg/reaction-TPM.tsv -m $md/reaction-to-rclass.map -n $md/rclass_name.txt -o $output/kegg/rclass-TPM.tsv
#woltka tools collapse -i $output/kegg/reaction-TPM.tsv -m $md/reaction-to-module.map -n $md/module_name.txt -o $output/kegg/module-TPM.tsv
#woltka tools collapse -i $output/kegg/module-TPM.tsv -m $md/module-to-class.map -o $output/kegg/module_class-TPM.tsv
#woltka tools collapse -i $output/kegg/module-TPM.tsv -m $md/module-to-pathway.map -n $md/pathway_name.txt -o $output/kegg/pathway-TPM.tsv
#woltka tools collapse -i $output/kegg/pathway-TPM.tsv -m $md/pathway-to-class.map -o $output/kegg/pathway_class-TPM.tsv

#md=$fundir/metacyc
#woltka tools collapse -i $input2 -m $md/orf-to-protein.map.xz -n $md/protein_name.txt -o $output/metacyc/protein-TPM.tsv
#woltka tools collapse -i $output/metacyc/protein-TPM.tsv -m $md/protein-to-gene.map -n $md/gene_name.txt -o $output/metacyc/gene-TPM.tsv
#woltka tools collapse -i $output/metacyc/protein-TPM.tsv -m $md/protein-to-enzrxn.map -n $md/enzrxn_name.txt -o $output/metacyc/enzrxn-TPM.tsv
#woltka tools collapse -i $output/metacyc/enzrxn-TPM.tsv -m $md/enzrxn-to-reaction.map -n $md/reaction_name.txt -o $output/metacyc/reaction-TPM.tsv
#woltka tools collapse -i $output/metacyc/reaction-TPM.tsv -m $md/reaction-to-pathway.map -n $md/pathway_name.txt -o $output/metacyc/pathway-TPM.tsv
#woltka tools collapse -i $output/metacyc/pathway-TPM.tsv -m $md/pathway-to-super_pathway.map -n $md/pathway_name.txt -o $output/metacyc/super_pathway-TPM.tsv
#woltka tools collapse -i $output/metacyc/super_pathway-TPM.tsv -m $md/pathway_type.txt -n $md/all_class_name.txt -o $output/metacyc/pathway_type-TPM.tsv

