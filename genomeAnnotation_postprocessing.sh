
## To generate final annotations

	Step 1: Merge annotations from both mikado runs (blast vs noblast) only incorporate unique noblast annotations.

		sed -i '/\tsuperlocus\t/d' blastmikado.gff
		sed -i '/\tsuperlocus\t/d' noblastmikado.gff
		sed -i 's/mikado./bmikado./g' blastmikado.gff
		sed -i 's/geneid=/yoid:/g' blastmikado.gff
		sed -i 's/geneid=/yoid:/g' noblastmikado.gff
		sed -i 's/Name=/name:/g' noblastmikado.gff
		sed -i 's/Name=/name:/g' blastmikado.gff

		gffcompare -R -r blastmikado.gff -p mikmerged -o mikmerged noblastmikado.gff

		awk '{ if ($4 == "u") print $0}' mikmerged.tracking > noblastunique.txt

		sed  -i 's/|/\:/g' noblastunique.txt
		awk -F":" '{ print $3 }' noblastunique.txt > ID_noblastunique.txt

		cat ID_noblastunique.txt | xargs -I % grep % noblastmikado.gff >> mRNA-noblastmikado_unique.gff

		awk -F":" '{ print $2 }' noblastunique.txt > geneID_noblastunique.txt
		cat geneID_noblastunique.txt | sort | uniq > geneID_noblast_uniq.txt
		sed -i 's/$/;/' geneID_noblast_uniq.txt
		awk '{ if ($3 ~ "gene") print $0 }' noblastmikado.gff > noblastmikado_geneonly.gff
		cat geneID_noblast_uniq.txt | xargs -I % grep % noblastmikado_geneonly.gff >> gene-noblastmikado_unique.gff

		cat gene-noblastmikado_unique.gff mRNA-noblastmikado_unique.gff > noblastmikado_uniqueall.gff
		./gff3sort.pl --precise --chr_order natural noblastmikado_uniqueall.gff > noblastmikado_uniqueall_sorted.gff
		cat blastmikado.gff noblastmikado_uniqueall_sorted.gff > mikado_all.gff3
		./gff3sort.pl --precise --chr_order natural mikado_all.gff3 > mikado_all_sorted.gff3


	Step 2: Merge yo and fb annotations (yovsfb) only incorporate unique fb annotations.

		sed -i 's/geneID=/yoid:/g' yo.gff
		sed -i 's/Name=/name:/g' fb.gff

		gffcompare -R -r yo.gff -p mikmerged -o mikmerged fb.gff

		awk '{ if ($4 == "u") print $0}' mikmerged.tracking > fbunique.txt

		sed  -i 's/|/\:/g' fbunique.txt
		awk -F":" '{ print $3 }' fbunique.txt > ID_fbunique.txt

		cat ID_fbunique.txt | xargs -I % grep % fb.gff >> mRNA-fb_unique.gff

		awk -F":" '{ print $2 }' fbunique.txt > geneID_fbunique.txt
		cat geneID_fbunique.txt | sort | uniq > geneID_fbunique_uniq.txt
		sed -i 's/$/;/' geneID_fbunique_uniq.txt
		awk '{ if ($3 ~ "gene") print $0 }' fb.gff > fb_geneonly.gff
		cat geneID_fbunique_uniq.txt | xargs -I % grep % fb_geneonly.gff >> gene-fb_unique.gff

		cat gene-fb_unique.gff mRNA-fb_unique.gff > fb_uniqueall.gff
		./gff3sort.pl --precise --chr_order natural fb_uniqueall.gff > fb_uniqueall_sorted.gff
		cat yo.gff fb_uniqueall_sorted.gff > fbyo_all.gff3
		./gff3sort.pl --precise --chr_order natural fbyo_all.gff3 > fbyo_all_sorted.gff3

	Step 3: Merge final mikado and fbyo annotations (mikado vs fbyo) only incorporate fbyo annotations.

		gffcompare -R -r mikado_all_sorted.gff3 -p mikmerged -o mikmerged fbyo_all_sorted.gff3

		awk '{ if ($4 == "u") print $0}' mikmerged.tracking > fbyounique.txt

		sed  -i 's/|/\:/g' fbyounique.txt
		awk -F":" '{ print $3 }' fbyounique.txt > ID_fbyounique.txt

		cat ID_fbyounique.txt | xargs -I % grep % fbyo_all_sorted.gff3 >> mRNA-fbyo_unique.gff

		cat mikado_all_sorted.gff3 mRNA-fbyo_unique.gff > dper.gff3
		./gff3sort.pl --precise --chr_order natural dper.gff3 > dper_sorted.gff3

	Step 4: Incorporate ncRNAs from the Machado lab.

		###sed -i 's/\tncRNA\t/\tnoncoding_transcript\t/' dper_ncRNA_machado.gff ** this did not work

		sed -i 's/\tnoncoding_transcript\t/\tncRNA\t/' dper_ncRNA_machado.gff
		sed  's/\tncRNA_gene\t/\tncRNA\t/' dper_sorted.gff > dper_sorted_mod.gff
		* agat_sp_merge_annotations.pl --gff dper_sorted_mod.gff --gff dper_ncRNA_machado.gff --out dper_final_2.gff
		* ./gff3sort.pl --precise --chr_order natural dper_final_2.gff > dper_final_sorted_2.gff

		###cat dper_sorted.gff3 dper_ncRNA_machado.gff > dper_final.gff3
		###

	Step 5: fix cds phases

		agat_sp_fix_cds_phases.pl -g dper_final_sorted_2.gff -f m40.fasta -o dper_final_fixedphases_2.gff3
		./gff3sort.pl --precise --chr_order natural dper_final_fixedphases_23.gff > dper_final_fixedphases_sorted_2.gff3
		agat_sp_manage_IDs.pl -gff dper_final_fixedphases_sorted_2.gff3 -o dperCM_2.gff
		./gff3sort.pl --precise --chr_order natural dperCM_2.gff > dperCM_2_sorted.gff

	Step 6: get statistics and prepare for ortholog clustering

		agat_sp_keep_longest_isoform.pl -gff dperCM_2_sorted.gff -o dperCM_2_longestisoforms.gff
		agat_sp_statistics.pl --gff dperCM_2_sorted.gff --gs m40.fasta -d -o dper_2_annotationStats
		agat_sp_extract_sequences.pl -g dperCM_2_longestisoforms.gff -f m40.fasta -p --clean_final_stop -o dperCM_2_proteins.fasta
		agat_sp_filter_incomplete_gene_coding_models.pl -gff dperCM_2_longestisoforms.gff --fasta m40.fasta -o dperCM_2_completecds.gff
		agat_sp_extract_sequences.pl -g dperCM_2_completecds.gff -f m40.fasta -p --clean_final_stop -o dperCM_2_completeproteins.fasta

		Step 7: get transcript sequences

		agat_sp_extract_sequences.pl -g dperCM_2_sorted.gff -f m40.fasta -t mrna --clean_final_stop -o dperCM_2_mrna.fasta
		agat_sp_extract_sequences.pl -g dperCM_2_sorted.gff -f m40.fasta -t ncrna --clean_final_stop -o dperCM_2_ncrna.fasta
		agat_sp_extract_sequences.pl -g dperCM_2_sorted.gff -f m40.fasta -t transcript --clean_final_stop -o dperCM_2_transcript.fasta

		agat_sp_extract_sequences.pl -g dperCM_2_longestisoforms.gff -f m40.fasta -t cds --clean_final_stop -o dperCM_2_cdslongestisoforms.fasta



		** ** agat_sp_fix_features_locations_duplicated.pl  -f dperCM_2_longestisoforms.gff -o dperCM_2_longestisoforms_nodup.gff ** **
