##### FINAL ANALYSIS #####

### DEL in Dper supported by INS in Dpse. ##################################################################################################################################

	# Format svim output to a simplified bed file
	awk '{ len = $3 - $2; if ($5 > 1) print $1"\t"$2"\t"$3"\t"$4" "$5" "len }' candidates_deletions.bed > alldeletions.bed

	# Sort the resulting bedfile
	bedtools sort -i alldeletions.bed > alldeletions_sorted.bed

	# Copy the sinteny blocks
	cp /media/DATA5/Javier/Dpse_Dper_annotation_ALL/svmu/all_newsvmu/5_test/pseref2per_collinearblocks_sorted.bed .

	# Copy the perl script to add IDs
	cp /media/DATA5/Javier/Dpse_Dper_annotation_ALL/svim/new_ver/notmasked/INDELpsevsper/DELpseasreference/addIDs_tocollinear.pl .

	# Add IDs to synteny blocks
	perl addIDs_tocollinear.pl pseref2per_collinearblocks_sorted.bed > pseref2per_collinearblocks_sorted_withID.bed

	# Find indel coordinates on collinear blocks between Dper and Dpse allowing partial overlap to include indels located at the break ends of the synteny blocks
	bedtools intersect -wo -a alldeletions_sorted.bed -b pseref2per_collinearblocks_sorted_withID.bed > alldeletions_collinearblocks.bed

	# Filter variants matching collinearblocks by the svim score (>10 is suggested to be considered a true variant)
	awk '{ if ($5 > 10) print $0}' alldeletions_collinearblocks.bed > passed_deletions_collinearblocks.bed

	# Sort filtered variants
	bedtools sort -i passed_deletions_collinearblocks.bed > passed_deletions_collinearblocks_sorted.bed

	# Transform synteny blocks to get Dper coordinates
	awk '{ if ($5 > $6) print $4"\t"$6"\t"$5"\t"$1" "$2" "$3" "$7 ; else print $4"\t"$5"\t"$6"\t"$1" "$2" "$3" "$7 }'  pseref2per_collinearblocks_sorted_withID.bed >  dPERpseref2per_collinearblocks_sorted_withID.bed

	# Sort dper blocks
	bedtools sort -i dPERpseref2per_collinearblocks_sorted_withID.bed > dPERpseref2per_collinearblocks_sorted_withID_sorted.bed

	# Copy insertions in Dpse
	cp /media/DATA5/Javier/Dpse_Dper_annotation_ALL/svim/new_ver/notmasked/pse2per/output/candidates/candidates_novel_insertions.bed .

	# Format the insertions file to a simplified bed file
	awk '{ len = $3 - $2; if ($5 > 1) print $1"\t"$2"\t"$3"\t"$4" "$5" "len }' candidates_novel_insertions.bed > allinsertions.bed

	# Sort the insertions bed file
	bedtools sort -i allinsertions.bed > allinsertions_sorted.bed

	# Find indel coordinates on synteny blocks but now on Dper allowing partial overlap to include indels located at the break ends of the synteny blocks
	bedtools intersect -wo -a allinsertions_sorted.bed -b dPERpseref2per_collinearblocks_sorted_withID_sorted.bed > allinsertions_collinearblocks.bed

	# Sort the insertions sorted bed file
	bedtools sort -i allinsertions_collinearblocks.bed > allinsertions_collinearblocks_sorted.bed

	# Copy the get_block_INDEL.pl script
	cp /media/DATA5/Javier/Dpse_Dper_annotation_ALL/svim/new_ver/notmasked/INDELpsevsper/INSpseasreference/get_block_INDEL.pl .

	# Get INDELs for each synteny block and create a file for each block (for both INS and DEL).
	perl get_block_INDEL.pl passed_deletions_collinearblocks_sorted.bed allinsertions_collinearblocks_sorted.bed  pseref2per_collinearblocks_sorted_withID.bed

	#create a new directory and move indel files
	mkdir block_INDELS
	mv block-*_indels.txt block_INDELS

	# Copy the transform_coordinates_v1.0.pl script
	cp ../../INSpseasreference/block_INDELS/transform_coordinates_v1.0.pl .

	# Transform coordinates to compare dels to ins and to find matching variants in the two genomes
	for i in $(ls *_indels.txt); do perl transform_coordinates_v1.0.pl $i;done

	# Move all blocks with no Dper variants to blocks-nodper/
	ls -l per*.bed | awk '{ if ($5 == "0") print $9 }' | xargs -n 1 -I % mv % blocks-nodper/

	# Sort remaining bed files
	for i in $(ls per_block*.bed); do bedtools sort -i $i > $i.sorted; done
	for i in $(ls pse_block*.bed); do bedtools sort -i $i > $i.sorted; done

	# create a list containing the blocks with dper variants and replace the 'per'
	ls per*.sorted > blocks.txt
	sed -i 's/per_//' blocks.txt

	# Find collinear variants using bedtools closest
	while read i; do bedtools closest -d -k 2 -a pse_$i -b per_$i > $i.closest ; done < blocks.txt

	# create directory for collinear variants and a list of the closest files
	mkdir matching_blocks_0dist
	ls *.sorted.closest > closests_blocks.txt

	# Copy the find_matching_indels.pl script
	cp ../../INSpseasreference/block_INDELS/find_matching_indels.pl .

	# find matching indels using the find_matching_indels.pl script
	perl find_matching_indels.pl closests_blocks.txt

	# create a directory and move the blocks with no matching variants
	mkdir matching_blocks_far
	mv *closest matching_blocks_far

	# Copy the get_matching_indels.pl script
	cp ../../../INSpseasreference/block_INDELS/matching_blocks_0dist/get_matching_indels.pl .

	# For the matching variants in matching_blocks_0dist, find the exact match using the script get_matching_indels.pl
	for i in $(ls *.closest); do perl get_matching_indels.pl $i > $i.exact;done

	# Concatenate all exact files
	cat *.exact > perdels_exact_match_perdels.txt

	# Copy the get_final_exact.pl script
	cp ../../../INSpseasreference/block_INDELS/matching_blocks_0dist/get_final_exact.pl .

	# Get matching variants that are similiar not only on the position but also in the size allowing 25% error, i.e. del: 70 bp then ins: ~70 bp
	perl get_final_exact.pl perdels_exact_match_perdels.txt 0.25

	# Get duplicated per dels from the variants passing the 0.25 cutoff
	awk '{ print $1$2$3 }' pseind_passing_0.25_cutoff.txt | uniq -d > pseind_passing_0.25_repeated.txt
	cat pseind_passing_0.25_repeated.txt | xargs -n 1 -I % awk '{ rep = $1$2$3; if (rep == "%") print $0}' pseind_passing_0.25_cutoff.txt > pseind_passing_0.25_repeated_full.txt

	#manually curate the file pseind_passing_0.25_repeated_full.txt

	# Get uniq variants from the variants passing the 0.25 cutoff
	awk '{ print $1$2$3 }' pseind_passing_0.25_cutoff.txt | uniq -u > pseind_passing_0.25_uniq.txt
	cat pseind_passing_0.25_uniq.txt | xargs -n 1 -I % awk '{ rep = $1$2$3; if (rep == "%") print $0}' pseind_passing_0.25_cutoff.txt > pseind_passing_0.25_uniq_full.txt

	# Get uniq variants not passing the 0.25 cutoff but unique in their respective synteny block
	awk '{ print $1$2$3 }' pseind_notpassing_0.25_cutoff.txt | uniq -u > pseind_notpassing_0.25_cutoff_butuniq_id.txt
	cat pseind_notpassing_0.25_cutoff_butuniq_id.txt | xargs -n 1 -I % awk '{ rep = $1$2$3; if (rep == "%") print $0}' pseind_notpassing_0.25_cutoff.txt > pseind_notpassing_0.25_cutoff_butuniq_id_full.txt

	## Get the final dpse dels

	# Concatenate uniq variants and the curated repeated variants
	cat pseind_passing_0.25_repeated_full.txt pseind_passing_0.25_uniq_full.txt > pseref_indels_cross-validated_final.txt

	# Add the uniq not passing the cutoff
	cat pseind_notpassing_0.25_cutoff_butuniq_id_full.txt pseref_indels_cross-validated_final.txt > pseref_indels_cross-validated_including_notpassing0.25.txt

	## Final processing to remove duplicates and transform to bed file
	## there are variants that overlap more than one collinear block.
	awk '{ print $4"\t"$5"\t"$6"\t"$7" "$8" "$1 }' pseref_indels_cross-validated_including_notpassing0.25.txt > dpse_delrefs_validatedwith_dperins.bed
	bedtools sort -i dpse_delrefs_validatedwith_dperins.bed > dpse_delrefs_validatedwith_dperins_sorted.bed
	awk '{ print $0 }' dpse_delrefs_validatedwith_dperins_sorted.bed | uniq -d > final_dpse_dels_repeated.txt
	awk '{ print $0 }' dpse_delrefs_validatedwith_dperins_sorted.bed | uniq -u > dpse_delrefs_validatedwith_dperins_sorted_uniq.txt
	cat dpse_delrefs_validatedwith_dperins_sorted_uniq.txt final_dpse_dels_repeated.txt > dpse_delrefs_validatedwith_dperins_norepeats.bed

	####FINAL FILE for DPSE as reference sorted: bedtools sort -i dpse_delrefs_validatedwith_dperins_norepeats.bed > per2dpse_dels_validated_with_ins_sorted.bed


	##Get the final dper ins

	## Final processing to remove duplicates and transform to bed file
	## there are variants that overlap more than one collinear block.
	awk '{ print $12"\t"$13"\t"$14"\t"$15" "$16" "$1}' pseref_indels_cross-validated_including_notpassing0.25.txt > dper_ins_dpsedelsref.bed
	bedtools sort -i dper_ins_dpsedelsref.bed > dper_ins_dpsedelsref_sorted.bed
	awk '{ print $0 }' dper_ins_dpsedelsref_sorted.bed | uniq -d > final_dper_ins_repeated.txt
	awk '{ print $0 }' dper_ins_dpsedelsref_sorted.bed | uniq -u > dper_ins_dpsedelsref_sorted_uniq.txt
	cat final_dper_ins_repeated.txt dper_ins_dpsedelsref_sorted_uniq.txt > dper_ins_dpsedelsref_norepeats.bed


	####FINAL FILE for DPER dels sorted: bedtools sort -i dper_ins_dpsedelsref_norepeats.bed > per2dpse_ins_usedtovalidate_dels_sorted.bed


### INS in Dper supported by DEL in Dpse. ##################################################################################################################################

	# Format svim output to a simplified bed file
	awk '{ len = $3 - $2; if ($5 > 1) print $1"\t"$2"\t"$3"\t"$4" "$5" "len }' candidates_novel_insertions.bed > allinsertions.bed

	# Sort the resulting bedfile
	bedtools sort -i allinsertions.bed > allinsertions_sorted.bed

	# Copy the sinteny blocks
	cp ../perDEL_analysis_for_paper/pseref2per_collinearblocks_sorted_withID.bed .
	cp ../perDEL_analysis_for_paper/dPERpseref2per_collinearblocks_sorted_withID_sorted.bed .

	# Find indel coordinates on collinear blocks between Dper and Dpse allowing partial overlap to include indels located at the break ends of the synteny blocks
	bedtools intersect -wo -a allinsertions_sorted.bed -b pseref2per_collinearblocks_sorted_withID.bed > allinsertions_collinearblocks.bed

	# Filter variants matching collinearblocks by the svim score (>10 is suggested to be considered a true variant)
	awk '{ if ($5 > 10) print $0}' allinsertions_collinearblocks.bed > passed_insertions_collinearblocks.bed

	# Sort filtered variants
	bedtools sort -i passed_insertions_collinearblocks.bed > passed_insertions_collinearblocks_sorted.bed

	# Copy deletions in Dpse
	cp /media/DATA5/Javier/Dpse_Dper_annotation_ALL/svim/new_ver/notmasked/pse2per/output/candidates/candidates_deletions.bed .

	# Format the insertions file to a simplified bed file
	awk '{ len = $3 - $2; if ($5 > 1) print $1"\t"$2"\t"$3"\t"$4" "$5" "len }' candidates_deletions.bed > alldeletions.bed

	# Sort the insertions bed file
	bedtools sort -i alldeletions.bed > alldeletions_sorted.bed

	# Find indel coordinates on synteny blocks but now on Dper allowing partial overlap to include indels located at the break ends of the synteny blocks
	bedtools intersect -wo -a alldeletions_sorted.bed -b dPERpseref2per_collinearblocks_sorted_withID_sorted.bed > alldeletions_collinearblocks.bed

	# Sort the insertions sorted bed file
	bedtools sort -i alldeletions_collinearblocks.bed > alldeletions_collinearblocks_sorted.bed

	# Copy the get_block_INDEL.pl script
	cp ../perDEL_analysis_for_paper/get_block_INDEL.pl .

	# Get INDELs for each synteny block and create a file for each block (for both INS and DEL).
	perl get_block_INDEL.pl passed_insertions_collinearblocks_sorted.bed alldeletions_collinearblocks_sorted.bed  pseref2per_collinearblocks_sorted_withID.bed

	#create a new directory and move indel files
	mkdir block_INDELS
	mv block-*_indels.txt block_INDELS

	# Copy the transform_coordinates_v1.0.pl script
	cp ../../perDEL_analysis_for_paper/block_INDELS/transform_coordinates_v1.0.pl .

	# Transform coordinates to compare dels to ins and to find matching variants in the two genomes
	for i in $(ls *_indels.txt); do perl transform_coordinates_v1.0.pl $i;done

	# Move all blocks with no Dper variants to blocks-nodper/
	ls -l per*.bed | awk '{ if ($5 == "0") print $9 }' | xargs -n 1 -I % mv % blocks-nodper/

	# Sort remaining bed files
	for i in $(ls per_block*.bed); do bedtools sort -i $i > $i.sorted; done
	for i in $(ls pse_block*.bed); do bedtools sort -i $i > $i.sorted; done

	# create a list containing the blocks with dper variants and replace the 'per'
	ls per*.sorted > blocks.txt
	sed -i 's/per_//' blocks.txt

	# Find collinear variants using bedtools closest
	while read i; do bedtools closest -d -k 2 -a pse_$i -b per_$i > $i.closest ; done < blocks.txt

	# create directory for collinear variants and a list of the closest files
	mkdir matching_blocks_0dist
	ls *.sorted.closest > closests_blocks.txt

	# Copy the find_matching_indels.pl script
	cp ../../perDEL_analysis_for_paper/block_INDELS/find_matching_indels.pl .

	# find matching indels using the find_matching_indels.pl script
	perl find_matching_indels.pl closests_blocks.txt

	# create a directory and move the blocks with no matching variants
	mkdir matching_blocks_far
	mv *closest matching_blocks_far

	# Copy the get_matching_indels.pl script
	cp ../../../perDEL_analysis_for_paper/block_INDELS/matching_blocks_0dist/get_matching_indels.pl .

	# For the matching variants in matching_blocks_0dist, find the exact match using the script get_matching_indels.pl
	for i in $(ls *.closest); do perl get_matching_indels.pl $i > $i.exact;done

	# Concatenate all exact files
	cat *.exact > perdels_exact_match_perdels.txt

	# Copy the get_final_exact.pl script
	cp ../../../perDEL_analysis_for_paper/block_INDELS/matching_blocks_0dist/get_final_exact.pl .

	# Get matching variants that are similiar not only on the position but also in the size allowing 25% error, i.e. del: 70 bp then ins: ~70 bp
	perl get_final_exact.pl perdels_exact_match_perdels.txt 0.25

	# Get duplicated per dels from the variants passing the 0.25 cutoff
	awk '{ print $1$2$3 }' pseind_passing_0.25_cutoff.txt | uniq -d > pseind_passing_0.25_repeated.txt
	cat pseind_passing_0.25_repeated.txt | xargs -n 1 -I % awk '{ rep = $1$2$3; if (rep == "%") print $0}' pseind_passing_0.25_cutoff.txt > pseind_passing_0.25_repeated_full.txt

	#manually curate the file pseind_passing_0.25_repeated_full.txt

	# Get uniq variants from the variants passing the 0.25 cutoff
	awk '{ print $1$2$3 }' pseind_passing_0.25_cutoff.txt | uniq -u > pseind_passing_0.25_uniq.txt
	cat pseind_passing_0.25_uniq.txt | xargs -n 1 -I % awk '{ rep = $1$2$3; if (rep == "%") print $0}' pseind_passing_0.25_cutoff.txt > pseind_passing_0.25_uniq_full.txt

	# Get uniq variants not passing the 0.25 cutoff but unique in their respective synteny block
	awk '{ print $1$2$3 }' pseind_notpassing_0.25_cutoff.txt | uniq -u > pseind_notpassing_0.25_cutoff_butuniq_id.txt
	cat pseind_notpassing_0.25_cutoff_butuniq_id.txt | xargs -n 1 -I % awk '{ rep = $1$2$3; if (rep == "%") print $0}' pseind_notpassing_0.25_cutoff.txt > pseind_notpassing_0.25_cutoff_butuniq_id_full.txt

	## Get the final dpse dels

	# Concatenate uniq variants and the curated repeated variants
	cat pseind_passing_0.25_repeated_full.txt pseind_passing_0.25_uniq_full.txt > pseref_indels_cross-validated_final.txt

	# Add the uniq not passing the cutoff
	cat pseind_notpassing_0.25_cutoff_butuniq_id_full.txt pseref_indels_cross-validated_final.txt > pseref_indels_cross-validated_including_notpassing0.25.txt

	## Final processing to remove duplicates and transform to bed file
	## there are variants that overlap more than one collinear block.
	awk '{ print $4"\t"$5"\t"$6"\t"$7" "$8" "$1 }' pseref_indels_cross-validated_including_notpassing0.25.txt > dpse_delrefs_validatedwith_dperins.bed
	bedtools sort -i dpse_delrefs_validatedwith_dperins.bed > dpse_delrefs_validatedwith_dperins_sorted.bed
	awk '{ print $0 }' dpse_delrefs_validatedwith_dperins_sorted.bed | uniq -d > final_dpse_dels_repeated.txt
	awk '{ print $0 }' dpse_delrefs_validatedwith_dperins_sorted.bed | uniq -u > dpse_delrefs_validatedwith_dperins_sorted_uniq.txt
	cat dpse_delrefs_validatedwith_dperins_sorted_uniq.txt final_dpse_dels_repeated.txt > dpse_delrefs_validatedwith_dperins_norepeats.bed

	####FINAL FILE for DPSE as reference sorted: bedtools sort -i dpse_delrefs_validatedwith_dperins_norepeats.bed > per2dpse_ins_validated_with_dels_sorted.bed


	##Get the final dper ins

	## Final processing to remove duplicates and transform to bed file
	## there are variants that overlap more than one collinear block.
	awk '{ print $12"\t"$13"\t"$14"\t"$15" "$16" "$1 }' pseref_indels_cross-validated_including_notpassing0.25.txt > dper_ins_dpsedelsref.bed
	bedtools sort -i dper_ins_dpsedelsref.bed > dper_ins_dpsedelsref_sorted.bed
	awk '{ print $0 }' dper_ins_dpsedelsref_sorted.bed | uniq -d > final_dper_ins_repeated.txt
	awk '{ print $0 }' dper_ins_dpsedelsref_sorted.bed | uniq -u > dper_ins_dpsedelsref_sorted_uniq.txt
	cat final_dper_ins_repeated.txt dper_ins_dpsedelsref_sorted_uniq.txt > dper_ins_dpsedelsref_norepeats.bed


	####FINAL FILE for DPER dels sorted: bedtools sort -i dper_ins_dpsedelsref_norepeats.bed > per2dpse_del_usedtovalidate_ins_sorted.bed
