# format svim variants
awk '{ print $1"\t"$2"\t"$3"\t"$4" "$5 }' SVper2pse_svmuconfirmed.txt > SV_svmuper2pse.txt

# Transform coordinates for perDel_per2pse.bed and pseIns_per2pse.bed
bash transform_coordinates_dper.sh variants_dper.bed
bash transform_coordinates_dpse.sh variants_dpse.bed

# generate final polarization file
perl polarization.pl Transformed_pseSVs_per2pse.bed Transformed_perSVs_per2pse.bed > perpseSVs_dmirpolarization.txt
