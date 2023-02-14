# Last update: Fri Feb 10 19:30:01 2023
# Last update: Fri Feb 10 19:30:01 2023

gunzip setA.fas

# bin directory ----
bin_dir=""
# dir for output
cp=""

time python3 $bin_dir/Step5_makeblastdb.py \
-c makeblastdb \
-i $cp"/GENOMES" \
-o $cp"/DBS" \
-u 16 \
-t nucl


# before start ----
"-outfmt", str("6 qseqid sseqid score"),

"-outfmt", str("6"),


time python3 $bin_dir/Step5_reciprocal_blast.py \
-c blastn \
-i $cp"/VFDB_setA_nt.fas" \
-d $cp"/DBS" \
-o $cp"/blast_op" \
-e 1e-5 \
-u 16 \
-m off
