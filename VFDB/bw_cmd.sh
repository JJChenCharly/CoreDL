# index building
time bowtie2-build -f /home/docker_qiime2_share_D/biodb/VFDB/VFDB_setA_nt.fas \
/home/docker_qiime2_share_D/91_EC/VFDB/bw/ind_bw/set_A_nt_inde --threads 36 > /home/docker_qiime2_share_D/91_EC/VFDB/bw/bw_ind_build_log.txt 2>&1

# test
time bowtie2 -q -p 36 --no-unal -N 0 -L 22 -i S,1,0.50 --n-ceil L,1,0 --end-to-end -a -D 22 -R 3 -I 200 -X 700 -t \
-x /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/ind_base_75/test_index_75 \
-1 /home/docker_qiime2_share_D/91_EC/fastp/server/op/A_1/A_1-fastp.1.fastq \
-2 /home/docker_qiime2_share_D/91_EC/fastp/server/op/A_1/A_1-fastp.2.fastq \
--met-file /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_met.txt \
-S /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op.sam > /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op_log.txt 2>&1

# another way index building
    foo="Hello"
    foo="${foo},World"
    echo "${foo}"


    ind_lst=""

    for x in `ls /home/docker_qiime2_share_D/91_EC/VFDB/bw_ind/VF75`
    do ind_lst="${ind_lst}'/home/docker_qiime2_share_D/91_EC/VFDB/bw_ind/VF75/$x',"
    done

    ind_lst=${ind_lst::-1}

    echo ${ind_lst:0:85}

    $ind_lst > 75.txt 2>&1

    #!!!!!
    echo "bowtie2-build --threads 36 -f $ind_lst /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/ind_base/test_index" > /home/docker_qiime2_share_D/91_EC/VFDB/bw/haha.sh 2>&1

    sh /home/docker_qiime2_share_D/91_EC/VFDB/bw/haha.sh
    !!!!
    bowtie2-build -f '/home/docker_qiime2_share_D/91_EC/VFDB/bw_ind/VF75/VFG000358(gb|WP_002212885).fasta' \
    /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/ind_base/test_index --threads 36 --no-unal

# bowtie2 ----

# samtools view
samtools view --threads 36 -bS /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op.sam > /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op.bam
samtools view --threads 36 -bS /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op.sam --output-fmt bam -o \
/home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op2.bam > /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op.bam


samtools sort -u --output-fmt bam --threads 36 /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op.bam -o /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op_sorted.bam
samtools coverage -l 100 /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op_sorted.bam -o /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/test_op_coverage.txt

####################actuall for all
date ; time for x in `ls /home/docker_qiime2_share_D/91_EC/fastp/server/op`
do 
echo "<<<" $x
time bowtie2 -q -p 36 --no-unal -N 0 -L 22 -i S,1,0.50 --n-ceil L,1,0 --end-to-end -a -D 22 -R 3 -I 200 -X 700 -t \
-x /home/docker_qiime2_share_D/91_EC/VFDB/bw/ind_bw/set_A_nt_inde \
-1 /home/docker_qiime2_share_D/91_EC/fastp/server/op/$x/$x-fastp.1.fastq \
-2 /home/docker_qiime2_share_D/91_EC/fastp/server/op/$x/$x-fastp.2.fastq \
-S /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_of_91/$x.sam
done > /home/docker_qiime2_share_D/91_EC/VFDB/bw/log_of_91_bw.txt 2>&1 ; date

--met-file /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_all/$x-met.txt \

# for all unpaired after fastp
    date ; time for x in `ls /home/docker_qiime2_share_D/91_EC/fastp/server/op/`
    do 
    bowtie2 -q -p 36 --no-unal -N 0 -L 22 -i S,1,0.50 --n-ceil L,1,0 --end-to-end -a -D 22 -R 3 -I 200 -X 700 -t \
    -x /home/docker_qiime2_share_D/91_EC/VFDB/bw/test/ind_base_75/test_index_75 \
    -U /home/docker_qiime2_share_D/91_EC/fastp/server/op/$x/$x-fastp.cated-un.fastq \
    -S /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_all/$x-un.sam
    done > /home/docker_qiime2_share_D/91_EC/VFDB/bw/log_of_unpaired_in.txt 2>&1 ; date

95 min

--met-file /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_all/$x-met.txt \
##### veiw all
for x in `ls /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_of_91 | cut -d "." -f 1`
do
samtools view -bS --threads 36 /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_of_91/$x.sam > \
/home/docker_qiime2_share_D/91_EC/VFDB/bw/op_view_all/$x-view.bam
done

# sort all
# for x in `ls /home/docker_qiime2_share_E/91_EC/VFDB/bw/op_all_view/ | cut -d "." -f 1`
# do
# samtools sort -u /home/docker_qiime2_share_E/91_EC/VFDB/bw/op_all_view/$x.bam \
# -o /home/docker_qiime2_share_E/91_EC/VFDB/bw/op_all_sort/$x-sort.bam
# done

parallel -j 36 samtools sort -u /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_view_all/{}-view.bam \
-o /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_sort_all/{}-sort.bam \
::: `ls /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_view_all | cut -d "-" -f 1`

# coverage
# for x in `ls /home/docker_qiime2_share_E/91_EC/VFDB/bw/op_all_sort/ | cut -d "." -f 1`
# do
# samtools coverage -l 100 /home/docker_qiime2_share_E/91_EC/VFDB/bw/op_all_sort/$x.bam \
# -o /home/docker_qiime2_share_E/91_EC/VFDB/bw/op_all_coverage/$x-cov.txt
# done

# -l  ignore reads shorter than INT bp [0]

parallel -j 36 samtools coverage -l 100 /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_sort_all/{}-sort.bam \
-o /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_cov_all/{}-cov.txt \
::: `ls /home/docker_qiime2_share_D/91_EC/VFDB/bw/op_sort_all/ | cut -d "-" -f 1`
