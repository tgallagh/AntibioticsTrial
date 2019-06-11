#! /bin/bash
#$ -N fungal_blast
#$ -q pub64,bio
#$ -m e
#$ -pe make 8 
#$ -R y
#$ -t 1-25
module load  blast/2.8.1
#Andrew's fungal + phage + other things database
DB="/dfs3/bio/tgallagh/databases/fungal_database/fb_database.fna"
#STARTING DIRECTORY WITH INPUT DATA
START_DIR="/dfs3/bio/tgallagh/AntibioticsTrial/data/processed/quality_filtered/"
#output directory
OUT_DIR="/dfs3/bio/tgallagh/AntibioticsTrial/data/processed/fungal_blast_output/"


cd $START_DIR
#make prefix file
find . -type f -name "*qc_decon_se.fastq" > temp.txt
sed -i "s/"_qc_decon_se.fastq"//g"  temp.txt
cut -c 3-  temp.txt > prefix.txt
rm temp.txt

input=$(head -n $SGE_TASK_ID prefix.txt | tail -n 1)

#put mates and singletons into one file
cat $input\_qc_decon*.fastq > $input\.cat.fastq
# convert to fasta
sed '/^@/!d;s//>/;N'  $input\.cat.fastq  > $input\.cat.fasta


blastn -query $START_DIR$input\.cat.fasta  -db $DB -out $OUT_DIR$input\output.txt -outfmt "6 qseqid sseqid evalue pident stitle"
