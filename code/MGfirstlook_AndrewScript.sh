#!/bin/bash

cd /dfs3/bio/tgallagh/AntibioticsTrial/data/processed
mkdir quality_filtered
mkdir midas
mkdir anvio
mkdir metaphlan

BASEDIR=/dfs3/bio/tgallagh/AntibioticsTrial/data/raw/
QUALITY=/dfs3/bio/tgallagh/AntibioticsTrial/data/processed/quality_filtered/
MIDAS=/dfs3/bio/tgallagh/AntibioticsTrial/data/processed/midas/
METAPHLAN=/dfs3/bio/tgallagh/AntibioticsTrial/data/processed/metaphlan/

cd $BASEDIR
for f in *_R2.fastq.gz; do
g=$(basename $f _R2.fastq.gz);

echo $g

echo "
#!/bin/bash
#$ -N sample_${g}
#$ -q bio,abio*,free*
#$ -pe openmp 8
#$ -ckpt restart

gunzip ${g}_R1.fastq.gz
gunzip ${g}_R2.fastq.gz

# Quality filter your reads
module load prinseq-lite/0.20.4

prinseq-lite.pl -fastq ${BASEDIR}${g}_R1.fastq \
-fastq2 ${BASEDIR}${g}_R2.fastq \
-out_bad null \
-min_qual_mean 30 \
-ns_max_n 0 \
-out_good ${QUALITY}${g}_qc

# Use Bowtie2 to decontaminate the human genome
module load bowtie2

bowtie2 -q -p 8 -x /dfs3/bio/whitesonlab/hg38/hg38 \
-1 ${QUALITY}${g}_qc_1.fastq \
-2 ${QUALITY}${g}_qc_2.fastq \
-U ${QUALITY}${g}_qc_1_singletons.fastq,${QUALITY}${g}_qc_2_singletons.fastq \
--un ${QUALITY}${g}_qc_decon_se.fastq \
--un-conc ${QUALITY}${g}_qc_decon_pe.fastq \
1>${QUALITY}${g}_out 2>${QUALITY}${g}_alignment_stats.txt

# Clean up any big files bowtie makes
rm ${QUALITY}${g}_out

# Metaphlan analysis
module load metaphlan2

export HDF5_DISABLE_VERSION_CHECK=2

cat ${QUALITY}${g}_qc_decon_*.fastq >> ${METAPHLAN}${g}_tmp.fastq

metaphlan2.py ${METAPHLAN}${g}_tmp.fastq --input_type fastq \
--mpa_pkl /dfs3/bio/whitesonlab/metaphlan_database/db_v20/mpa_v20_m200.pkl \
--bowtie2db /dfs3/bio/whitesonlab/metaphlan_database/db_v20/mpa_v20_m200 \
--nproc 8 > ${METAPHLAN}${g}.txt

rm ${METAPHLAN}${g}_tmp.fastq

# Midas analysis

'''module purge
module load enthought_python/7.3.2
export PYTHONPATH=$PYTHONPATH:/data/users/wengland/software/MIDAS
export PATH=$PATH:/data/users/wengland/software/MIDAS/scripts

run_midas.py species \
${MIDAS}${g} \
-d /dfs3/bio/wengland/MIDAS-dbs/midas_db_v1.2 \
-1 ${QUALITY}${g}_qc_decon_pe.1.fastq \
-2 ${QUALITY}${g}_qc_decon_pe.2.fastq \
-t 8 --remove_temp

run_midas.py genes \
${MIDAS}${g} \
-d /dfs3/bio/wengland/MIDAS-dbs/midas_db_v1.2 \
-1 ${QUALITY}${g}_qc_decon_pe.1.fastq \
-2 ${QUALITY}${g}_qc_decon_pe.2.fastq \
-t 8 --remove_temp \
--species_cov 2
'''

gzip ${g}_R1.fastq
gzip ${g}_R2.fastq
gzip ${QUALITY}${g}_qc_decon_pe.1.fastq
gzip ${QUALITY}${g}_qc_decon_pe.2.fastq
gzip ${QUALITY}${g}_qc_decon_se*fastq

" > ${g}_redo.sh

qsub ${g}_redo.sh
done
