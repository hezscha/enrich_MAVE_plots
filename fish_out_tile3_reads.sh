#check which shell is running
echo $SHELL

#I need to source in order to be able to use conda activate. Either source bashrc if calling this script with 'bash run_domain_analysis.sh' or source .zshrc if running with 'zsh run_domain_analysis.sh'
#source /home/hezscha/.zshrc
source /home/hezscha/.bashrc
#reference: https://trello.com/c/foieoBrq/72-current-procedure
#conda activate p3-env-cyvcf2
conda activate bwa

sel=$1

for tp in {0..3..1}
	do
	sample=3${sel}${tp}

	tmp_dir=/storage1/shared/data/dhfr_mave/fulldhfr_oct2021_quality_controll/tile3_map/temp/
	dest_dir=/storage1/shared/data/dhfr_mave/fastq_fulldhfr_oct2021/${sample}_merge/
	#cd /storage1/shared/data/dhfr_mave/fulldhfr_oct2021_quality_controll/tile3_map/

	echo $sample
	echo $dest_dir

	#map
	echo 'Aligning reads'
	bwa mem /storage1/shared/data/dhfr_mave/fulldhfr_oct2021_quality_controll/tile3_map/tile3_wt_fw.fasta /storage1/shared/data/dhfr_mave/fastq_fulldhfr_oct2021/${sample}_merge/${sample}_R1.fastq.gz /storage1/shared/data/dhfr_mave/fastq_fulldhfr_oct2021/${sample}_merge/${sample}_R2.fastq.gz > ${tmp_dir}/${sample}_both.sam

	#convert to bam
	echo 'covert to bam'
	samtools view -b ${tmp_dir}/${sample}_both.sam > ${tmp_dir}/${sample}_both.bam
	rm ${tmp_dir}/${sample}_both.sam

	#filter to only mapped
	echo 'filter'
	samtools view -u -f 2 -F 4 ${tmp_dir}/${sample}_both.bam | samtools sort -n -o ${tmp_dir}/${sample}_both.mapped.sort.bam

	#make to fastq
	echo 'making fastq'
	samtools fastq -1 ${dest_dir}/${sample}_R1.mapped.fastq -2 ${dest_dir}/${sample}_R2.mapped.fastq -s ${tmp_dir}/${sample}_singleton.mapped.fastq -0 ${tmp_dir}/${sample}_ambiguos.mapped.fastq ${tmp_dir}/${sample}_both.mapped.sort.bam

	#zip
	gzip ${dest_dir}/*fastq

	#trimmomatic
	conda activate cutadapt

	trimmomatic PE ${dest_dir}/${sample}_R1.fastq.gz ${dest_dir}/${sample}_R2.fastq.gz ${dest_dir}/${sample}_R1.trimmomatic.paired.fastq.gz ${dest_dir}/${sample}_R1.trimmomatic.unpaired.fastq.gz ${dest_dir}/${sample}_R2.trimmomatic.paired.fastq.gz ${dest_dir}/${sample}_R2.trimmomatic.unpaired.fastq.gz MINLEN:151 CROP:151 2> /storage1/shared/data/dhfr_mave/fulldhfr_oct2021_quality_controll/read_lengths/${sample}.trimmomatic.stats

	cutadapt -a XXX -A XXX --action none -m 151 -M 151 -o ${dest_dir}/${sample}_R1_length_filtered.fastq.gz -p ${dest_dir}/${sample}_R2_length_filtered.fastq.gz ${dest_dir}/${sample}_R1.mapped.fastq.gz ${dest_dir}/${sample}_R2.mapped.fastq.gz
done
	
