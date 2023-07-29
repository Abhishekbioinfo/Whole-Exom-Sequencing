########################################################################################################################
#	This shell script takes input as folder name of Tumor.								#
#	e.g WES_pipeline_Tumor_hg38_v0.1.sh Tumor 									#
#	Last Modified on 20/07/2023			
#    WES pipeline for somatic
#	 Karkinos-Bioinformatics
#    Date of Implementation 19th July 2023.							                	#
#	
#	
#															#
#															#
#															#
#	Tools Used in this pipeline											#
#	1.  Fastqc													#
#	2.  Trimmomatic													#
#	3.  BWA	mem													#
#	4.  MarkDuplicates							       																		#
#       5.  Samtools                                                                                            	#
#       6.  BQSR													#
#	7.  Mutect2													#	
#	8.  Varscan													#
#	9. VarDict													#
#	10. Lofreq
#	
#
#########################################################################################################################


## Taking user inputs and setting resource paths
input1="$1"     # Taking input of sample folder name as argument
workdir=$(pwd)  # storing path of Working directory
raw_tumor="$workdir/$input1"   # Assigning path of sample
out="$workdir/Output_WES_hg38_somatic_final_test_${input1}"       # Creating and Assigning path of output folder
db="/mnt/database"              # Assigning database folder path
hg38="$db/GRCH38.P14/hg38.fa"   # Assigning path of human genome reference file

#./email.sh.x "$input1 Run Submitted" "WES panel job Sumitted for $input1" 

# Merging fastq files and Creating list file
cd $input1

fn=$(ls | head -1| cut -f 1-2 -d"_")

echo $fn

zcat *R1* > ${fn}_All_R1_001.fastq
zcat *R2* > ${fn}_All_R2_001.fastq

ls *All_R1* > ../list_${input1}

cd $workdir

input2="list_${input1}"

echo "Test list filename" $input2
#Dynamic allocation of cpus

t="$(nproc --all)"      # Fetching total number of cpus
tnp=$((`expr $t - 44`)) # Maximum number of cpus will be utilized
val=$((`expr $tnp / 2`)) # Arithmatic operation on Max cpu count
np=$((`printf "%.*f\n" 0 $val`)) # Number of shared cpu

echo "test" $t $tnp $val $np

eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs_gatk4                # Activating conda environment
# Removing pre-exiting files if any

rm $out/WES_somatic_analysis.log
rm $out/parallel_variantCaller
#rm $out/parallel_variantCaller_germline

echo "First $input1 Second $input2" # Print inputs as a test

# Creating folders for the intermediate results
mkdir -p $out/$input1/FastQC_output
mkdir -p $out/$input1/trimmomatic_output
mkdir -p $out/$input1/FastQC_output_after-trimmomatic
mkdir -p $out/$input1/alignments_stats
mkdir -p $out/$input1/bqsr
mkdir -p $out/$input1/markdup
mkdir -p $out/variantCaller/Varscan
mkdir -p $out/variantCaller/Vardict
mkdir -p $out/variantCaller/Mutect2
mkdir -p $out/variantCaller/lofreq


# Fastqc and Trimmomatics run for Normal
while IFS= read -r line			# While loop to add scripts for each lane in "parallel_fastqc_Normal" 
do					# "parallel_trimmomatic_Normal" and "parallel_after_trimmomatic_Normal"
	echo $line
	header=$(cat $raw_tumor/$line | head -n 1);  # Extracting header from fastq file. Will not work if fastq file
                                                        # is not gunzip file. Change command to cat
	
	id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');  # Extracting run ID from fastq file
	
	echo "header content" $header # printing header
	sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-2 -d"_");  # Extracting file name from file list
	echo $id
	echo $sm

done < $input2
# QC checking

echo "Fastqc for $input1 Started" >> $out/WES_somatic_analysis.log   # Writing in log file
date >> $out/WES_somatic_analysis.log                                # Writing in log file

fastqc $raw_tumor/${sm}_All_R1*  -t $tnp -o $out/$input1/FastQC_output/

fastqc $raw_tumor/${sm}_All_R2*  -t $tnp -o $out/$input1/FastQC_output/

echo "Fastqc for $input1 Completed" >> $out/WES_somatic_analysis.log         # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

# Adaptor Trimming and cleaning

echo "Trimmomatic for $input1 Started" >> $out/WES_somatic_analysis.log      # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file

trimmomatic PE -threads $tnp -phred33 $raw_tumor/${sm}_All_R1* $raw_tumor/${sm}_All_R2* $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50 -trimlog $out/$input1/trimmomatic_output/${sm}_trimlog.txt

echo "Trimmomatic for $input1 Completed" >> $out/WES_somatic_analysis.log      # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

# QC-rechecking after trimming
echo "Fastqc after trimmomatic for $input1 Started" >> $out/WES_somatic_analysis.log         # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file

fastqc -t $tnp $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input1/FastQC_output_after-trimmomatic/

echo "Fastqc after trimmomatic for $input1 Completed" >> $out/WES_somatic_analysis.log # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

#Fastqc_summary

unzip $out/$input1/FastQC_output_after-trimmomatic/${sm}_R1-trimmed_P_fastqc.zip -d $out/$input1/FastQC_output_after-trimmomatic/
unzip $out/$input1/FastQC_output_after-trimmomatic/${sm}_R2-trimmed_P_fastqc.zip -d $out/$input1/FastQC_output_after-trimmomatic/

if
        grep 'PASS	Basic Statistics' $out/$input1/FastQC_output_after-trimmomatic/${sm}_R1-trimmed_P_fastqc/summary.txt && grep 'PASS	Basic Statistics' $out/$input1/FastQC_output_after-trimmomatic/${sm}_R2-trimmed_P_fastqc/summary.txt
then
        echo Fastqc is pass   >> $out/FASTQC_validation.log 
else
        echo Fastqc failed    >> $out/FASTQC_validation.log
        exit 1
fi


# BWA and read group attachment

while IFS= read -r line			# While loop to add scripts for each lane in "parallel_BWA" 
do	

	echo $line
	header=$(cat $raw_tumor/$line | head -n 1);	# Extracting header from fastq file. Will not work if fastq file
                                                        # is not gunzip file. Change command to cat

	id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');	# Extracting run ID from fastq file
	echo "header content" $header	# printing header
	sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-2 -d"_");	# Extracting file name from file 
	echo $id
	echo $sm

echo "BWA for $input1 Started" >> $out/WES_somatic_analysis.log      # Writing in log file
date >> $out/WES_somatic_analysis.log                                # Writing in log file

	if [ -z "$(ls -A  $out/$input1/trimmomatic_output/)" ]		# If statement
	then

	bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:WES\tSM:$sm\tPI:200" $hg38 $raw_tumor/${sm}_All_R1_*.fastq $raw_tumor/${sm}_All_R2_*.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam					                                                            # Will execute when 'trimmomatic_output' folder is empty

	else

	bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:WES\tSM:$sm\tPI:200" $hg38 $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam	# Will execute when condition is false

	fi

done < "$input2"

echo "BWA for $input1 Completed" >> $out/WES_somatic_analysis.log      # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

##location of sorted bam
echo "bwa" >> $out/File_location.log
echo "$out/$input1/${sm}.sorted.bam" >> $out/File_location.log


# Markduplicate
echo "Markduplicate for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log

gatk --java-options "-Xmx30g" MarkDuplicates -I $out/$input1/${sm}.sorted.bam -O $out/$input1/markdup/${sm}_dedup.sorted.bam -M $out/$input1/markdup/${sm}_metrics.txt --MAX_RECORDS_IN_RAM 1000000000 

echo "Markduplicate for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

echo "Samtool Alignment Statistics for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file

sambamba flagstat -t $tnp $out/$input1/${sm}.sorted.bam   > $out/$input1/alignments_stats/${sm}_sorted.txt   

sambamba flagstat -t $tnp  $out/$input1/markdup/${sm}_dedup.sorted.bam > $out/$input1/alignments_stats/${sm}_dedup.sorted.txt	

# Getting alignment statistics

echo "Samtool Alignment Statistics for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

# Generating on_target bam and indexing it
samtools view -b -h -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed $out/$input1/markdup/${sm}_dedup.sorted.bam > $out/$input1/markdup/${sm}_on_target.bam 

samtools index -@ $tnp $out/$input1/markdup/${sm}_dedup.sorted.bam

#Baserecaliberation
echo "Baserecalibrator for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar BaseRecalibrator --input $out/$input1/markdup/${sm}_dedup.sorted.bam  --reference $hg38 --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.vcf --output $out/$input1/bqsr/${sm}_recal_before_data.table -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed  -CPB 50000

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar ApplyBQSR -R $hg38 -I $out/$input1/markdup/${sm}_dedup.sorted.bam -bqsr $out/$input1/bqsr/${sm}_recal_before_data.table -O $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed -CPB 50000

echo "Baserecalibrator for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

##location of bqsr bam
echo "bqsr" >> $out/File_location.log
echo "$out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam" >> $out/File_location.log

echo "Samtool indexing for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file

samtools index -@ $tnp $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam

echo "Samtool indexing for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

echo "Samtools mpileup for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file

samtools mpileup -B -f $hg38 $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam > $out/$input1/bqsr/${sm}_mpileup.bqsr.dedup.sorted.bam

echo "Samtools mpileup for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file
#Variant caller
echo "somatic variant calling for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file 
date >> $out/WES_somatic_analysis.log

echo $sm > sample_list.txt

#VarScan
echo "Varscan for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file 
date >> $out/WES_somatic_analysis.log

echo "java  -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar  /home/genomics/anaconda3/envs/wgs_gatk4/share/varscan-2.3.7-4/VarScan.jar mpileup2cns $out/$input1/bqsr/${sm}_mpileup.bqsr.dedup.sorted.bam --vcf-sample-list sample_list.txt --min-coverage 10 --min-reads2 5 --output-vcf 1 --variants --regions-file /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome_varscan.bed > $out/variantCaller/Varscan/${input1}_varscan.vcf" >> $out/parallel_variantCaller

echo "Varscan for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log
echo "##############################" >> $out/WES_somatic_analysis.log

#Mutect2
echo "Mutect2 for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file 
date >> $out/WES_somatic_analysis.log

echo "java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp  -Xms30g -Xmx30g  -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar Mutect2 -R $hg38 -I $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.interval_list --native-pair-hmm-threads $tnp --germline-resource /mnt/database/GRCH38.P14/gnomAD/af-only-gnomad.hg38.vcf.gz -O $out/variantCaller/Mutect2/${input1}_Mutect2.vcf.gz" >> $out/parallel_variantCaller

echo "Mutect2 for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log
echo "##############################" >> $out/WES_somatic_analysis.log

#VarDict
echo "Vardict for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file 
date >> $out/WES_somatic_analysis.log

echo "java  -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/software/VarDictJava/build/libs/VarDict-1.8.3.jar -r 5 -G $hg38 -f 0.01 -N WES -b $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam -c 1 -S 2 -E 3 -g 4 /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed | /home/genomics/software/VarDict/var2vcf_valid.pl -A -N $sm -d 10 > $out/variantCaller/Vardict/${input1}_vardict.vcf" >> $out/parallel_variantCaller

echo "Vardict for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log
echo "##############################" >> $out/WES_somatic_analysis.log

#lofreq
echo "lofreq for $input1 Started" >> $out/WES_somatic_analysis.log    # Writing in log file 
date >> $out/WES_somatic_analysis.log

echo "lofreq call-parallel --pp-threads $tnp  --call-indels -l /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed -s -S /mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -f $hg38 -o $out/variantCaller/lofreq/${input1}_lofreq.vcf $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam" >> $out/parallel_variantCaller

echo "lofreq for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log
echo "##############################" >> $out/WES_somatic_analysis.log

parallel -j 1 < $out/parallel_variantCaller

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp  -Xms30g -Xmx30g  -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar FilterMutectCalls -R $hg38 -V $out/variantCaller/Mutect2/${input1}_Mutect2.vcf.gz -O $out/variantCaller/Mutect2/${input1}_filtered_Mutect2.vcf.gz

#Reformatting lofreq vcf
Rscript lofreqReformat.R $out/variantCaller/lofreq/${input1}_lofreq.vcf 

echo "somatic variant calling for $input1 Completed" >> $out/WES_somatic_analysis.log    # Writing in log file
date >> $out/WES_somatic_analysis.log
echo "##############################" >> $out/WES_somatic_analysis.log       # Writing in log file

##location of vcf of somatic variant callers
#Mutect2
echo "Mutect2" >> $out/File_location.log
echo "$out/variantCaller/Mutect2/${input1}_filtered_Mutect2.vcf.gz" >> $out/File_location.log
#Vardict
echo "Vardict" >> $out/File_location.log
echo "$out/variantCaller/Vardict/${input1}_vardict.vcf"  >> $out/File_location.log
#Varscan 
echo "Varscan" >> $out/File_location.log
echo "$out/variantCaller/Varscan/${input1}_varscan.vcf" >> $out/File_location.log
#lofreq
echo "lofreq" >> $out/File_location.log
echo "$out/variantCaller/lofreq/${input1}_lofreq_rf.vcf" >> $out/File_location.log

#print version
fastqc --version >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo "Trimmomatic"   >> $out/WES_version.log
trimmomatic -version >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo "BWA 0.7.17-r1188" >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo "gatk4-4.2.6.1-1" >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo "sambamba 0.8.2" >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo "samtools 1.6" >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo "VarScan v2.3" >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo "VarDict 1.8" >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log
echo  "LoFreq version 2" >> $out/WES_version.log
echo "##############################" >> $out/WES_version.log

#./email2.sh.x "$input1 Run Completed" "WES panel job Completed for $input1" "$out/WES_somatic_analysis.log" "$out/File_location.


