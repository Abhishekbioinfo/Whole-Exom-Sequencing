## Taking user inputs and setting resource paths
input1="$1"     # Taking input of sample folder name as argument
workdir=$(pwd)  # storing path of Working directory
raw_blood="$workdir/$input1"   # Assigning path of sample
out="$workdir/Output_WES_germline_hg38_Final_test_${input1}"       # Creating and Assigning path of output folder
db="/mnt/database"              # Assigning database folder path
hg38="$db/GRCH38.P14/hg38.fa"   # Assigning path of human genome reference file

#./email.sh.x "$input1 Run Submitted" "K75 panel job Sumitted for $input1"

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

# Removing pre-exiting files if any

rm $out/WES_germline_analysis.log
rm $out/parallel_variantCaller

echo "First $input1 Second $input2" # Print inputs as a test

# Creating folders for the intermediate results
mkdir -p $out/$input1/FastQC_output
mkdir -p $out/$input1/trimmomatic_output
mkdir -p $out/$input1/FastQC_output_after-trimmomatic
mkdir -p $out/$input1/alignments_stats
mkdir -p $out/$input1/markdup
mkdir -p $out/$input1/bqsr
mkdir -p $out/variantCaller/HaplotypeCaller
mkdir -p $out/variantCaller/Freebayes
mkdir -p $out/variantCaller/Platypus
mkdir -p $out/variantCaller/Strelka2

eval "$(conda shell.bash hook)"		# Setting bash for conda environment
conda activate wgs_gatk4		# Activating conda environment

# Fastqc and Trimmomatics run 
while IFS= read -r line			# While loop to add scripts for each lane in "parallel_fastqc_Normal" 
do					# "parallel_trimmomatic_Normal" and "parallel_after_trimmomatic_Normal"
	echo $line
	header=$(cat $raw_blood/$line | head -n 1);  # Extracting header from fastq file. Will not work if fastq file

                                                        # is not gunzip file. Change command to cat
	
	id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');  # Extracting run ID from fastq file
	
	echo "header content" $header # printing header
	sm=$(echo $raw_blood/$line | xargs -n 1 basename | cut -f 1-2 -d"_");  # Extracting file name from file 
        echo $id
	echo $sm

done < $input2

# QC checking

echo "Fastqc for $input1 Started" >> $out/WES_germline_analysis.log   # Writing in log file
date >> $out/WES_germline_analysis.log                                # Writing in log file

fastqc $raw_blood/${sm}_All_R1*  -t $tnp -o $out/$input1/FastQC_output/

fastqc $raw_blood/${sm}_All_R2*  -t $tnp -o $out/$input1/FastQC_output/

echo "Fastqc for $input1 Completed" >> $out/WES_germline_analysis.log         # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log       # Writing in log file


#Adaptor Trimming and cleaning

echo "Trimmomatic for $input1 Started" >> $out/WES_germline_analysis.log      # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file

/home/genomics/anaconda3/envs/wgs_gatk4/bin/java -Xmx64g -Xmx64g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads $tnp -phred33 $raw_blood/${sm}_All_R1* $raw_blood/${sm}_All_R2* $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50 -trimlog $out/$input1/trimmomatic_output/${sm}_trimlog.txt

echo "Trimmomatic for $input1 Completed" >> $out/WES_germline_analysis.log      # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log       # Writing in log file

# QC-rechecking after trimming

echo "Fastqc after trimmomatic for $input1 Started" >> $out/WES_germline_analysis.log         # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file

fastqc -t $tnp $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input1/FastQC_output_after-trimmomatic/

echo "Fastqc after trimmomatic for $input1 Completed" >> $out/WES_germline_analysis.log # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log       # Writing in log file

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

# BWA 

while IFS= read -r line			# While loop to add scripts for each lane in "parallel_BWA" 
do	

	echo $line
	header=$(cat $raw_blood/$line | head -n 1);	# Extracting header from fastq file. Will not work if fastq file
                                                        # is not gunzip file. Change command to cat

	id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');	# Extracting run ID from fastq file
	echo "header content" $header	# printing header
	sm=$(echo $raw_blood/$line | xargs -n 1 basename | cut -f 1-2 -d"_");	# Extracting file name from file 
	echo $id
	echo $sm

echo "BWA for $input1 Started" >> $out/WES_germline_analysis.log      # Writing in log file
date >> $out/WES_germline_analysis.log                                # Writing in log file

	if [ -z "$(ls -A  $out/$input1/trimmomatic_output/)" ]		# If statement
	then

	bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:WES\tSM:$sm\tPI:200" $hg38 $raw_blood/${sm}_All_R1_*.fastq $raw_blood/${sm}_All_R2_*.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam      
                            # Will execute when 'trimmomatic_output' folder is empty

        else

      bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:WES\tSM:$sm\tPI:200" $hg38 $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam    # Will execute when condition is false

	fi

done < "$input2"

echo "BWA for $input1 Completed" >> $out/WES_germline_analysis.log      # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log       # Writing in log file

##location of sorted bam 
echo "bwa" >> $out/File_location.log
echo "$out/$input1/${sm}.sorted.bam" >> $out/File_location.log

# Markduplicate

echo "Markduplicate for $input1 Started" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log

gatk --java-options "-Xmx64g" MarkDuplicates -I $out/$input1/${sm}.sorted.bam -O $out/$input1/markdup/${sm}_dedup.sorted.bam -M $out/$input1/markdup/${sm}_metrics.txt --MAX_RECORDS_IN_RAM 1000000000

echo "Markduplicate for $input1 Completed" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log       # Writing in log file

#Samtool Alignment Statistics

echo "Samtool Alignment Statistics for $input1 Started" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file

sambamba flagstat -t $tnp $out/$input1/${sm}.sorted.bam   > $out/$input1/alignments_stats/${sm}_sorted.txt   

sambamba flagstat -t $tnp  $out/$input1/markdup/${sm}_dedup.sorted.bam > $out/$input1/alignments_stats/${sm}_dedup.sorted.txt	# Getting alignment statistics

echo "Samtool Alignment Statistics for $input1 Completed" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log       # Writing in log file

#Generating on_target bam and indexing it

samtools view -b -h -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed $out/$input1/markdup/${sm}_dedup.sorted.bam > $out/$input1/markdup/${sm}_on_target.bam

samtools index -@ $tnp $out/$input1/markdup/${sm}_dedup.sorted.bam

#Baserecaliberation

echo "Baserecalibrator for $input1 Started" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar BaseRecalibrator --input $out/$input1/markdup/${sm}_dedup.sorted.bam  --reference $hg38 --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.vcf --output $out/$input1/bqsr/${sm}_recal_before_data.table -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed -CPB 50000

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar ApplyBQSR -R $hg38 -I $out/$input1/markdup/${sm}_dedup.sorted.bam -bqsr $out/$input1/bqsr/${sm}_recal_before_data.table -O $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed -CPB 50000

echo "Baserecalibrator for $input1 Completed" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log                                        # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log       # Writing in log file

##location of bqsr bam 

echo "BQSR" >> $out/File_location.log
echo "$out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam" >> $out/File_location.log

#Variant_Caller

echo "germline variant calling for $input1 started" >> $out/WES_germline_analysis.log
date >> $out/WES_germline_analysis.log

#HaplotypeCaller
echo "HaplotypeCaller for $input1 Started" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar HaplotypeCaller -R $hg38 -I $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam --dbsnp $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -O $out/variantCaller/HaplotypeCaller/${sm}_HaplotypeCaller.vcf  --native-pair-hmm-threads $tnp -L /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.interval_list

echo "HaplotypeCaller for $input1 Completed" >> $out/WES_germline_analysis.log  # Writing in log file
date >> $out/WES_germline_analysis.log                           # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log

conda deactivate

eval "$(conda shell.bash hook)"
conda activate rna_preprocessing

#Platypus
echo "Platypus for $input1 Started" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log

platypus callVariants --nCPU $tnp --bamFiles=$out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam --refFile=$hg38 --regions=/home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed  --output=$out/variantCaller/Platypus/${sm}_platypus.vcf.gz

echo "Platypus for $input1 Completed" >> $out/WES_germline_analysis.log  # Writing in log file
date >> $out/WES_germline_analysis.log                           # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log

#Strelka2
echo "Strelka2 germline analysis Started" >> $out/WES_germline_analysis.log  # Writing in log file
date >> $out/WES_germline_analysis.log                           # Writing in log file

/home/genomics/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam --referenceFasta $hg38 --runDir $out/variantCaller/Strelka2 --exome --callRegions /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed.gz

python $out/variantCaller/Strelka2/runWorkflow.py -m local -j $tnp
mv $out/variantCaller/Strelka2/results/variants/variants.vcf.gz $out/variantCaller/Strelka2/results/variants/${sm}_Strelka2.vcf.gz

echo "Strelka2 germline analysis Completed" >> $out/WES_germline_analysis.log  # Writing in log file
date >> $out/WES_germline_analysis.log                           # Writing in log file
echo "##############################" >> $out/WES_germline_analysis.log

conda deactivate

#Freebayes
echo "Freebayes for $input1 Started" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log

freebayes -f $hg38 -t /home/genomics/bedfile/WES_Clinical_IDT_Panel/exome.bed -b $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam  > $out/variantCaller/Freebayes/${sm}_freebayes.vcf

echo "Freebayes for $input1 Completed" >> $out/WES_germline_analysis.log    # Writing in log file
date >> $out/WES_germline_analysis.log
echo "##############################" >> $out/WES_germline_analysis.log

echo "germline variant calling for $input1 completed" >> $out/WES_germline_analysis.log
date >> $out/WES_germline_analysis.log
echo "##############################" >> $out/WES_germline_analysis.log

##location of vcf of germline variant callers

#Haplotypecaller
echo "HaplotypeCaller" >> $out/File_location.log	
echo "$out/variantCaller/HaplotypeCaller/${sm}_HaplotypeCaller.vcf" >> $out/File_location.log
#Platypus
echo "Platypus" >> $out/File_location.log
echo "$out/variantCaller/Platypus/${sm}_platypus.vcf.gz"  >> $out/File_location.log
#Freebayes
echo "Freebayes" >> $out/File_location.log
echo "$out/variantCaller/Freebayes/${sm}_freebayes.vcf" >> $out/File_location.log
#Strelka
echo "Strelka2" >> $out/File_location.log
echo "$out/variantCaller/Strelka2/results/variants/${sm}_Strelka2.vcf.gz" >> $out/File_location.log


eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs_gatk4                # Activating conda environment

# print the version of the tools in the log file 
fastqc --version >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "Trimmomatic"   >> $out/WES_germline_version.log
trimmomatic -version >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "BWA 0.7.17-r1188" >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "gatk4-4.2.6.1-1" >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "sambamba 0.8.2" >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "samtools 1.6" >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "Haplotype_gatk-4.2.6.1" >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "Platypus-0.8.1.2" >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "strelka2-2.9.2" >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log
echo "Freebayes-0.9.21 " >> $out/WES_germline_version.log
echo "##############################" >> $out/WES_germline_version.log

#./email2.sh.x "$input1 Run Completed" "K75 panel job Completed for $input1" "$out/WES_germline_analysis.log" "$out/File_location.log"

