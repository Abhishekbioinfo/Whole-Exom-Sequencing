## Taking user inputs and setting resource paths
input1="$1"     # Taking input of sample folder name as argument
workdir=$(pwd)  # storing path of Working directory
raw_tumor="$workdir/$input1"   # Assigning path of sample
out="$workdir/Output_K75_CNV_hg38_${input1}"       # Creating and Assigning path of output folder
db="/mnt/database"              # Assigning database folder path
hg38="$db/GRCH38.P14/hg38.fa"   # Assigning path of human genome reference file

# Merging fastq files and Creating list file
cd $input1

fn=$(ls | head -1| cut -f 1-3 -d"_")
echo $fn

zcat *R1* > ${fn}_All_R1_001.fastq
zcat *R2* > ${fn}_All_R2_001.fastq

ls *All_R1* > ../list_${input1}

cd $workdir

input2="list_${input1}"

#echo "Test list filename" $input2
#Dynamic allocation of cpus

t="$(nproc --all)"      # Fetching total number of cpus
tnp=$((`expr $t - 3`)) # Maximum number of cpus will be utilized
val=$((`expr $tnp / 2`)) # Arithmatic operation on Max cpu count
np=$((`printf "%.*f\n" 0 $val`)) # Number of shared cpu

echo "test" $t $tnp $val $np
# Removing pre-exiting files if any

rm $out/K75_CNV_analysis.log

echo "First $input1 Second $input2" # Print inputs as a test

# Creating folders for the intermediate results
mkdir -p $out/$input1/FastQC_output
mkdir -p $out/$input1/trimmomatic_output
mkdir -p $out/$input1/FastQC_output_after-trimmomatic
mkdir -p $out/$input1/alignments_stats
mkdir -p $out/$input1/CNV_calling

eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs_gatk4                # Activating conda environment

# Fastqc and Trimmomatics run for Normal
while IFS= read -r line                 # While loop to add scripts for each lane in "parallel_fastqc_Normal" 
do                                      # "parallel_trimmomatic_Normal" and "parallel_after_trimmomatic_Normal"
        echo $line
        header=$(cat $raw_tumor/$line | head -n 1);  # Extracting header from fastq file. Will not work if fastq file
                                                        # is not gunzip file. Change command to cat

        id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');  # Extracting run ID from fastq file

        echo "header content" $header # printing header
        sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");  # Extracting file name from file 
	#nam=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-2 -d"_"); # list
        echo $id
        echo $sm
       # echo $nam

# QC checking

echo "Fastqc for $input1 Started" >> $out/K75_CNV_analysis.log   # Writing in log file
date >> $out/K75_CNV_analysis.log                                # Writing in log file

fastqc $raw_tumor/$line -t $tnp -o $out/$input1/FastQC_output/

echo "Fastqc for $input1 Completed" >> $out/K75_CNV_analysis.log         # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file

# Adaptor Trimming and cleaning

echo "Trimmomatic for $input1 Started" >> $out/K75_CNV_analysis.log      # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file

trimmomatic PE -threads $tnp -phred33 $raw_tumor/${sm}_All_R1* $raw_tumor/${sm}_All_R2* $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50
echo "Trimmomatic for $input1 Completed" >> $out/K75_CNV_analysis.log      # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file

# QC-rechecking after trimming

echo "Fastqc after trimmomatic for $input1 Started" >> $out/K75_CNV_analysis.log         # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file

fastqc -t $tnp $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input1/FastQC_output_after-trimmomatic/

done < "$input2"

echo "Fastqc after trimmomatic for $input1 Completed" >> $out/K75_CNV_analysis.log # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file

echo "BWA for $input1 Started" >> $out/K75_CNV_analysis.log      # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file

# BWA

while IFS= read -r line                 # While loop to add scripts for each lane in "parallel_BWA"
do

        echo $line
        header=$(cat $raw_tumor/$line | head -n 1);     # Extracting header from fastq file. Will not work if fastq file
                                                        # is not gunzip file. Change command to cat

        id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');   # Extracting run ID from fastq file
        echo "header content" $header   # printing header
        sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");   # Extracting file name from file
        #nam=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");  # list
        echo $id
        echo $sm
        #echo $nam

echo "BWA for $input1 Started" >> $out/K75_CNV_analysis.log      # Writing in log file
date >> $out/K75_CNV_analysis.log                                # Writing in log file

        if [ -z "$(ls -A  $out/$input1/trimmomatic_output/)" ]          # If statement
        then

        bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:K75\tSM:$sm\tPI:200" $hg38 $raw_tumor/${sm}_R1_001.fastq $raw_tumor/${sm}_R2_001.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam                                  # Will execute when 'trimmomatic_output' folder is empty

        else

        bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:K75\tSM:$sm\tPI:200" $hg38 $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam    #
Will execute when condition is false
   fi

done < "$input2"

echo "BWA for $input1 Completed" >> $out/K75_CNV_analysis.log      # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file


echo "Samtool Alignment Statistics for $input1 Started" >> $out/K75_CNV_analysis.log    # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file

sambamba flagstat -t $tnp $out/$input1/${sm}.sorted.bam  > $out/$input1/alignments_stats/${sm}.txt	# Getting alignment statistics

echo "Samtool Alignment Statistics for $input1 Completed" >> $out/K75_CNV_analysis.log    # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file


echo "CNV calling for $input1 Started" >> $out/K75_CNV_analysis.log    # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file

# Activating conda env for cnvkit
conda activate cnvkit

#cnvkit.py batch --method hybrid --normal Bam_files/${sm}.sorted.bam  --targets K75_finalbait.bed --fasta /mnt/database/GRCH38.P14/hg38.fa  --output-reference results_hybrid/my_reference.cnn  --output-dir results_hybrid/
cnvkit.py batch $out/$input1/${sm}.sorted.bam -r /home/genomics/bedfile/K75_CNV/my_reference.cnn -d $out/$input1/CNV_calling/

cnvkit.py segmetrics $out/$input1/CNV_calling/${sm}.sorted.cnr -s $out/$input1/CNV_calling/${sm}.sorted.cns --ci -o $out/$input1/CNV_calling/${sm}_segment.cns

cnvkit.py call $out/$input1/CNV_calling/${sm}_segment.cns --filter ci -m threshold -o $out/$input1/CNV_calling/${sm}_segment_call.cns

awk '{ if($6 >= 3) { print }} ' $out/$input1/CNV_calling/${sm}_segment_call.cns > $out/$input1/CNV_calling/${sm}_cnv_call_duplication

awk '{ if($6 <= 1) { print }} ' $out/$input1/CNV_calling/${sm}_segment_call.cns > $out/$input1/CNV_calling/${sm}_cnv_call_deletion

awk '!($4=="-")' $out/$input1/CNV_calling/${sm}_cnv_call_duplication  > $out/$input1/CNV_calling/${sm}_cnv_call_duplication_final

awk '!($4=="-")' $out/$input1/CNV_calling/${sm}_cnv_call_deletion  > $out/$input1/CNV_calling/${sm}_cnv_call_deletion_final

echo "CNV calling for $input1 Completed" >> $out/K75_CNV_analysis.log    # Writing in log file
date >> $out/K75_CNV_analysis.log                                        # Writing in log file
echo "##############################" >> $out/K75_CNV_analysis.log       # Writing in log file

