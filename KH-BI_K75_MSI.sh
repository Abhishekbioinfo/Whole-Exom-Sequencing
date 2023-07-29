## Taking user inputs and setting resource paths
input1="$1"     # Taking input of sample folder name as argument
workdir=$(pwd)  # storing path of Working directory
raw="$workdir/$input1"   # Assigning path of sample
out="$workdir/Output_TSO500_MSI"       # Creating and Assigning path of output folder
db="/mnt/database"              # Assigning database folder path
hg38="$db/GRCH38.P14/hg38_homopolymer_microsatelittes.txt"   # Assigning path of human genome reference file
#msi="/home/genomics/software/msisensor2/./msisensor2 msi"              # Assigning software folder path
bed="/home/genomics/bedfile/K75/k75_complete_bed/K75_hg38_vardict.bed"

# Creating folders for the intermediate results
cd $workdir
mkdir -p $out/MSI_${input1}

#MSI
/home/genomics/software/msisensor2/./msisensor2 msi -d $hg38 -t ${input1}/bqsr/${input1}_bqsr.dedup.sorted.bam -e $bed -o $out/MSI_${input1}/${input1}_msi.txt -b 30
