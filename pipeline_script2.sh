#!/bin/sh

##colours 
RED='\033[0;31m'
BLUE='\033[0;34m'
LBLUE='\033[1;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

##initial description 
echo "${RED}T${BLUE}H${GREEN}I${YELLOW}S ${LBLUE}I${NC}S ${RED}T${BLUE}H${GREEN}E ${YELLOW}S${LBLUE}T${NC}A${RED}R${BLUE}T ${GREEN}O${YELLOW}F ${LBLUE}T${NC}H${RED}E ${BLUE}P${GREEN}I${YELLOW}P${LBLUE}E${NC}L${RED}I${BLUE}N${GREEN}E${NC}" 
echo " "
echo "This is a pipeline for returning a spoligotype from raw MinION data. This will consist of basecalling, producing a draft assembly, ${NC}creating an accurate consensus sequence and Spoligotype."
echo " "
##Listing all programmes, and storing their variable paths! INCLUDE FINAL / ON ALL PATHS
echo "${BLUE}The following programmes will be used:${NC}"
echo " "
echo "${LBLUE}Guppy${NC}"
guppy_path= #assumed to be in bashrc

echo "canu 1.8"
canu_path=/opt/canu-1.8/Linux-amd64/bin/ 

echo "Nanopolish"
nanopolish_path=/opt/nanopolish/

echo "Minimap2 (minimap2-2.17_x64-linux)"
minimap2_path=/opt/minimap2-2.17_x64-linux/

echo "${LBLUE}Samtools"
samtools_path= #assumed to be in bashrc

echo "${LBLUE}Blast"
blast_path= #assumed to be in bashrc

echo "Python2.7${NC}"
python_path= #should be in bashrc

echo "Spotyping-2.1 (uses Blast which is assumed to be in bashrc file)"
spotyping_path=/opt/SpoTyping-2.1/SpoTyping-v2.1-commandLine/ #PATH FOR SPOTYPING!

echo "${LBLUE}autocomplete.sh${NC}"
autocomplete_path=/home/phil/ #should be in bashrc

echo " "

echo "Currently the listed programmes have been assumed to have the path "/opt/", except for the ones in light blue which are assumed to be in your bashrc file or current working directory. If this is correct press ${GREEN}enter${NC}, if not please terminate this programme ${RED}(Ctrl^c)${NC} and input the paths manually into this script. Open in text editor or with the command touch pipeline_script.sh. Path variables for each should start at line 20."
read ENTER1

echo "All files and directories made will be placed in the current working directory. Are you running this script from the directory you intend to store the results in? Press ${GREEN}enter${NC} for yes or ${RED}(Ctrl^c)${NC} to end the programme, then please move into the correct directory"
read ENTER2
echo " "
##Now onto getting inputs#
###add path to autocomplete script###

echo "Please input the full path that contains your fastfile data files (the directory)" #give an example here?
read fast5
echo " "
#echo $fast5 #for debugging inputs

echo "All intermediate and final files will need their own name tag in order to allow you to easily search for them, please give a name tag to identify this sequencing data"
read nametag
echo " "
#echo $nametag #for debugging inputs

##The Script itself##
#####################
###Basecalling
echo "Now running guppy basecalling set to config to dna_r9.4.1_450bps_hac.cfg ; qscore_filtering on with min_qscore set to 7 ; cpu_threads_per_caller are 8 with 1 caller. Output will be in the guppy_nametag directory."
echo "${YELLOW}If a progress bar doesn't appear soon its likely that the path to your fast5 files is incorrect, if so press ctrl^c and do this again.${NC}"
$guppy_path\guppy_basecaller -r --input_path $fast5 --save_path ./guppy_$nametag --config dna_r9.4.1_450bps_hac.cfg --qscore_filtering --min_qscore 7 --cpu_threads_per_caller 8 --num_callers 1
echo "${GREEN}Basecalling step complete.${NC}"
echo " "
###Draft assembly production: Canu

###canu script
echo "Now running Canu. genomeSize is set to 5k (as we are working with Direct repeat region), and the draft assembly will be in a canu_nametag folder."

$canu_path\canu -p draft_assembly_$nametag -d ./canu genomeSize=5k -nanopore-raw ./guppy_$nametag/pass/*.fastq
echo "${GREEN}A draft assembly has been created.${NC}"
echo " "


####Nanopolish####
echo "Now a nanopolish consensus will be created using the canu draft"
echo " "

##index with nanopolish##
echo "First the fastq files are being index with the fast5 files"
$nanopolish_path\nanopolish index -d $fast5 ./guppy_$nametag/pass/*.fastq
echo "${GREEN}Fast5 indexing has occurred.${NC}"
echo " "

###minimap2 alignment
echo "Minimap2 alignment is now occuring..."
echo " "

##alignment with minimap2 to sam file#
echo "Alignment in progress..."
$minimap2_path\minimap2 -ax map-ont -t 8 ./canu/draft_assembly_$nametag.contigs.fasta ./guppy_$nametag/pass/*.fastq > alignment_$nametag.sam 
echo "${GREEN}sam file produced.${NC}"
echo " "

##conversion to bam file#
echo "conversion to bam format..."
$samtools_path\samtools view -S -b alignment_$nametag.sam > alignment_$nametag.bam
echo "${GREEN}bam file conversion succesfull.${NC}"
echo " "

##sorting of bam file#
echo "bam file is now being sorted"
$samtools_path\samtools sort -mof alignment_$nametag.bam alignment_$nametag.sorted
echo "${GREEN}Sorted.bam file now made.${NC}"
echo " "

##creating bai index#
echo "indexing the sorted.bam file..."
$samtools_path\samtools index alignment_$nametag.sorted.bam
echo "${GREEN}bai index file made.${NC}"
echo " "

###nanopolish step##

##polishing#
echo "Now creating a polished consensus sequence..."
$nanopolish_path\nanopolish variants --consensus -o polished_$nametag.vcf -r ./guppy_$nametag/pass/*.fastq -b alignment_$nametag.sorted.bam -g ./canu/draft_assembly_$nametag.contigs.fasta
echo " "
echo "${GREEN}vcf polished consensus sequence has been made!${NC}"

##making fasta format##
echo "Converting from vcf file to fasta file"
$nanopolish_path\nanopolish vcf2fasta --skip-checks -g ./canu/draft_assembly_$nametag.contigs.fasta polished_$nametag.vcf > polished_$nametag.fasta
echo "${GREEN}Polished consensus sequence is now in your current working directory and is called polished_$nametag.fasta.${NC}"
echo " "

#####future goal to add a MUMmer step to compare Canu to nanopolish


###Spoligotyping###
echo "Now performing spoligotyping..."   
$python_path\python2.7 $spotyping_path\SpoTyping.py --seq polished_$nametag.fasta -o $nametag.spoligotype
##displaying spoligotype
echo "${GREEN}Spoligotyping complete.${NC}"
echo "spolgiotype is..."
cat $nametag.spoligotype
echo "Spoligotyping logs and output file should be in your current working directory called $nametag.spoligotype.logs and $nametag.spoligotype respectively"
echo " "

#END STATEMENT
echo "${RED}F${BLUE}I${GREEN}N${YELLOW}I${LBLUE}S${NC}H${RED}E${BLUE}D${NC}"