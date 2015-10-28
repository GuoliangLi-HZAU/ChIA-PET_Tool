PROGRAM_DIRECTORY='/home/local/ChIA-PET/program'
INPUT_INFO_FILE='MCF7.input.information.txt'
#INPUT_ANCHOR_FILE='/home/data/givenAnchor.cluster.txt' ### The path and file name of given anchors' clusters. If you don't have this file, please input 'null'. 
INPUT_ANCHOR_FILE='null'
OUTPUT_DIRECTORY='/home/data/results'
OUTPUT_PREFIX='MCF7'

SPECIES='1' ##1:human ; 2:mouse
CYTOBAND_DATA='hg19_cytoBandIdeo.txt'
CHROM_SIZE_INFO='hg19.chromSize.txt'
GENOME_LENGTH='3E9'  ### human
#GENOME_LENGTH='2.7E9'  ### mouse
GENOME_COVERAGE_RATIO='0.8' ### the proportion of the genome covered by the reads
GENOME_INDEX='/home/data/genome/hg19.genome.fa'

BWA='/root/Downloads/bwa-0.7.10/bwa'
NTHREADS='6' ### number of threads used in mapping reads to a reference genome
SAMTOOLS='/usr/local/bin/samtools'
BAM2BEDPE='/root/software/bedtools-2.17.0/bin/bamToBed -bedpe'

MAPPING_CUTOFF='20' ### cutoff of mapping quality score for filtering out low-quality or multiply-mapped reads

MERGE_DISTANCE='2'
SELF_LIGATION_CUFOFF='8000'
EXTENSION_LENGTH='500'
PVALUE_CUTOFF_PEAK='0.00001'
PVALUE_CUTOFF_INTERACTION='0.05'

PEAK_MODE='2' ### value: 1: the peak is an enriched region; 2: the peak is a local summit
MIN_DISTANCE_BETWEEN_PEAK='500' ### minimum distance between two different peaks
MIN_COVERAGE_FOR_PEAK='5' ### minimum coverage for peaks by extended reads


####### linker filtering
time=`date '+%s'` 
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.LinkerFiltering_FastQ_PET ${INPUT_INFO_FILE} ${OUTPUT_DIRECTORY} ${OUTPUT_PREFIX}
cat ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.1_1.R1.fastq ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.2_2.R1.fastq |wc -l |sed 's/ /\t/g' |  cut -f1 | awk '{print "Same-linker PETs after linker filtering\t"$1/4}'>> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt 
echo ${time} > ${OUTPUT_DIRECTORY}/time.txt
####### parameters
printf "parameters" > temp.txt
awk -v NTHREADS=${NTHREADS} -v GENOME_INDEX=${GENOME_INDEX} '{print "Number of threads\t"NTHREADS"\nGENOME_INDEX directory and prefix\t"GENOME_INDEX}' temp.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Mapping.txt
awk -v MAPPING_CUTOFF=${MAPPING_CUTOFF} -v MERGE_DISTANCE=${MERGE_DISTANCE} '{print "Mapping quality score cutoff\t"MAPPING_CUTOFF"\nMerge distance\t"MERGE_DISTANCE}' temp.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.RemovingRedundancy.txt
awk -v SELF_LIGATION_CUFOFF=${SELF_LIGATION_CUFOFF} '{print "Self-ligation cutoff\t"SELF_LIGATION_CUFOFF}' temp.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Categories.txt
awk -v EXTENSION_LENGTH=${EXTENSION_LENGTH} -v PVALUE_CUTOFF=${PVALUE_CUTOFF_INTERACTION} '{print "Tag extension length\t"EXTENSION_LENGTH"\nP-value cutoff\t"PVALUE_CUTOFF}' temp.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.clustering.txt
awk -v SELF_LIGATION_CUFOFF=${SELF_LIGATION_CUFOFF} -v MIN_DISTANCE_BETWEEN_PEAK=${MIN_DISTANCE_BETWEEN_PEAK} -v PEAK_MODE=${PEAK_MODE} -v MIN_COVERAGE_FOR_PEAK=${MIN_COVERAGE_FOR_PEAK} -v GENOME_LENGTH=${GENOME_LENGTH} -v GENOME_COVERAGE_RATIO=${GENOME_COVERAGE_RATIO} -v PVALUE_CUTOFF=${PVALUE_CUTOFF_PEAK} '{print "Self-ligation cutoff\t"SELF_LIGATION_CUFOFF"\nMin distance between peak\t"MIN_DISTANCE_BETWEEN_PEAK"\nPeak mode (1: region, 2: summit)\t"PEAK_MODE"\nMin coverage for peak\t"MIN_COVERAGE_FOR_PEAK"\nGenome length\t"GENOME_LENGTH"\nGenome coverage ratio\t"GENOME_COVERAGE_RATIO"\nP-value cutoff\t"PVALUE_CUTOFF}' temp.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Peakcalling.txt
awk -v CYTOBAND_DATA=${CYTOBAND_DATA} -v SPECIES=${SPECIES} '{print "Cytoband data\t"CYTOBAND_DATA"\nSpecies (1: human, 2: mouse)\t"SPECIES}' temp.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Report.txt
rm temp.txt
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt

####### mapping reads to a reference genome
for x in  1_1 1_2 2_1 2_2
do
    ${BWA} aln -n 2 -t ${NTHREADS} ${GENOME_INDEX} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R1.fastq 1> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R1.sai 2>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R1.sai.output.info.txt
    ${BWA} aln -n 2 -t ${NTHREADS} ${GENOME_INDEX} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R2.fastq 1> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R2.sai 2>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R2.sai.output.info.txt
    ${BWA} sampe -o 1 ${GENOME_INDEX} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R1.sai ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R2.sai ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R1.fastq ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R2.fastq 1>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.sam 2>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.sam.output.info.txt
    java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.util.MappingStatistics ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.sam ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x} ${MAPPING_CUTOFF}
    ${SAMTOOLS} view -Sb ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.temp.sam >  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bam
    ${BAM2BEDPE} -i ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bam > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe
    cut -f8    < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe | LANG=C sort -n | uniq -c > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.qc.dist.txt
    cut -f9,10 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe | LANG=C sort -n | uniq -c > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.strand.dist.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.sam
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.temp.sam
done
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt

####### purifying the mapped reads
for x in 1_1 1_2 2_1 2_2
do
    awk -v mapping_cutoff=${MAPPING_CUTOFF} '{if(($1!=".") && ($4!=".")){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t"mapping_cutoff"\t"$9"\t"$10}}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.txt
    LANG=C sort < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.sorted.txt
    uniq -c < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.sorted.txt | sed 's/^[ \t][ \t]*//g;s/[ \t][ \t]*/\t/g' | cut -f1 | LANG=C sort -n | uniq -c > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.dist.txt
    uniq < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.sorted.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.txt
## merging similar PETs according to chromosome and strand
    awk '{if($9=="+"){printf $1"\t"$3"\t"$9"\t"}else{printf $1"\t"$2"\t"$9"\t"}; if($10=="+"){print $4"\t"$6"\t"$10}else{print $4"\t"$5"\t"$10}}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.txt | LANG=C sort -k1,1 -k4,4 -k3,3 -k6,6 > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.pet.txt
    java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.MergeSimilarPETs2 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.pet.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.merged.pet.txt ${MERGE_DISTANCE}
    awk '{if($1==$4)print}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.txt | cut -f9,10 | LANG=C sort | uniq -c > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.intra-chrom.strand.dist.txt
done
cat ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.1_1.bedpe.selected.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.2_2.bedpe.selected.txt | wc -l |sed 's/ /\t/g' | awk '{print "Uniquely Mapped same-linker PETs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
cat ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.1_1.bedpe.selected.unique.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.2_2.bedpe.selected.unique.txt | wc -l |sed 's/ /\t/g' | awk '{print "Merging same same-linker PETs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt 

####### combining data
LANG=C sort -m ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.1_1.bedpe.selected.merged.pet.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.2_2.bedpe.selected.merged.pet.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.pet.txt
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.pet.txt |sed 's/ /\t/g' | awk '{print "Merging similar same-linker PETs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt

####### span distribution
awk '{if($1==$4){print $5-$2"\t"$3"\t"$6}}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.pet.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.txt
awk '{if(($2=="+")&&($3=="+"))print $1}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.plusplus.txt
awk '{if(($2=="+")&&($3=="-"))print $1}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.plusminus.txt
awk '{if(($2=="-")&&($3=="+"))print $1}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.minusplus.txt
awk '{if(($2=="-")&&($3=="-"))print $1}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.minusminus.txt

date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt

####### divide the PETs into different categories
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.PetClassification ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.pet.txt  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.opet ${SELF_LIGATION_CUFOFF}
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet | sed 's/ /\t/g' | awk '{print "Self-ligation PETs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet | sed 's/ /\t/g' | awk '{print "Inter-ligation PETs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.opet | sed 's/ /\t/g' | awk '{print "Other PETs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt

####### interaction calling
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.file.Pet2Cluster1 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt ${EXTENSION_LENGTH} ${PROGRAM_DIRECTORY}/${CHROM_SIZE_INFO}
LANG=C sort -k1,1 -k4,4  < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt >  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted
if [ ${INPUT_ANCHOR_FILE} != 'null' ];then
java -Xmx10G -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.PetClusterWithGivenAnchors ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted ${INPUT_ANCHOR_FILE} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster 1
else
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.PetCluster2 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster
fi
awk '{if($13>=2)print}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered
cut -f1-3,7-9,13-15 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.txt

###### calculation of p-values
cut -f1-3 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln
cut -f4-6 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet >>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln

cut -f1-3,13 <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1
cut -f7-9,13 <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2
for y in anchor1 anchor2
do
    java -Xmx10G -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.shortReads.TagCountInGivenRegions ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.${y} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.${y}.tagCount.txt 1 2
done

## generate the global tag count for global density
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln | sed 's/ /\t/g' |  cut -f1 >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nTags.txt
## calculate p-value
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nTags.txt nTags.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1.tagCount.txt  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2.tagCount.txt |  cut -f4,5,10  >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt data.txt
R --vanilla < ${PROGRAM_DIRECTORY}/hypergeometric.r
mv -f result.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.hypergeo.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.hypergeo.txt | cut -f1-3,7-9,13-15,19-24 > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.withpvalue.txt
awk -v pvalue_cutoff=${PVALUE_CUTOFF_INTERACTION} '{if($13<pvalue_cutoff+0.0)print}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.withpvalue.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt
awk '{print $7"\t"$8}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt |LANG=C sort |uniq -c | awk '{if($2>=10){a+=$1;b+=$1*$3;f+=$1;g+=$1*$3}else{c[$2]+=$1;d[$2]+=$1*$3;f+=$1;g+=$1*$3}}END{for(i in c){print c[i]"\t"i"\t"d[i]};print a"\t10\t"b;print f"\t11\t"g}' |LANG=C sort -k2,2n | awk 'BEGIN{print "PET counts\tNo. of clusters\tNo.intra-chrom clusters\tNo.inter-chrom clusters\tPercent of intra-chrom clusters"}{if($2==10){print ">="$2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else if($2==11){print "Total\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else{print $2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}}'> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.PET_count_distribution.txt

############################################################
####### peak calling
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.BindingSitesFromPETs  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak  ${EXTENSION_LENGTH}  ${SELF_LIGATION_CUFOFF}  ${MIN_COVERAGE_FOR_PEAK}  ${PEAK_MODE}  ${MIN_DISTANCE_BETWEEN_PEAK}

###### p-value calculation for peaks
# generate local tag counts for local density
awk '{center=int(($2+$3)/2); print $1"\t"$2"\t"$3"\t"center-5000"\t"center+5000"\t"$6}'   < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak | awk '{if($4 < 0){a=0}else{a=$4};if($2<0){b=0}else{b=$2}{print $1"\t"b"\t"$3"\t"a"\t"$5"\t"$6}}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K.temp
awk '{center=int(($2+$3)/2); print $1"\t"$2"\t"$3"\t"center-10000"\t"center+10000"\t"$6}' < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak | awk '{if($4 < 0){a=0}else{a=$4};if($2<0){b=0}else{b=$2}{print $1"\t"b"\t"$3"\t"a"\t"$5"\t"$6}}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K.temp
awk '{if(ARGIND==1){a[$1]=$2}else{b=a[$1];if($5>b){c=b}else{c=$5};if($3>b){d=b}else{d=$3}{print $1"\t"$4"\t"c"\t"$2"\t"d"\t"$6}}}' ${PROGRAM_DIRECTORY}/${CHROM_SIZE_INFO} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K.temp > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K
awk '{if(ARGIND==1){a[$1]=$2}else{b=a[$1];if($5>b){c=b}else{c=$5};if($3>b){d=b}else{d=$3}{print $1"\t"$4"\t"c"\t"$2"\t"d"\t"$6}}}' ${PROGRAM_DIRECTORY}/${CHROM_SIZE_INFO} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K.temp > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K
rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K.temp
rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K.temp

for y in 5K_5K 10K_10K
do
    java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.spetCountForPeaks ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.${y}  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.${y}
done

# generate the global spet count for global density
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.spet | sed 's/ /\t/g' |  cut -f1 >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nSpets.txt
# calculate p-value
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nSpets.txt nSpets.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.5K_5K  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.10K_10K |  cut -f5,6,13  >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.txt
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.txt data.txt
R --no-save --no-readline --args genomeLengthStr=${GENOME_LENGTH}  genomeCoverageRatioStr=${GENOME_COVERAGE_RATIO} extensionLengthStr=${EXTENSION_LENGTH} < ${PROGRAM_DIRECTORY}/pois.r
mv -f result.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt

# intra- and inter-chromosomal pet counts for peaks
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.TagCountForPeaks ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.withTagCounts
awk '{print $6"\t"$7"\t"$6+$7}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.withTagCounts >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.temp.txt

paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt | awk -v pvalue_cutoff=${PVALUE_CUTOFF_PEAK} '{if($8<pvalue_cutoff+0.0) print}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.tsv
awk '{print $1"\t"$2"\t"$3"\t.\t"0-100*log($8)/log(10)"\t.\t"$4"\t"$5}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.tsv >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.bed
awk '{print $1"\t"int(($2+$3)/2)"\t"$6}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.tsv  >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.aln

awk '{print $1":"$2"-"$3"\t"$6}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.compact
paste  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.compact  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.temp.txt | awk -v pvalue_cutoff=${PVALUE_CUTOFF_PEAK} '{if($4<pvalue_cutoff+0.0) print}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.long
cat ${PROGRAM_DIRECTORY}/peakHeader.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.long >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.xls

paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt |awk -v pvalue_cutoff=${PVALUE_CUTOFF_PEAK} '{if($8<pvalue_cutoff+0.0){print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.FDRfiltered.txt

wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.FDRfiltered.txt | sed 's/ /\t/g' | awk '{print "Peaks from self-ligation\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt | sed 's/ /\t/g' | awk '{print "Interacting pairs\t"$1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt

## running time 
awk '{printf $0"\t"}' ${OUTPUT_DIRECTORY}/time.txt |awk '{printf "Linker filtering\t";printf("%.2f",($2-$1)/60);printf "\nMapping\t";printf ("%.2f",($3-$2)/60); printf "\nRemoving redundancy\t"; printf ("%.2f",($4-$3)/60); printf "\nCategorization of PETs\t"; printf ("%.2f",($5-$4)/60); printf "\nPeak calling\t"; printf ("%.2f",($7-$6)/60); printf "\nClustering\t"; printf ("%.2f",($6-$5)/60); printf "\nTotal\t"; printf ("%.2f",($7-$1)/60); printf "\n"}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.running_time.txt
rm ${OUTPUT_DIRECTORY}/time.txt

## summary
sed -n '2p' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.linker_composition_distribution.txt | awk '{print "same-linker PETs/Total PETs\t"($1+$4)/$6}' > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk -F"\t" '{printf $2"\t"}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt | awk '{print "PETs after removing redundancy/Total PETs\t"$5/$1"\ninter-ligation PETs/PETs after removing redundancy\t"$7/$5}' >>  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk '{if($8==1){print}}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt | wc -l | awk '{printf $1"\t"}' >>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
awk '{if($8==1 && $9<1000000){print}}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt | wc -l | awk '{printf $1"\t"}' >>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt | awk '{print $1}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt 
sed -n '4p' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt |awk '{print "intra-chromosomal inter-ligation PETs/inter-ligation PETs\t"$1/$3"\nintra-chromosomal inter-ligation PETs within 1Mb/intra-chromosomal inter-ligation PETs\t"$2/$3}' >> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt
sed -i '4d' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt

sh ${PROGRAM_DIRECTORY}/../deletion.sh ${OUTPUT_DIRECTORY} ${OUTPUT_PREFIX} ${PROGRAM_DIRECTORY}

# Statistics Report
cp -Rf ${PROGRAM_DIRECTORY}/ChIA-PET_Tool_Report ${OUTPUT_DIRECTORY}/files_for_report
cp ${PROGRAM_DIRECTORY}/Rscript_and_genome_data/${CYTOBAND_DATA} ${OUTPUT_DIRECTORY}/files_for_report/${CYTOBAND_DATA}
Rscript --verbose ${PROGRAM_DIRECTORY}/Rscript_and_genome_data/ChIA-PET_Tool_Report.r ${PROGRAM_DIRECTORY} ${OUTPUT_DIRECTORY}/files_for_report ${OUTPUT_PREFIX} ${CYTOBAND_DATA} ${SPECIES}

if [ -d ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ChIA-PET_Tool_Report ];then
	    rm -rf ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ChIA-PET_Tool_Report
fi

mv ${OUTPUT_DIRECTORY}/files_for_report/ChIA-PET_Tool_Report/ ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ChIA-PET_Tool_Report

# Genomic Browser
## UCSC
awk -v OUTPUT_PREFIX=${OUTPUT_PREFIX} 'BEGIN{print "browser position chr1:9997500-10922500\nbrowser hide all\ntrack name=ChIAPET description=\"ChIA-PET "OUTPUT_PREFIX" human\""}$1==$4{print $1"\tChIAPET\t"OUTPUT_PREFIX"\t"$2"\t"$3"\t.\t.\t.\ttouch"NR"\n"$4"\tChIAPET\t"OUTPUT_PREFIX"\t"$5"\t"$6"\t.\t.\t.\ttouch"NR}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.toUCSC.gff
awk -v OUTPUT_PREFIX=${OUTPUT_PREFIX} 'BEGIN{print "browser position chr17:73626165-73912870\nbrowser hide all\nbrowser pack refGene encodeRegions\nbrowser full altGraph\n#300 base wide bar graph, autoScale is on by default == graphing\n#limits will dynamically change to always show full range of data\n#in viewing window, priority = 20 positions this as the second graph\n#Note, zero-relative, half-open coordinate system in use for bedGraph format\ntrack type=bedGraph name=\""OUTPUT_PREFIX" ChIAPET peak\" description=\""OUTPUT_PREFIX" ChIA-PET peaks\" visibility=full color=32,178,170 altColor=0,255,127 priority=20"}{print $1"\t"$2"\t"$3"\t"$4}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.FDRfiltered.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.toUCSC.bedGraph 
