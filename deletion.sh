OUTPUT_DIRECTORY=$1
OUTPUT_PREFIX=$2
PROGRAM_DIRECTORY=$3
    
##move
if [ ! -d ${OUTPUT_DIRECTORY}/files_for_report ];then
    mkdir ${OUTPUT_DIRECTORY}/files_for_report
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.LinkerFiltering_FastQ_PET.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.runningInformation.LinkerFiltering_FastQ_PET.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.tag_length_distribution.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.tag_length_distribution.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.linker_alignment_score_distribution.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.linker_alignment_score_distribution.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.linker_alignment_score_difference_distribution.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.linker_alignment_score_difference_distribution.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.linker_composition_distribution.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.linker_composition_distribution.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Mapping.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.runningInformation.Mapping.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.RemovingRedundancy.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.runningInformation.RemovingRedundancy.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Categories.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.runningInformation.Categories.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.clustering.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.runningInformation.clustering.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Peakcalling.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.runningInformation.Peakcalling.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.runningInformation.Report.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.runningInformation.Report.txt
    for x in 1_1 1_2 2_1 2_2
    do
        mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.mapping_statistics.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.${x}.mapping_statistics.txt
        mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.qc.dist.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.${x}.bedpe.qc.dist.txt
        mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.strand.dist.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.${x}.bedpe.strand.dist.txt
        mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.dist.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.${x}.bedpe.selected.dist.txt
        mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.intra-chrom.strand.dist.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.${x}.bedpe.selected.unique.intra-chrom.strand.dist.txt
    done
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.txt 
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.plusplus.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.plusplus.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.plusminus.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.plusminus.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.minusplus.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.minusplus.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.minusminus.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.bedpe.selected.intra-chrom.distance.minusminus.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.PET_count_distribution.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.PET_count_distribution.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.basic_statistics.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.basic_statistics.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.running_time.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.running_time.txt
    mv ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.summary.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.summary.txt
    
    ## copy
    cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt
    cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.FDRfiltered.txt ${OUTPUT_DIRECTORY}/files_for_report/${OUTPUT_PREFIX}.peak.FDRfiltered.txt
    
    ## remove
    for x in 1_1 1_2 2_1 2_2
    do
        rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R1.sai
        rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R1.sai.output.info.txt
        rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R2.sai
        rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.R2.sai.output.info.txt
        rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.sam.output.info.txt
        rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bam
        rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.${x}.bedpe.selected.sorted.txt
    done
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1.tagCount.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2.tagCount.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nTags.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.hypergeo.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.withpvalue.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.5K_5K
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.10K_10K
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.5K_5K
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.10K_10K
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nSpets.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.spetCounts.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.pois.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.withTagCounts
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.temp.txt
    rm ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.peak.aln
    rm ${PROGRAM_DIRECTORY}/../nTags.txt
    rm ${PROGRAM_DIRECTORY}/../data.txt
    rm ${PROGRAM_DIRECTORY}/../nSpets.txt
fi    
