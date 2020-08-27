
params.genome = "$HOME/ensamble/genome/Homo_sapiens.GRCh38.dna.chromosome.1.fa"
params.gtf = "$HOME/ensamble/gtf/Homo_sapiens.GRCh38.101.gtf"
params.reads = "$HOME/task/reads/SRR*_{1,2}_sub.fastq"
params.result = "$HOME/task/result"
//params.picard = "$HOME/task/picard"
params.picard = "/NGStools/"
params.bam="$HOME/task/result/SRR1/SRR1_Aligned.sortedByCoord.out.bam"

process STAR_index {

    input:
    file genome from file(params.genome)
    file gtf from file(params.gtf)
     
    output:
    file 'index' into index_channel

    script:       
    """
    STAR --runThreadN 6 \
        --runMode genomeGenerate \
        --genomeDir index \
        --genomeFastaFiles $genome \
        --sjdbGTFfile $gtf \
        --sjdbOverhang 99 \
        --genomeSAindexNbases 12
    """
}

Channel 
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .into { channel_allign; channel_qc; channel_metrics } 



process STAR_allign {
    publishDir params.result, mode:'copy'

    input:
    file index from index_channel
    set pair_id, file(reads) from channel_allign
    file output from file(params.result)

    output:
    set pair_id, "${pair_id}_Aligned.sortedByCoord.out.bam" into bam_channel
 
    script:       
    """
    STAR --genomeDir $index \
        --runThreadN 6 \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --quantMode GeneCounts \
        --outFileNamePrefix $output/$pair_id/${pair_id}_
    
    mv $output/$pair_id/${pair_id}_Aligned.sortedByCoord.out.bam .
    """
}

process fastqc {
    publishDir params.result, mode:'copy'

    input:
    set sample_id, file(reads) from channel_qc

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}

process picard_metrics {
    publishDir params.result, mode:'copy'

    input:
    path picard_path from params.picard
    file reference_genome from file(params.genome)
    set pair_id, file(bam) from bam_channel

    output:
    file("${pair_id}_picard_metrics") into picard_channel

    script:
    """
    mkdir ${pair_id}_picard_metrics
    java -jar $picard_path/picard.jar CollectAlignmentSummaryMetrics \
          R=$reference_genome \
          I=$bam \
          O=${pair_id}_picard_metrics/metrics.txt
    """  
} 

