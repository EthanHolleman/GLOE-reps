include: "rules/download_data.smk"
include: "rules/bowtie2.smk"
include: "rules/trimmomatic.smk"
include: "rules/process_alignment.smk"
include: "rules/compare_replicates.smk"



replicate_output_path = 'output/compare_replicates/{}-{}.closest.bed'
replicate_comparison_targets = [
    replicate_output_path.format(key, value) for key, value in rep_dict.items()
    ]


rule all:
    input:
        # Download all GLOE-seq reads from SRA
        expand('rawdata/GLOE-seq/{sample_name}.sra', sample_name=GLOE_SAMPLES['Sample Name']),
        # Trimm GLOE reads
        expand(
            'output/{sample}/trimmomatic/{sample}.trimmed.fastq.gz', 
            sample=GLOE_SAMPLES['Sample Name']
        ),
        #replicate depth comparisons
        expand(replicate_comparison_targets)

