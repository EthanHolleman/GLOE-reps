include: "rules/download_data.smk"
include: "rules/bowtie2.smk"
include: "rules/trimmomatic.smk"


rule all:
    input:
        # Download all GLOE-seq reads from SRA
        expand('rawdata/GLOE-seq/{sample_name}.sra', sample_name=GLOE_SAMPLES['Sample Name']),
        # Trimm GLOE reads
        expand(
            'output/{sample}/trimmomatic/{sample}.trimmed.fastq.gz', 
            sample=GLOE_SAMPLES['Sample Name']
        )

