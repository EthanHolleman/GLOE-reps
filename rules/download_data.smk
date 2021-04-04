import pandas as pd

GLOE_SAMPLES = pd.read_table(
    'samples/GLOE_samples.csv', sep=','
).set_index('Sample Name', drop=False)


# Download GLOE-seq data and process into fastq

rule expand_gloe_samples:
    input:
        expand('rawdata/GLOE-seq/{sample_name}.sra', sample_name=GLOE_SAMPLES['Sample Name'])


rule download_all_gloe_samples:
    conda:
        '../envs/sra-toolkit.yml'
    params:
        sra_accession = lambda wildcards: GLOE_SAMPLES.loc[wildcards.sample_name]['Run'],
    output:
       temp('rawdata/GLOE-seq/{sample_name}.sra')
    shell:'''
    prefetch {params.sra_accession} --output-file {output}
    '''


rule dump_gloe_fastq:
    input:
        'rawdata/GLOE-seq/{sample}.sra'
    output:
        'rawdata/GLOE-seq/{sample}.fastq.gz'
    shell:'''
    fastq-dump -Z {input} | gzip > {output}
    '''

# Download primers




