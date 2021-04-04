include: "rules/download_data.smk"
include: "rules/bowtie2.smk"
include: "rules/trimmomatic.smk"
include: "rules/process_alignment.smk"
include: "rules/compare_replicates.smk"

import itertools


def group_replicates(samples):
    all_samples = []
    for _, sample in samples.iterrows():
        modification = sample['modification']
        replicates = samples.loc[GLOE_SAMPLES['modification'] == modification]
        replicates = tuple([rep['Sample Name'] for _, rep in replicates.iterrows()])
        all_samples.append(replicates)
    return all_samples


def make_replicate_key_dict(replicate_groups):
    permuted_groups = []
    for each_rep_group in replicate_groups:
        permute = itertools.permutations(each_rep_group)
        permuted_groups += permute
    return dict(permuted_groups)


group_reps = group_replicates(GLOE_SAMPLES)
rep_dict = make_replicate_key_dict(group_reps)
replicate_output_path = 'output/compare_replicates/{}-{}/closest.footloop.bed'
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

