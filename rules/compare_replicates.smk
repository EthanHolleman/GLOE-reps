import itertools
# first need to figure out which samples are replicates of each other
# the modification column is the same for biological replicates
# could group these into tuples and expand the tuples then expand again?

GLOE_SAMPLES = pd.read_table(
    'samples/GLOE_samples.csv', sep=','
).set_index('Sample Name', drop=False)

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


rule compare_replicate_closeness:
    conda:
        '../envs/bedtools.yml'
    input:
        rep1='output/{sample_1}/mapped/{sample_1}.sorted.trim.depth.bed',
        rep2='output/{sample_2}/mapped/{sample_2}.sorted.trim.depth.bed'
    output:
        'output/compare_replicates/{sample_1}-{sample_2}.closest.bed'
    shell:'''
    bedtools closest -a {input.rep1} -b {input.rep2} > {output}
    '''

