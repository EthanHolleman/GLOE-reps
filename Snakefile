include: "rules/call_peaks.smk"
include: "rules/download_data.smk"
include: "rules/bowtie2.smk"
include: "rules/trimmomatic.smk"
include: "rules/process_alignment.smk"
include: "rules/compare_replicates.smk"
include: "rules/reorient_aligments.smk"
include: "rules/depth.smk"
include: "rules/remove_okazaki.smk"
include: "rules/compare_footloop.smk"
include: "rules/metaplot.smk"
include: "rules/rfd.smk"
include: "rules/fingerprints.smk"
include: "rules/general.smk"



wildcard_constraints:
    strand="fwd|rev",
    mode="indirect|direct",
    peak_call_method='normal|swapped',
    
    control="GSM[A-Za-z0-9]+",
    sample_a="GSM[A-Za-z0-9]+",
    sample_b="GSM[A-Za-z0-9]+",
    sample_c="GSM[A-Za-z0-9]+",
    sample_d="GSM[A-Za-z0-9]+"


import itertools

MODES = ['direct']
STRANDS = ['fwd', 'rev']
REGIONS = ['footloop', 'all']
PEAK_CALL_METHODS = ['normal', 'swapped']

# swapped peak calling is calling peaks with the control as the treatment

def group_replicates(samples):
    all_samples = []
    for _, sample in samples.iterrows():
        modification = sample['modification']
        replicates = samples.loc[samples['modification'] == modification]
        replicates = tuple([rep['Sample Name'] for _, rep in replicates.iterrows()])
        all_samples.append(replicates)
    return all_samples


replicate_groups = group_replicates(GLOE_SAMPLES)



multibamsummary = expand(
    'output/fingerprint/multibamsummary/{rep[0]}_{rep[1]}.all.summ.npz',
    rep=replicate_groups
)

paired_replicates_boxplot_output = expand(
    'output/compare_replicates/{rep[0]}-{rep[1]}/plots/closest/{mode}/{strand}/boxplot_{region}.png',
    rep=replicate_groups, allow_missing=True
)

treatment=list(GLOE_SAMPLES.loc[GLOE_SAMPLES['modification'] == 'siRNA LIG1']['Sample Name'])
control=list(GLOE_SAMPLES.loc[GLOE_SAMPLES['modification'] == 'siRNA control']['Sample Name'])
macs2_peaks = expand(
    'output/call_peaks/normal/{treatment}.vs.{control}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed',
    treatment=treatment, control=control,
    mode=MODES, strand=STRANDS, region=REGIONS,
    allow_missing=True
)


rfd = expand(
    'output/rfd/rfd_bw/{sample}.{mode}.all.1000.rfd.bw',
    sample=treatment, mode=MODES
)

all_peaks_stranded_average_normal = expand(
    'output/call_peaks/normal/{treatment_a}_{treatment_b}.vs.{control_a}_{control_b}_macs2/{mode}/all/macs2_peak_call_summits.strand.all.bed',
    zip, treatment_a=treatment, treatment_b=treatment, mode=MODES,
    control_a=control, control_b=control
)

all_peaks_stranded_average_normal = expand(
    'output/call_peaks/normal/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.all.macs2_peak_call_summits.strand.all.bed',
    zip, sample_a=treatment, sample_b=treatment[::-1], mode=MODES,
    sample_d=control[::-1], sample_c=control
)

all_peaks_stranded_average_swapped = expand(
    'output/call_peaks/swapped/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}_macs2.{mode}.all.macs2_peak_call_summits.strand.all.bed',
    zip, sample_a=control, sample_b=control[::-1], mode=MODES,
    sample_d=treatment[::-1], sample_c=treatment
)

intersection_hg19_swapped = expand(
    expand('output/metaplot/intersect/swapped/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.all.{strand}.intersect.bed', allow_missing=True, strand=STRANDS),
     zip, sample_a=control, sample_b=control[::-1], mode=MODES, strand=STRANDS,
    sample_d=treatment[::-1], sample_c=treatment
)


intersection_hg19_normal = expand(
    expand('output/metaplot/intersect/normal/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.all.{strand}.intersect.bed', allow_missing=True, strand=STRANDS),
    zip, sample_a=treatment, sample_b=treatment[::-1], mode=MODES, strand=STRANDS,
    sample_d=control[::-1], sample_c=control
)


metaplot_metabed_files_normal = expand(
    expand('output/metaplot/metabed/normal/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.all.{strand}.meta.bed', strand=STRANDS, allow_missing=True),
    zip, sample_a=treatment, sample_b=treatment[::-1], mode=MODES, strand=STRANDS,
    sample_d=control[::-1], sample_c=control
)

metaplot_metabed_files_swapped = expand(
    expand('output/metaplot/metabed/swapped/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.all.{strand}.meta.bed', allow_missing=True, strand=STRANDS),
    zip, sample_a=control, sample_b=control[::-1], mode=MODES, strand=STRANDS,
    sample_d=treatment[::-1], sample_c=treatment
)


# metaplots_normal = expand(
#     'output/metaplot/metabed/normal/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.all.{strand}.meta.bed',
#     zip, sample_a=treatment, sample_b=treatment[::-1], mode=MODES, strand=STRANDS,
#     sample_d=control[::-1], sample_c=control
# )

# metaplots_swapped = expand(
#     'output/metaplot/plots/swapped/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.all.all.meta.png',
#     zip, sample_a=control, sample_b=control[::-1], mode=MODES, strand=STRANDS,
#     sample_d=treatment[::-1], sample_c=treatment
# )


# all_peak_stranded_average_plots = expand(
#     'output/call_peaks/plots/{sample_c}_{sample_d}.{sample_a}_{sample_b}.{mode}.all.anno.all.png',
#     zip, sample_a=control, sample_b=control[::-1], mode=MODES,
#     sample_d=treatment[::-1], sample_c=treatment
# )

# plot_correlation = expand(
#     expand('output/fingerprint/plots/{sample_a}_{sample_b}.{sample_c}_{sample_d}.{mode}.all.{strand}.correlation.png',allow_missing=True, strand=STRANDS),
#     zip, sample_a=treatment, sample_b=treatment[::-1], mode=MODES, strand=STRANDS,
#     sample_d=control[::-1], sample_c=control
# )

plot_profiles = expand(
    expand('output/fingerprint/plots/{sample_a}_{sample_b}.vs.{sample_c}_{sample_d}.{mode}.all.{strand}.profile.png',allow_missing=True, strand=STRANDS),
    zip, sample_a=treatment, sample_b=treatment[::-1], mode=MODES, strand=STRANDS,
    sample_d=control[::-1], sample_c=control
)

plot_aligned_reads_profiles = expand(
    'output/fingerprints/plots/{sample}.{mode}.all.profile.{profile_type}.png',
    sample=(treatment + control), mode=MODES, profile_type='genes'
)

coverage_oriented_reads_genes = expand(
    'output/fingerprint/bed_coverage_genes/{sample}.{mode}.all.genes.coverage.hist.bed',
    sample=(treatment + control), mode=MODES)



print(coverage_oriented_reads_genes)




rule all:
    input:
        multibamsummary,
        coverage_oriented_reads_genes


