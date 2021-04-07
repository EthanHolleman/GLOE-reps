# GLOE-seq paper focused on mapping okazaki fragments, but we are not
# actually interested in these for the purpose of R-loops. O-frags where called
# by comparing samples treated with siRNA knockout for ligase-1 against
# untreated controls. Calling peaks on these samples would then in theory
# reveal the location of breaks produced only in the absesne of ligase-1


# rule remove_called_peaks:
#     input:
#         peaks='output/call_peaks/{treatment}.vs.{control}_macs2/{mode}/{region}/{strand}/macs2_peak_call_summits.bed',
#         treatment='output/{treatment}/reorient_alignments/{mode}/{strand}/seperated.{region}.bed'
#     output:
#         'output/{treatment}/ok_free/{mode}/{strand}/ok_free.{region}.bed'
#     shell:'''
#     bedtools intersect -v -a {input.treatment} -b {input.peaks} > {output}
#     '''