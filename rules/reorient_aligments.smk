# Use the direct and indirect perl scripts provided in the GLOE-pipe pipline
# to correct orrientation or aligned reads. Different orrientation is due
# to GLOE-seq protocol


rule direct_mode:
    input:
        'output/{sample}/process_alignment/{sample}.sorted.trim.{region}.bed'
    output:
        sites=temp('output/{sample}/reorient_alignments/direct/{sample}.direct.trim.{region}.bed')
    params:
        index_dir='output/{sample}/direct'
    shell:'''
    module load perl
    mkdir -p {params.index_dir}
    perl scripts/direct_mode.pl {input} > {output.sites}
    '''


rule indirect_mode:
    input:
        'output/{sample}/process_alignment/{sample}.sorted.trim.{region}.bed'
    output:
        sites=temp('output/{sample}/reorient_alignments/indirect/{sample}.indirect.trim.{region}.bed')
    params:
        index_dir='output/{sample}/indirect'
    shell:"""
    module load perl
    mkdir -p {params.index_dir}
    perl scripts/indirect_mode.pl {input} > {output.sites}
    """


rule get_second_column:
    input:
        'output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.trim.{region}.bed'
    output:
        temp('output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.col2.trim.{region}.bed')
    shell:"""
    awk '($2 >= 0)' {input} > {output}
    """


rule perl_mode_big_awk:
    input:
        'output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.col2.trim.{region}.bed'
    output:
        'output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.bed'
    shell:"""
    awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6}}'  {input} > {output}
    """


rule seperate_forward_strand:
    input:
        'output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.bed'
    output:
        'output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.fwd.bed'
    
    shell:'''
    grep "+" {input} > {output}
    '''


rule seperate_reverse_strand:
    input:
        'output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.bed'
    output:
        'output/{sample}/reorient_alignments/{mode}/{sample}.{mode}.sorted.trim.{region}.rev.bed'
    shell:'''
    grep "-" {input} > {output}
    '''