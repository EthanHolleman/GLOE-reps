# Use the direct and indirect perl scripts provided in the GLOE-pipe pipline
# to correct orrientation or aligned reads. Different orrientation is due
# to GLOE-seq protocol


rule direct_mode:
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bed'
    output:
        sites='output/{sample}/direct/{sample}.sites.direct.trim.bed'
    params:
        index_dir='output/{sample}/direct'
    shell:'''
    module load perl
    mkdir -p {params.index_dir}
    perl scripts/direct_mode.pl {input} > {output.sites}
    '''


rule indirect_mode:
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bed'
    output:
        sites='output/{sample}/indirect/{sample}.sites.indirect.trim.bed'
    params:
        index_dir='output/{sample}/indirect'
    shell:"""
    module load perl
    mkdir -p {params.index_dir}
    perl scripts/indirect_mode.pl {input} > {output.sites}
    """


rule get_second_column:
    input:
        'output/{sample}/{mode}/{sample}.sites.{mode}.sorted.trim.bed'
    output:
        temp('output/{sample}/{mode}/{sample}.sites.{mode}.sorted.col2.trim.bed')
    shell:"""
    awk '($2 >= 0)' {input} > {output}
    """


rule perl_mode_big_awk:
    input:
        'output/{sample}/{mode}/{sample}.sites.{mode}.sorted.col2.trim.bed'
    output:
        'output/{sample}/{mode}/{sample}.{mode}.sorted.trim.bed'
    shell:"""
    awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6}}'  {input} > {output}
    """


rule seperate_forward_strand:
    input:
        'output/{sample}/{mode}/{sample}.{mode}.sorted.trim.bed'
    output:
        'output/{sample}/{mode}/{sample}.fwd.{mode}.sorted.trim.bed'
    
    shell:'''
    grep "+" {input} > {output}
    '''


rule seperate_reverse_strand:
    input:
        'output/{sample}/{mode}/{sample}.{mode}.sorted.trim.bed'
    output:
        'output/{sample}/{mode}/{sample}.rev.{mode}.sorted.trim.bed'
    shell:'''
    grep "-" {input} > {output}
    '''