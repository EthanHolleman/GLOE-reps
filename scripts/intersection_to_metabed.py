import argparse
import os
import pandas as pd
import multiprocessing as mp
import math


DESCRIP = '''
Create a metabed (bed) file that describes the amount of signal in a given
percentage of a gene body. Created for data processsing for metaplots.
Should supply an intersected bed file containing genes (or other features)
intersected with a signal of interest and the bed file that contains the
features intersected (this would have been the file passed to bedtools intersect
-a parameter). 

The bedfile produced contains the total signal (sum of scores of intersections)
for each percentage region of the genes (features). 
'''

def parse_args():
    args = argparse.ArgumentParser(description=DESCRIP)
    args.add_argument('intersection_bed', metavar='I', 
    help='Path to gene-feature intersection file.')
    args.add_argument('gene_coordinates_bed', metavar='G', 
    help='Path to bedfile containing all genes features were intersected with.')
    args.add_argument('output_path', metavar='O', help='Path to output file.')
    args.add_argument('--threads', default=1, type=int, help='Number of threads.')
    args.add_argument('--n_bins', default=100, type=int, help='Number of bins per gene.')

    return args.parse_args()


def setup_mp_pool(number_threads):
    return mp.Pool(number_threads)


def read_gene_coords_as_df(gene_coords_bed):
    gene_coords_df = pd.read_csv(gene_coords_bed, sep='\t', header=None)
    gene_coords_df.columns = ['chromosome', 'start', 'end', 'name', 'score', 'strand']
    gene_coords_df[['name', 'id', 'name_extra']] = gene_coords_df['name'].str.split(';', expand=True)
    assert len(gene_coords_df.columns) == 8, 'Gene coordinate file has unexpected length'
    return gene_coords_df


def read_intersection_as_df(intersection_bed):
    intersect_df = pd.read_csv(intersection_bed, sep='\t', header=None)
    intersect_df[[3]] = intersect_df[3].str.split(';').str[0]
    # only take the gene name
    return intersect_df


def make_gene_intersections(intersect_df, gene_coords_df, pool, number_bins=100):

    args = [
            (
            intersect_df.loc[intersect_df[3] == row['name']],  # intersections for specific gene
            str(row['chromosome']),  
            int(row['start']),
            int(row['end']),
            str(row['name']),
            str(row['strand']),
            number_bins
        ) for _, row in gene_coords_df.iterrows()]
    print('made args')
    
    results = pool.starmap(bin_gene_intersections, args)
    return results
    
    
def bin_gene_intersections(gene_intersect_df, gene_chr, gene_start, gene_end, 
                           gene_name, gene_strand, number_bins=100):
    # smaller dataframe just for a specific gene
    bin_size = round(abs(gene_start - gene_end) / number_bins)
    bin_signal_df = pd.DataFrame(
        columns=('chr', 'bin_start', 'bin_end', 'name', 'signal', 'strand')
        )
    bins = [
        (gene_start + (bin_index * bin_size), 
        gene_start + ((bin_index + 1) * bin_size)
        ) 
        for bin_index in range(number_bins)]


    bins[len(bins)-1] = (bins[len(bins)-1][0], gene_end)
    bin_counter = 0
    for bin_start, bin_end in bins:
        intersecting_bin = gene_intersect_df.loc[
            (gene_intersect_df[7] >= bin_start) & (gene_intersect_df[8] <= bin_end)
            ]
        intersecting_bin_signal = sum(intersecting_bin[len(intersecting_bin.columns)-1])

        # the index of the bin is dependent on the strand, if fwd strand no
        # change. If reverse strand then the gene actually starts at the "end"
        # so need to flip the bin index
        if gene_strand == '-':
            bin_index = number_bins - (bin_counter + 1)
        else:
            bin_index = bin_counter

        bin_name = f'{gene_name}-{bin_index}-{bin_size}'
        bin_signal_df.loc[bin_counter] = (
            gene_chr, bin_start, bin_end, bin_name, intersecting_bin_signal, gene_strand)
        bin_counter += 1

    assert bin_signal_df['bin_start'][0] == gene_start, 'Gene start does not match bin start'
    assert bin_signal_df['bin_end'][len(bin_signal_df.index)-1] == gene_end
    assert len(bin_signal_df.index) == number_bins
    return bin_signal_df


def write_binned_genes_to_tsv(binned_genes, output_path):
    binned_genes_df = pd.concat(binned_genes)
    binned_genes_df.to_csv(output_path, header=False, index=False, sep='\t')


def main():

    args = parse_args()

    gene_coords_df = read_gene_coords_as_df(args.gene_coordinates_bed)
    print('read coords')
    intersect_df = read_intersection_as_df(args.intersection_bed)
    print('read intersect')
    pool = setup_mp_pool(number_threads=args.threads)
    binned_genes = make_gene_intersections(
        intersect_df, gene_coords_df, pool, number_bins=args.n_bins)
    write_binned_genes_to_tsv(binned_genes, args.output_path)


if __name__ == '__main__':
    main()