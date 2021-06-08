'''
This script will take in file from tgca with the following column names:
1. cance = cancer
2. gene = gene
3. sample = sample id
4. tumor_methylation = gene methylation in tumors
5. normal_methylation = gene methylation in normal cells
6. tumor_expression = gene expression in tumors
7. normal_expression = gene expression in normal cells

We want to generate:
1. a histogram of gene methylation in tumors and non tumors
2. a histogram of gene expression in tumors and non tumors
3. a t-test between tumor and normal samples for methylation and expression to determine if the gene is related to a cancer
4. have the numerical data present somewhere in figure
'''
##### import modules #####
import argparse
import pandas as pd
from scipy.stats import ttest_ind # to do the t-test
import matplotlib.pyplot as plt # to plot the figures


### create functions ###
def main(in_path, out_path):
    '''
    This function will read in the data from input and:
    1. Remove columns for which there's missing data for either cancer/normal For methylation and gene expression
    2. Calculate t-test pvalue for each
    3. Plot 2 figures (methylation/expression)--make it a histogram
        - the histogram should have the mean of each data plotted as well has have the p value listed
        - make sure to label the axis
    4. save the plots to the file specified in output
    '''
    cancer_data             = pd.read_csv(in_path, sep = '\t')
    
    ### get expression data first ###
    expression              = cancer_data.dropna(axis = 0, subset = ['tumor expression', 'normal expression']) # filter data
    expression_tumor_mean   = round(expression['tumor expression'].mean(), 4) # get mean for tumor
    expression_normal_mean  = round(expression['normal expression'].mean(), 4) # get mean for normal
    expression_pval         = ttest_ind(expression['tumor expression'], expression['normal expression']).pvalue # round numbers
    
    ### then do methylation data ###
    methylation             = cancer_data.dropna(axis = 0, subset = ['tumor methylation', 'normal methylation'])
    methylation_tumor_mean  = round(methylation['tumor methylation'].mean(), 4)
    methylation_normal_mean = round(methylation['normal methylation'].mean(), 4)
    methylation_pval        = ttest_ind(methylation['tumor methylation'], methylation['normal methylation']).pvalue
    
    ### plot the data ###
    fig, axs = plt.subplots(nrows = 2, ncols = 1) # initialize figure
    fig.suptitle('{}'.format(cancer_data.gene.unique()[0])) # Assumes only 1 gene, and it will be the figure title
    fig.subplots_adjust(hspace=0.4, wspace=0.4) # when I initially plotted, there wasn't enough space between subplots, so I add this

    # do expression first
    axs[0].hist(expression['tumor expression'], bins = 20, alpha = 0.5, label = 'tumor mean = {}'.format(expression_tumor_mean)) # plots tumor expression
    axs[0].hist(expression['normal expression'], bins = 20, alpha = 0.5, label = 'normal mean = {}'.format(expression_normal_mean)) # plots normal expression
    axs[0].legend(loc='upper left')
    axs[0].set_title('expression p-value = {}'.format(expression_pval))
    
    # then methylation
    axs[1].hist(methylation['tumor methylation'], bins = 20, alpha = 0.5, label = 'tumor mean = {}'.format(methylation_tumor_mean)) # plots tumor methylation
    axs[1].hist(methylation['normal methylation'], bins = 20, alpha = 0.5, label = 'normal mean = {}'.format(methylation_normal_mean)) # plots normal methylation
    axs[1].legend(loc='upper left')
    axs[1].set_title('methylation p-value = {}'.format(methylation_pval))
    
    plt.savefig(out_path)


if __name__ == '__main__':
    ### set arguments ###
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help = 'input data file')
    parser.add_argument('-o', '--output', help = 'output file png')
    argv = parser.parse_args() # we can call each argument using argv.input and argv.output
    main(argv.input, argv.output)
