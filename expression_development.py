from tool_box import exclude_cds
from tool_box import make_set_PTC_genes
from tool_box import compute_mean_std_error


# define excluded as a global variable
excluded = exclude_cds('briWS242.errors.txt')


def briggsae_expression_development(expression_development_file, CB_gene_IDs_file, PTC_file):
    '''
    (file, file, file) -> dict
    Returns a dictionary with gene expression data and PTC status for all genes
    '''

    CBIDs = open(CB_gene_IDs_file, 'r')
    expression_file = open(expression_development_file, 'r')
    
    # make a set of PTC genes
    from tool_box import make_set_PTC_genes
    PTC = make_set_PTC_genes(PTC_file)

    # make a dictionnary with cbr-CBGID as key and cbr-WBGID as value
    ID = {}
    for line in CBIDs:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            if line[-2].startswith('CBG'):
                ID[line[-2]] = line[1]

    # make a dictionnary for WBG gene containing average expression across development
    expression = {}
    expression_file.readline()
    expression_file.readline()
    for line in expression_file:
        line = line.rstrip()
        if line != '':
            line  = line.split()
            gene = line[0]
            total = 0
            for i in range(2, 12):
                total += float(line[i])
            average = total / 10
            expression[gene] = average
                

    # make a dictionnary of expression for PTC and non PTC genes using the CBG ID as key
    expression_PTC = {}
    for gene in ID:
        if ID[gene] in expression:
            expression_PTC[gene] = [expression[ID[gene]]]

    # delete genes that were removed from the analysis of SNP calling
    for gene in excluded:
        if gene in expression_PTC:
            del expression_PTC[gene]

    # add the PTC status
    for gene in expression_PTC:
        if gene in PTC:
            expression_PTC[gene].append('yes')
        else:
            expression_PTC[gene].append('no')

    expression_file.close()
    CBIDs.close()
    return expression_PTC
            


def ttest_PTC_expression_development(expression_development_file, CB_gene_IDs_file, PTC_file):
    '''
    (file, file, file) -> None
    Performs a t-test of independance between PTC genes and non-PTC genes for expression level
    Returns a tuple with the t-test statistics and the p-value assuming unequal variance
    '''

    expression_PTC = briggsae_expression_development(expression_development_file, CB_gene_IDs_file, PTC_file)

    from scipy import stats

    # make lists of expression level for PTC and non_PTC genes
    PTC_express = []
    non_PTC_express = []
    for gene in expression_PTC:
        if expression_PTC[gene][-1] == 'yes':
            PTC_express.append(expression_PTC[gene][0])
        else:
            non_PTC_express.append(expression_PTC[gene][0])
            
    
    expression = stats.ttest_ind(PTC_express, non_PTC_express, equal_var = False)
    
    expression_ttest = (abs(round(float(expression[0]), 4)), float(expression[1]))

    
    
    PTC_stats = compute_mean_std_error(PTC_express)
    non_PTC_stats = compute_mean_std_error(non_PTC_express)
   
    print('expression differences: t-test =\t{0}\tp-value =\t{1}'.format(expression_ttest[0], expression_ttest[1]))
    print('PTC: mean expression =\t%6.4f\tstandard error =\t%6.4f' % PTC_stats)
    print('non PTC: mean expression =\t%6.4f\tstandard error =\t%6.4f' % non_PTC_stats)
    




            
            
        
