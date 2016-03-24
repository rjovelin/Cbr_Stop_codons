def X_linked_constraints(divergence_file, genome_file):
    '''
    (file, file) -> None
    Print the results of t-tests comparing the mean divergence for X-linked genes and autosomes
    '''

    from PTC_chromosomes import gene_chromosome
    from scipy import stats
    from tool_box import compute_mean_std_error
    genes = gene_chromosome(genome_file)
    from PTC_genes_features import excluded
    
    # make a dictionnary with divergence data
    divergence = open(divergence_file, 'r')
    divergence_data = {}
    divergence.readline()
    for line in divergence:
        line = line.rstrip()
        if line != '':
            line  = line.split()
            divergence_data[line[0]] = []
            for i in range(1, 4):
                if line[i] == 'NA':
                    divergence_data[line[0]].append(line[i])
                else:
                    divergence_data[line[0]].append(float(line[i]))

    
    # create a new divergence dict by removing the letter at the end of the CBG gene
    constraints = {}
    for gene in divergence_data:
        if gene[-1].isalpha():
            name = gene[:-1]
            constraints[name] = divergence_data[gene]
        else:
            constraints[gene] = divergence_data[gene]
                    
    # delete genes with dS = 0 (omega = NA)
    impossible_ratio = []
    zero_dS = 0
    for gene in constraints:
        if constraints[gene][2] == 'NA':
            impossible_ratio.append(gene)
            zero_dS += 1
    for gene in impossible_ratio:
        del constraints[gene]

    # delete genes with dN/dS > 1
    large_omega = []
    omega_over_1 = 0
    for gene in constraints:
        if constraints[gene][2] > 1:
            large_omega.append(gene)
            omega_over_1 += 1
    for gene in large_omega:
        del constraints[gene]

    # delete genes in excluded
    for gene in excluded:
        if gene in constraints:
            del constraints[gene]

    # make a set of gene that have divergence and linkage
    rate = set()
    for gene in constraints:
        rate.add(gene)
    chromo = set()
    for gene in genes:
        chromo.add(gene)
    rate_chromo = rate.intersection(chromo)

    # add the chromosomal linkage to the genes in divergence data
    for gene in rate_chromo:
        if genes[gene] == 'X' or genes[gene] == 'X_random':
            constraints[gene].append('X')
        else:
            constraints[gene].append('auto')    

    # make a list of divergence values for X and autosomal genes
    dN_X = []
    dN_auto = []
    dS_X = []
    dS_auto = []
    omega_X = []
    omega_auto = []

    for gene in constraints:
        if gene in rate_chromo:
            if constraints[gene][-1] == 'X':
                dN_X.append(constraints[gene][0])
                dS_X.append(constraints[gene][1])
                omega_X.append(constraints[gene][2])
            else:
                dN_auto.append(constraints[gene][0])
                dS_auto.append(constraints[gene][1])
                omega_auto.append(constraints[gene][2])
              
    dN = stats.ttest_ind(dN_X, dN_auto, equal_var = False)
    dS = stats.ttest_ind(dS_X, dS_auto, equal_var = False)
    omega = stats.ttest_ind(omega_X, omega_auto, equal_var = False)

    dN_ttest = (abs(round(float(dN[0]), 4)), float(dN[1]))
    dS_ttest = (abs(round(float(dS[0]), 4)), float(dS[1]))
    omega_ttest = (abs(round(float(omega[0]), 4)), float(omega[1]))

    
    mean_dN_X = compute_mean_std_error(dN_X)
    mean_dN_auto = compute_mean_std_error(dN_auto)
    mean_dS_X = compute_mean_std_error(dS_X)
    mean_dS_auto = compute_mean_std_error(dS_auto)
    mean_omega_X = compute_mean_std_error(omega_X)
    mean_omega_auto = compute_mean_std_error(omega_auto)

    divergence.close()
    print('dN: t-test =\t {0}\t p-value =\t {1}'.format(dN_ttest[0], dN_ttest[1]))
    print('dS: t-test =\t {0}\t p-value =\t {1}'.format(dS_ttest[0], dS_ttest[1]))
    print('dN/dS: t-test =\t {0}\t p-value =\t {1}'.format(omega_ttest[0], omega_ttest[1]))
    print('X: mean dN =\t %6.4f,\t standard error =\t %6.4f' % mean_dN_X)
    print('auto: mean dN =\t %6.4f\t standard error =\t %6.4f' % mean_dN_auto)
    print('X: mean dS =\t %6.4f\t standard error =\t %6.4f' % mean_dS_X)
    print('auto: mean dS =\t %6.4f\t standard error =\t %6.4f' % mean_dS_auto)
    print('X: mean dN/dS =\t %6.4f\t standard error =\t %6.4f' % mean_omega_X)
    print('auto: mean dN/dS =\t %6.4f\t standard error =\t %6.4f' % mean_omega_auto)
    print('{0}\t genes with dN/dS > 1 were removed'.format(omega_over_1))
    print('{0}\t genes with dS = 0 were removed'.format(zero_dS))



def expression_all_genes(expression_file):
    '''
    (file) -> dict
    Return a dictionnary with the average male expression, and average female expression plus p-values of male and female comparisons
    '''
    from PTC_genes_features import excluded

    # make a dictionnary with expression data {gene: [avg_female, avg_male, PTC_status]}
    expression = open(expression_file, 'r')
    
    expression_data = {}
    expression.readline()
    for line in expression:
        line = line.rstrip()
        if line != '':
            line = line.split()
            expression_data[line[0]] = [float(line[4])]
            expression_data[line[0]].append(float(line[8]))
            expression_data[line[0]].append(float(line[9]))
            
    # delete genes that were removed from the analysis of SNP calling
    for gene in excluded:
        if gene in expression_data:
            del expression_data[gene]

    expression.close()
    return expression_data


def X_expression_difference(expression_file, genome_file):
    '''
    (file, file) -> None
    Print the results of t-tests comparing the mean expression for X-linked genes and autosomes for males and females
    '''

    from PTC_chromosomes import gene_chromosome
    from scipy import stats
    from tool_box import compute_mean_std_error
    from PTC_genes_features import excluded
    
    # make a dictionnary with gene information
    genes = gene_chromosome(genome_file)

    # make a dictionnary with expression data for all genes
    expression = expression_all_genes(expression_file)

    # make a set of gene that have expression and linkage
    expressed = set()
    for gene in expression:
        expressed.add(gene)
    chromo = set()
    for gene in genes:
        chromo.add(gene)
    expression_chromo = expressed.intersection(chromo)

    # add the chromosomal linkage to the genes in expression
    for gene in expression_chromo:
        if genes[gene] == 'X' or genes[gene] == 'X_random':
            expression[gene].append('X')
        else:
            expression[gene].append('auto')    

    # make a list of expression values for X and autosomal genes
    expression_female_X = []
    expression_female_auto = []
    expression_male_X = []
    expression_male_auto = []
    
    for gene in expression:
        if gene in expression_chromo:
            if expression[gene][-1] == 'X':
                expression_female_X.append(expression[gene][0])
                expression_male_X.append(expression[gene][1])
            else:
                expression_female_auto.append(expression[gene][0])
                expression_male_auto.append(expression[gene][1])
              
    female = stats.ttest_ind(expression_female_X, expression_female_auto, equal_var = False)
    male = stats.ttest_ind(expression_male_X, expression_male_auto, equal_var = False)
    
    female_ttest = (abs(round(float(female[0]), 4)), float(female[1]))
    male_ttest = (abs(round(float(female[0]), 4)), float(male[1]))
        
    mean_female_X = compute_mean_std_error(expression_female_X)
    mean_female_auto = compute_mean_std_error(expression_female_auto)
    mean_male_X = compute_mean_std_error(expression_male_X)
    mean_male_auto = compute_mean_std_error(expression_male_auto)
    
    print('female: t-test =\t {0}\t p-value =\t {1}'.format(female_ttest[0], female_ttest[1]))
    print('male: t-test =\t {0}\t p-value =\t {1}'.format(male_ttest[0], male_ttest[1]))
    print('female_X:\t %6.4f\t std_err\t %6.4f' % mean_female_X)
    print('female_auto:\t %6.4f\t std_err\t %6.4f' % mean_female_auto)
    print('male_X:\t %6.4f\t std_err\t %6.4f' % mean_male_X)
    print('male_auto:\t %6.4f\t std_err\t %6.4f' % mean_male_auto)


