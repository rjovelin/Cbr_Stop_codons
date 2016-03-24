# import tool_box functions
from tool_box import convert_fasta
from tool_box import exclude_cds
from tool_box import make_set_PTC_genes
from tool_box import compute_GC_content
from tool_box import cds_translate
from tool_box import compute_mean_std_error

# define excluded as a global variable
excluded = exclude_cds('briWS242.errors.txt')

def PTC_expression(expression_file, PTC_file):
    '''
    (file, file, file) -> dict
    
    Partition expression data in males and in females between PTC genes
    '''

    expression = open(expression_file, 'r')
    
    # make a set of PTC genes
    PTC_genes = make_set_PTC_genes(PTC_file)

    # make a dictionnary with expression data {gene: [avg_female, avg_male, PTC_status]}
    expression_data = {}
    expression.readline()
    for line in expression:
        line = line.rstrip()
        if line != '':
            line = line.split()
            expression_data[line[0]] = [float(line[4])]
            expression_data[line[0]].append(float(line[8]))
            if line[0] in PTC_genes:
                expression_data[line[0]].append('yes')
            else:
                expression_data[line[0]].append('no')

    # delete genes that were removed from the analysis of SNP calling
    for gene in excluded:
        if gene in expression_data:
            del expression_data[gene]

    expression.close()
    return expression_data


def ttest_PTC_expression_different(expression_file, PTC_file):
    '''
    (file, file, file) -> None
    Performs a t-test of independance between PTC genes and non-PTC genes for expression level
    Returns a tuple with the t-test statistics and the p-value assuming unequal variance
    '''

    expression_data = PTC_expression(expression_file, PTC_file)

    from scipy import stats

    # make lists of expression level in females for PTC and non_PTC genes
    female_PTC = []
    female_non_PTC = []
    for gene in expression_data:
        if expression_data[gene][-1] == 'yes':
            female_PTC.append(expression_data[gene][0])
        else:
            female_non_PTC.append(expression_data[gene][0])
            
    # make lists of expression level in males for PTC and non_PTC genes
    male_PTC = []
    male_non_PTC = []
    for gene in expression_data:
        if expression_data[gene][-1] == 'yes':
            male_PTC.append(expression_data[gene][1])
        else:
            male_non_PTC.append(expression_data[gene][1])

    female_expression = stats.ttest_ind(female_PTC, female_non_PTC, equal_var = False)
    male_expression = stats.ttest_ind(male_PTC, male_non_PTC, equal_var = False)

    female_ttest = (abs(round(float(female_expression[0]), 4)), float(female_expression[1]))
    male_ttest = (abs(round(float(male_expression[0]), 4)), float(male_expression[1]))

    female_stats_PTC = compute_mean_std_error(female_PTC)
    female_stats_non_PTC = compute_mean_std_error(female_non_PTC)
    male_stats_PTC = compute_mean_std_error(male_PTC)
    male_stats_non_PTC = compute_mean_std_error(male_non_PTC)

    print('expression in females: t-test =\t{0},\tp-value =\t{1}'.format(female_ttest[0], female_ttest[1]))
    print('expression in males: t-test =\t{0},\tp-value =\t{1}'.format(male_ttest[0], male_ttest[1]))
    print('PTC: mean male expression =\t%6.4f,\tstandard error =\t%6.4f' % male_stats_PTC)
    print('PTC: mean female expression =\t%6.4f,\tstandard error =\t%6.4f' % female_stats_PTC)
    print('non PTC: mean male expression =\t%6.4f,\tstandard error =\t%6.4f' % male_stats_non_PTC)
    print('non PTC: mean female expression =\t%6.4f,\tstandard error =\t%6.4f' % female_stats_non_PTC)
    
        
def PTC_expression_printer(expression, outputfile):
    '''
    (dict, filename) -> file
    Save a dictionnary of gene expression partitionned for PTC and nonPTC genes in a new outputfile
    '''

    myfile = open(outputfile, 'w')

    # write the header of the outputfile
    myfile.write('geneID' + '\t' + 'Aver_Female' + '\t' + 'Aver_Male' + '\t' + 'PTC_gene' + '\n')

    # sort the expression averages into a new file
    for gene in expression:
        myfile.write(gene + '\t' + str(expression[gene][0]) + '\t' + str(expression[gene][1]) +  '\t' + expression[gene][2] + '\n')
 
    myfile.close()


def PTC_sex_expression_bias(expression_file, PTC_file):
    '''
    (file, file, file) -> list
    Return a list of lists with the number of genes with significantly higher expression in females
    or in males for PTC genes and non-PTC genes
    '''

    expression = open(expression_file, 'r')

    # make a set of PTC genes
    PTC_genes = make_set_PTC_genes(PTC_file)

    # make a dictionnary with expression data
    expression_data = {}
    expression.readline()
    for line in expression:
        line = line.rstrip()
        if line != '':
            line = line.split()
            expression_data[line[0]] = line

    # delete genes that were removed from the analysis of SNP calling
    for gene in excluded:
        if gene in expression_data:
            del expression_data[gene]

    # make lists to update with gene counts
    # [higher_in male, higher_in_female]
    PTC_sex_expression = [0] * 2
    nonPTC_sex_expression = [0] * 2
    
    for gene in expression_data:
        male = float(expression_data[gene][8])
        female = float(expression_data[gene][4])
        if float(expression_data[gene][9]) < 0.05:
            if gene in PTC_genes:
                if female > male:
                    PTC_sex_expression[1] += 1
                elif female < male:
                    PTC_sex_expression[0] += 1
            else:
                if female > male:
                    nonPTC_sex_expression[1] += 1
                elif female < male:
                    nonPTC_sex_expression[0] += 1

    expression.close()
    return [PTC_sex_expression, nonPTC_sex_expression]


def PTC_sex_expression_differences(expression_file, PTC_file):
    '''
    (file, file, file) -> None
    Performs a  chi square test of uniform distribution of expression between PTC and nonPTC genes
    in males and females using the dictionary of expression sorted by sex and PTC status PTC_sex_expression
    and print the results to screen
    '''
    
    from scipy import stats

    sex_expression = PTC_sex_expression_bias(expression_file, PTC_file)
    
    # performs a chi_square test of contingency
    # unpack the tuple in variables    
    chi_2, p, df, expected = stats.chi2_contingency(sex_expression)

    print('chisquare:\t' +  '%6.4f' % chi_2)
    print('p_value:\t' + '{0}'.format(p))
    if sex_expression[0][1] != 0:
        print('PTC: male_bias/ female_bias:\t' + '{0}'.format(round(sex_expression[0][0] / sex_expression[0][1], 4)))
    if sex_expression[1][1] != 0:
        print('nonPTC: male_bias / female_bias:\t' + '{0}'.format(round(sex_expression[1][0] / sex_expression[1][1], 4)))





def no_reciprocal_blast_hits(blast_file):
    '''
    (file) -> list
    Make a list of genes that have no reciprocal blast hits to eliminate from the list of genes with divergence
    '''

    no_blast_genes = []

    no_blast = open(blast_file, 'r')
    header = no_blast.readline()
    for line in no_blast:
        line = line.rstrip()
        if line != '':
            line = line.split()
            ID = line[0]
            ID = ID.split('_')
            no_blast_genes.append(ID[-1])
    no_blast.close()
    return no_blast_genes



def PTC_divergence(divergence_file, PTC_file, blast_file):
    '''
    (file, file) -> (dict, int, int)
    Return a tuple with a dictionnary of divergence for PTC and non-PTC genes, along the with the number of genes excluded
    '''

    divergence = open(divergence_file, 'r')
    
    # make a set of PTC genes
    PTC_genes = make_set_PTC_genes(PTC_file)

    # make a list of genes that do not have reciprocal blast hits
    no_blast = no_reciprocal_blast_hits(blast_file)

    # make a dictionnary with divergence data
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
            if line[0] in PTC_genes:
                divergence_data[line[0]].append('yes')
            else:
                divergence_data[line[0]].append('no')

    # delete genes that were removed from the analysis of SNP calling
    for gene in excluded:
        if gene in divergence_data:
            del divergence_data[gene]

    # delete genes with no reciprocal blast
    for gene in no_blast:
        if gene in divergence_data:
            del divergence_data[gene]
          
    # delete genes with dS = 0 (omega = NA)
    impossible_ratio = []
    zero_dS = 0
    for gene in divergence_data:
        if divergence_data[gene][2] == 'NA':
            impossible_ratio.append(gene)
            zero_dS += 1
    for gene in impossible_ratio:
        del divergence_data[gene]

    # delete genes with dN/dS > 1
    large_omega = []
    omega_greater_1 = 0
    for gene in divergence_data:
        if divergence_data[gene][2] > 1:
            large_omega.append(gene)
            omega_greater_1 += 1
    for gene in large_omega:
        del divergence_data[gene]

    divergence.close()
    return divergence_data, zero_dS, omega_greater_1


def ttest_PTC_divergence_different(divergence_file, PTC_file, blast_file):
    '''
    (file, file) -> None
    Performs a t-test to compare the mean divergence between PTC and non-PTC genes, assuming unequal variance and print the results on screen
    '''

    divergence_data = PTC_divergence(divergence_file, PTC_file, blast_file)
    divergence = divergence_data[0]
    zero_dS = divergence_data[1]
    omega_greater_1 = divergence_data[2]

    from scipy import stats

    # make a list of divergence values for PTC and non-PTC genes
    dN_PTC = []
    dN_non_PTC = []

    dS_PTC = []
    dS_non_PTC = []

    omega_PTC = []
    omega_non_PTC = []

    for gene in divergence:
        if divergence[gene][-1] == 'yes':
            dN_PTC.append(divergence[gene][0])
            dS_PTC.append(divergence[gene][1])
            omega_PTC.append(divergence[gene][2])
        else:
            dN_non_PTC.append(divergence[gene][0])
            dS_non_PTC.append(divergence[gene][1])
            omega_non_PTC.append(divergence[gene][2])
              
    dN = stats.ttest_ind(dN_PTC, dN_non_PTC, equal_var = False)
    dS = stats.ttest_ind(dS_PTC, dS_non_PTC, equal_var = False)
    omega = stats.ttest_ind(omega_PTC, omega_non_PTC, equal_var = False)

    dN_ttest = (abs(round(float(dN[0]), 4)), float(dN[1]))
    dS_ttest = (abs(round(float(dS[0]), 4)), float(dS[1]))
    omega_ttest = (abs(round(float(omega[0]), 4)), float(omega[1]))

    # compute mean and standard error
    dN_stats_PTC = compute_mean_std_error(dN_PTC)
    dN_stats_non_PTC = compute_mean_std_error(dN_non_PTC)
    dS_stats_PTC = compute_mean_std_error(dS_PTC)
    dS_stats_non_PTC = compute_mean_std_error(dS_non_PTC)
    omega_stats_PTC = compute_mean_std_error(omega_PTC)
    omega_stats_non_PTC = compute_mean_std_error(omega_non_PTC)

    print('dN: t-test =\t {0},\t p-value =\t {1}'.format(dN_ttest[0], dN_ttest[1]))
    print('dS: t-test =\t {0},\t p-value =\t {1}'.format(dS_ttest[0], dS_ttest[1]))
    print('dN/dS: t-test =\t {0},\t p-value =\t {1}'.format(omega_ttest[0], omega_ttest[1]))
    print('PTC: mean dN =\t %6.4f,\t standard error =\t %6.4f' % dN_stats_PTC)
    print('nonPTC: mean dN =\t %6.4f,\t standard error =\t %6.4f' % dN_stats_non_PTC)
    print('PTC: mean dS =\t %6.4f,\t standard error =\t %6.4f' % dS_stats_PTC)
    print('nonPTC: mean dS =\t %6.4f,\t standard error =\t %6.4f' % dS_stats_non_PTC)
    print('PTC: mean dN/dS =\t %6.4f,\t standard error =\t %6.4f' % omega_stats_PTC)
    print('nonPTC: mean dN/dS =\t %6.4f,\t standard error =\t %6.4f' % omega_stats_non_PTC)
    print('{0}\t genes with dS = 0 were excluded'.format(zero_dS))
    print('{0}\t tgenes with dN/dS > 1 were excluded'.format(omega_greater_1))


def PTC_divergence_printer(divergence_data, outputfile):
    '''
    (dict, filename) -> file
    Save a dictionnary of divergence partitionned for PTC and nonPTC genes in a new outputfile
    '''

    myfile = open(outputfile, 'w')
    
    # write the header of the outputfile
    myfile.write('geneID' + '\t' + 'dN' + '\t' + 'dS' + '\t' + 'dN/dS' + '\n')

    # sort the divergence into a new file
    for gene in divergence_data:
        myfile.write(gene + '\t')
        for item in divergence_data[gene][:-1]:
            myfile.write(str(item) + '\t')
        myfile.write(divergence_data[gene][-1] + '\n')

    myfile.close()

def PTC_GC_content(fasta_file, PTC_file):
    '''
    (file, file) -> dict
    Returns a dictionnary with each briggsae gene as key and a list of value including the GC content of the coding sequence and the PTC status
    '''

    sequences = open(fasta_file, 'r')

    # make a set of PTC genes
    PTC_genes = make_set_PTC_genes(PTC_file)

    # convert the fasta_file into a dictionnary of fasta sequences
    cds = convert_fasta(fasta_file)

    # store the GC content of each gene in a dictionnary
    GC_content = {}

    for gene in cds:
        GC = compute_GC_content(cds[gene])
        GC_content[gene] = [GC]
        if gene in PTC_genes:
            GC_content[gene].append('yes')
        else:
            GC_content[gene].append('no')

    # delete genes that were removed from the analysis of SNP calling
    for gene in excluded:
        if gene in GC_content:
            del GC_content[gene]

    # delete genes that are in the SCL group
    for gene in SCL_genes:
        if gene in GC_content:
            del GC_content[gene]

    sequences.close()
    return GC_content  
    
def ttest_PTC_GC_content_different(fasta_file, PTC_file):
    '''
    (file, file, str) -> None
    Performs a t-test to compare the mean GC content between PTC and non-PTC genes, assuming unequal variance and print the results on screen
    '''

    GC_content = PTC_GC_content(fasta_file, PTC_file)

    from scipy import stats

    # make a list of GC values for PTC and non-PTC genes
    GC_PTC = []
    GC_non_PTC = []

    for gene in GC_content:
        if GC_content[gene][-1] == 'yes':
            GC_PTC.append(GC_content[gene][0])
        else:
            GC_non_PTC.append(GC_content[gene][0])
                          
    GC = stats.ttest_ind(GC_PTC, GC_non_PTC, equal_var = False)
    
    GC_ttest = (abs(round(float(GC[0]), 4)), float(GC[1]))

    GC_stats_PTC = compute_mean_std_error(GC_PTC)
    GC_stats_non_PTC = compute_mean_std_error(GC_non_PTC)
    
    print('GC: t-test = {0}, p-value = {1}'.format(GC_ttest[0], GC_ttest[1]))
    print('PTC: mean GC = %6.4f, standard error = %6.4f' % GC_stats_PTC)
    print('non PTC: mean GC = %6.4f, standard error = %6.4f' % GC_stats_non_PTC)

def PTC_GC_content_printer(GC_content, outputfile):
    '''
    (dict, filename) -> file
    Save a dictionnary of GC content partitionned for PTC and nonPTC genes in a new outputfile
    '''

    myfile = open(outputfile, 'w')
    
    # write the header of the outputfile
    myfile.write('geneID' + '\t' + 'GC' + '\n')

    # sort the divergence into a new file
    for gene in GC_content:
        myfile.write(gene + '\t')
        for item in GC_content[gene][:-1]:
            myfile.write(str(item) + '\t')
        myfile.write(GC_content[gene][-1] + '\n')

    myfile.close()


def codon_position_bias(PTC_file):
    '''
    (file) -> list
    Return a list with the frequency of PTC mutations at the 1st, 2nd and3rd codon position
    '''

    # convert the PTC_file into a dictionnary with gene as key and list of all entries as values
    PTC = {}
    stops = open(PTC_file, 'r')
    header = stops.readline()
    for line in stops:
        line = line.rstrip()
        if line != '':
            line = line.split()
            PTC[line[0]] = line

    # make a tuple to store the position counts
    codon_position = [0] * 3

    # go through each gene, then count the number of PTC at each codon position and update the list
    for gene in PTC:
        count1 = PTC[gene][9].count('1')
        count2 = PTC[gene][9].count('2')
        count3 = PTC[gene][9].count('3')
        codon_position[0] += count1
        codon_position[1] += count2
        codon_position[2] += count3

    # count the total number of snps
    total = 0
    for count in codon_position:
        total += count

    # compute the frequencies
    for i in range(0, 4):
        codon_position[i] = round((codon_position[i] / total), 4)
    
    return codon_position
        

