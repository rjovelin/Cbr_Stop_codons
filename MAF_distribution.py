def MAF_counts_distribution(SNP_file):
    '''
    (file) -> list
    Return the list of SNP counts in each MAF bin
    '''

    myfile = open(SNP_file, 'r')
    header = myfile.readline()

    # make a dictionnary with gene as key and list of MAF as value
    genes_MAF = {}
    for line in myfile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            genes_MAF[gene] = []
            if ';' in line[10]:
                strains = line[10].split(';')
                snps = line[11].split(';')
                for i in range(len(strains)):
                    MAF = (int(snps[i]) / int(strains[i])) * 100
                    genes_MAF[gene].append(MAF)
            else:
                genes_MAF[gene].append((int(line[11]) / int(line[10])) * 100)

    # make a list with all the MAF
    all_MAF = []
    for gene in genes_MAF:
        for MAF in genes_MAF[gene]:
            all_MAF.append(MAF)

    # make a list to update with MAF counts
    MAF_counts = [0] * 5

    # determine the index in the list where the MAF should be added
    # count the number of times a position appear within the window
    for MAF in all_MAF:
        which_range = int(MAF) // 20
        MAF_counts[which_range] += 1

    return MAF_counts


def compare_MAF_distribution(PTC_file, SNP_file):
    '''
    (file, file) -> tuple
    Performs a chi-square to compare the distribution of PTCs and other SNPs, and return a tuple with the chi-square, the p-value, and the degree of freedom
    '''
    from scipy import stats
    
    PTC = MAF_counts_distribution(PTC_file)
    SNP = MAF_counts_distribution(SNP_file)

    diff = stats.chi2_contingency([PTC, SNP])

    return(diff[0], diff[1], diff[2])


def histogram_MAF_distribution(SNP_file):
    '''
    (file) -> list
    Return a lit with frequencies of SNPs in each MAF bin to plot the histogram
    '''

    SNP = MAF_counts_distribution(SNP_file)
    total = 0
    for count in SNP:
        total += count

    MAF_freq = []
    for count in SNP:
        MAF_freq.append(count / total)

    return MAF_freq
            
            
def average_MAF_along_CDS(PTC_file, window):
    '''
    (file, int) -> list of list
    Return a list of lists, with each inner list containing the MAF value of PTC in bins of size window
    Precondition: only the most 5' PTC is considered
    '''

    from tool_box import compute_mean_std_error
   
    stops = open(PTC_file, 'r')
    header = stops.readline()

    # make a list of empty lists to be updated with the MAF values
    i = 100 // window
    PTC = []
    while i != 0:
        PTC.append([])
        i -= 1

    # add the MAF of the most 5' PTC to each bin
    for line in stops:
        line = line.rstrip()
        if line != '':
            line = line.split()
            position = int(line[8]) / int(line[3]) * 100
            if ';' in line[10]:
                strains = line[10].split(';')
                snps = line[11].split(';')
                MAF = (int(snps[0]) / int(strains[0])) * 100
            else:
                MAF = int(line[11]) / int(line[10]) * 100
            which_range = int(position) // window
            PTC[which_range].append(MAF)

    # compte MAF mean and std err for each bin
    average_MAF = []
    for i in range(len(PTC)):
        average = compute_mean_std_error(PTC[i])
        average_MAF.append(average)

    stops.close()
    return average_MAF

