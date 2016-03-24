def get_number_genes_alleles(PTC, position):
    '''
    (dict, int) -> tuple
    Return a tuple with the counts of all PTC genes, PTC alleles and non_singleton PTC genes, non_singleton alleles
    found in a given position of the list contained in dictionnary PTC
    '''

    all_genes_non_singleton = 0
    all_alleles_non_singleton = 0
    all_genes = 0
    all_alleles = 0
    for gene in PTC:
        SNPs = PTC[gene][position]
        SNPs = SNPs.split(';')
        null = SNPs.count('0')
        singleton = SNPs.count('1')
        all_alleles += (len(SNPs) - null)
        all_alleles_non_singleton += (len(SNPs) - null - singleton)
        if len(SNPs) != null:
            all_genes += 1
        if len(SNPs) != (null + singleton):
            all_genes_non_singleton += 1
    
    return (all_genes, all_alleles, all_genes_non_singleton, all_alleles_non_singleton)


def count_PTC_genes_alleles(PTC_file):
    '''
    (file) -> tuple
    Return the numbers of PTC genes and alleles in different groups of strains
    and paritioned by their frequency in each group
    '''

    myfile = open(PTC_file, 'r')
    header = myfile.readline()

    # make a dictionnary with gene as key and list of entries as value,
    # keeping the same order and content of entries as in the header

    PTC = {}
    for line in myfile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            PTC[line[0]] = line
    
    # count the total number of genes and alleles for non-singletons in all
    all_PTC = get_number_genes_alleles(PTC, 11)

    # count the number of PTC genes and alleles in the temperate strains
    temp_PTC = get_number_genes_alleles(PTC, 13)

    # count the number of non-singleton tropical genes and alleles
    trop_PTC = get_number_genes_alleles(PTC, 15)

    myfile.close()       
    return all_PTC, temp_PTC, trop_PTC
    


def print_PTC_count(PTC):
    '''
    (tuple) -> None
    Print the count of PTC genes and alleles in the different groups on screen
    '''
    print('all strains: the number of all PTC genes is:\t', PTC[0][0])
    print('all strains: the number of all PTC alleles is:\t', PTC[0][1])
    print('all strains: the number of non-singleton PTC genes is:\t', PTC[0][2])
    print('all strains: the number of non-singleton PTC alleles is:\t', PTC[0][3])
    print('temperate strains: the number of all PTC genes is:\t', PTC[1][0])
    print('temperate strains: the number of all PTC alleles is:\t', PTC[1][1])
    print('temperate strains: the number of non-singleton PTC genes is:\t', PTC[1][2])
    print('temperate strains: the number of non-singleton PTC alleles is:\t', PTC[1][3])
    print('tropical strains: the number of all PTC genes is:\t', PTC[2][0])
    print('tropical strains: the number of all PTC alleles is:\t', PTC[2][1])
    print('tropical strains: the number of non-singleton PTC genes is:\t', PTC[2][2])
    print('tropical strains: the number of non-singleton alleles is:\t', PTC[2][3])

