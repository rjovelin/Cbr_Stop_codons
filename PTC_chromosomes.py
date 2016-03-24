from PTC_genes_features import excluded

def gene_chromosome(genome_file):
    '''
    (file) -> dict
    Return a dictionnary with the CBG gene as key and the chromosomal linkage as value
    '''
    # make a dictionnary to store the gene:chromo key:value pairs
    genome = open(genome_file, 'r')
    genes = {}
    for line in genome:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[8][line[8].index('CBG'):line[8].find(';', line[8].index('CBG'))]
            chromosome = line[0]
            genes[gene] = chromosome

    # make a list of genes with unknown chromosome
    unknown = []
    for gene in genes:
        if genes[gene] == 'un':
            unknown.append(gene)

    # remove all genes with unknown linkage
    for gene in unknown:
        del genes[gene]

    # remove all genes in excluded
    for gene in excluded:
        if gene in genes:
            del genes[gene]

    genome.close()
    return genes


def PTC_X_autosome(PTC_file, genome_file):
    '''
    (file, file) -> list
    Return a list of lists with the number of number of genes that are:
    PTC on the X, PTC on autosomes, non_PTC on the X, non_PTC on autosomes
    '''

    # make a set of PTC genes
    from tool_box import make_set_PTC_genes
    PTC_genes = make_set_PTC_genes(PTC_file)

    # make a dictionnary with the CBG names as keys and chromosome as value
    genes = gene_chromosome(genome_file)

    # count the number of genes in the different categories
    PTC_X = 0
    PTC_auto = 0
    nonPTC_X = 0
    nonPTC_auto = 0
    for gene in genes:
        if gene in PTC_genes:
            if genes[gene] == 'X' or genes[gene] == 'X_random':
                PTC_X +=1
            else:
                PTC_auto += 1
        else:
            if genes[gene] == 'X' or genes[gene] == 'X_random':
                    nonPTC_X += 1
            else:
                nonPTC_auto += 1

    return [[PTC_X, PTC_auto], [nonPTC_X, nonPTC_auto]]


def difference_PTC_X_auto(PTC_file, genome_file):
    '''
    (file, file) -> None
    Returns the result of a chisquare test (chi, p-value), testing if PTC and non_PTC genes are equally distributed on the X and autosomes
    '''

    from scipy import stats
    
    # get the number of genes in the different categories
    gene_count = PTC_X_autosome(PTC_file, genome_file)

    # performs a chi_square test of contingency
    # unpack the tuple in variables    
    chi_2, p, df, expected = stats.chi2_contingency(gene_count)

    print('chisquare:\t %6.4f' % chi_2)
    print('p_value:\t {0}'.format(p))
    if gene_count[0][1] != 0:
        print('PTC X/ PTC auto:\t {0}'.format(round(gene_count[0][0] / gene_count[0][1], 4)))
    if gene_count[1][1] != 0:
        print('nonPTC X/ nonPTC auto:\t {0}'.format(round(gene_count[1][0] / gene_count[1][1], 4)))


def PTC_chromosomal_domains(PTC_file, genome_file, chr_domains_file, chromo):
    '''
    (file, file) -> list of lists
    Returns a list of lists with first list of PTC counts in each chromosomal domain of a given chromosome
    and second list of expected number of counts correcting for the gene density in each domain,
    and a third list of gene count in each domain
    '''

    if chromo not in {'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX'}:
        chromo = input('please encter a valid chromosome name: ')

    domains = open(chr_domains_file, 'r')
    stops = open(PTC_file, 'r')
    genome = open(genome_file, 'r')

    # make a dictionnary to store the gene: (chromo, start) key:value pairs
    genes = {}
    for line in genome:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[8][line[8].index('CBG'):line[8].find(';', line[8].index('CBG'))]
            if '_' in line[0]:
                chromosome = 'chr' + line[0][:line[0].index('_')]
            else:
                chromosome = 'chr' + line[0]
            if line[6] == '+':
                start = int(line[3])
            elif line[6] == '-':
                start == int(line[4])
            genes[gene] = (chromosome, start)

    # make a list of genes with unknown chromosome
    unknown = []
    for gene in genes:
        if genes[gene][0] == 'chrun':
            unknown.append(gene)

    # remove all genes with unknown linkage
    for gene in unknown:
        del genes[gene]

    # remove genes that were deleted durint the SNP calling
    from PTC_genes_features import excluded
    for gene in excluded:
        if gene in genes:
            del genes[gene]

    # make a set of PTC genes
    from tool_box import make_set_PTC_genes
    PTC_genes = make_set_PTC_genes(PTC_file)

    # make a dictionnary of chromosomal domains {chromo:{'tip': N, 'arm': N, 'center': N}}
    chromo_domains = {}
    domains.readline()
    for line in domains:
        line = line.rstrip()
        if line != '':
            line = line.split()
            chromo_domains[line[0]] = {}
            chromo_domains[line[0]]['tip'] = [int(line[2]), int(line[-1])]
            chromo_domains[line[0]]['arm'] = [range(int(line[3]), int(line[4])), range(int(line[7]), int(line[8]))]
            chromo_domains[line[0]]['center'] = [range(int(line[5]), int(line[6]))]

    # count the number of PTC genes in the different chromosomal domains
    # [PTC_tip, PTC_arm, PTC_center]
    gene_count = [0] * 3
    for gene in genes:
        start = genes[gene][1]
        if genes[gene][0] == chromo:
            if gene in PTC_genes:
                if start < chromo_domains[chromo]['tip'][0] or start > chromo_domains[chromo]['tip'][1]:
                    gene_count[0] += 1
                elif start in chromo_domains[chromo]['arm'][0] or start in chromo_domains[chromo]['arm'][1]:
                    gene_count[1] += 1
                elif start in chromo_domains[chromo]['center'][0]:
                    gene_count[2] += 1

    # count expected number of PTC given the gene densoty in the different domains
    # expected = # genes in domain / # genes on chromo * # PTC genes on chromo

    total_PTC_chromo = 0
    for count in gene_count:
        total_PTC_chromo += count

    total_gene_chromo = 0
    for gene in genes:
        if genes[gene][0] == chromo:
            total_gene_chromo += 1

    total_tip = 0
    total_arm = 0
    total_center = 0
    for gene in genes:
        start = genes[gene][1]
        if genes[gene][0] == chromo:
            if start < chromo_domains[chromo]['tip'][0] or start > chromo_domains[chromo]['tip'][1]:
                total_tip += 1
            elif start in chromo_domains[chromo]['arm'][0] or start in chromo_domains[chromo]['arm'][1]:
                total_arm += 1
            elif start in chromo_domains[chromo]['center'][0]:
                total_center += 1

    expected = [total_tip, total_arm, total_center]
    for i in range(len(expected)):
        expected[i] = round(expected[i] / total_gene_chromo * total_PTC_chromo)

    totals = [total_tip, total_arm, total_center]

    genome.close()
    domains.close()
    stops.close()
    return [gene_count, expected, totals]


def difference_PTC_chromosomal_domains(PTC_file, genome_file, chr_domains_file, chromo):
    '''
    (file, file) -> None
    Returns the result of a chisquare test (chi, p-value), testing if PTC gene are distributed equally
    among the 3 chromosomal domains of a given chromosome, correcting for differences in gene density
    '''

    from scipy import stats
    
    # get the number of genes in the different categories
    # with their expected counts given gene densities
    gene_count = PTC_chromosomal_domains(PTC_file, genome_file, chr_domains_file, chromo)

    # performs ch_square test
    observed = gene_count[0]
    expected = gene_count[1]
    chisquare_test = stats.chisquare(observed, f_exp = expected)
    print('chisquare:\t %6.4f' % chisquare_test[0])
    print('p_value:\t {0}'.format(chisquare_test[1]))
    print('tip: PTC / total_genes:\t {0}'.format(round(gene_count[0][0] / gene_count[2][0], 4)))
    print('arm: PTC / total_genes:\t {0}'.format(round(gene_count[0][1] / gene_count[2][1], 4)))
    print('center: PTC / total_genes:\t {0}'.format(round(gene_count[0][2] / gene_count[2][2], 4)))
              



def SNP_positions(SNP_file):
    '''
    (file) -> dict
    Returns a dictionnary with a list of positions of each SNP for each chromsome
    '''

    number = lambda x: int(x)

    linkage = {}
    SNP = open(SNP_file, 'r')
    header = SNP.readline()
    for line in SNP:
        line = line.rstrip()
        if line != '':
            line = line.split()
            chromo = line[4]
            mutations = line[7].split(';')
            if chromo in linkage:
                linkage[chromo].extend(mutations)
            else:
                linkage[chromo] = mutations

    for chromo in linkage:
        linkage[chromo] = list(map(number, linkage[chromo]))
    
    SNP.close()
    return linkage
    

def SNP_counts_chromosomal_domains(SNP_file, chr_domains_file, chromo):
    '''
    (file, file, str) -> list
    Returns a list of SNP counts located in the tips, arms and center of a given chromosome
    '''
    
    if chromo not in {'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX'}:
        chromo = input('please enter a valid chromosome name: ')

    domains = open(chr_domains_file, 'r')

    # get the SNP positions for each chromosome
    SNP_pos = SNP_positions(SNP_file)

    # make a dictionnary of chromosomal domains {chromo:{'tip': N, 'arm': N, 'center': N}}
    chromo_domains = {}
    domains.readline()
    for line in domains:
        line = line.rstrip()
        if line != '':
            line = line.split()
            chromo_domains[line[0]] = {}
            chromo_domains[line[0]]['tip'] = [int(line[2]), int(line[-1])]
            chromo_domains[line[0]]['arm'] = [range(int(line[3]), int(line[4])), range(int(line[7]), int(line[8]))]
            chromo_domains[line[0]]['center'] = [range(int(line[5]), int(line[6]))]

    # count the number of SNPs in the different chromosomal domains
    # [tip_counts, arm_counts, center_counts]
    SNP_count = [0] * 3
    for allele_position in SNP_pos[chromo]:
        if allele_position < chromo_domains[chromo]['tip'][0] or allele_position > chromo_domains[chromo]['tip'][1]:
            SNP_count[0] += 1
        elif allele_position in chromo_domains[chromo]['arm'][0] or allele_position in chromo_domains[chromo]['arm'][1]:
            SNP_count[1] += 1
        elif allele_position in chromo_domains[chromo]['center'][0]:
            SNP_count[2] += 1

    domains.close()
    return SNP_count


def SNP_counts_difference_arms_center(PTC_file, SNP_file, chr_domains_file, chromo):
    '''
    (file, file, file, str) -> None
    Returns the result of a chisquare test (chi, p-value), testing if the proportions of PTC SNPs in ther arm and centers of a given chromosome
    is significantly different than the proportion of another class of SNPs in the same domains on the same chromsome
    '''

    from scipy import stats
    
    # get the SNP counts 
    PTC_count = SNP_counts_chromosomal_domains(PTC_file, chr_domains_file, chromo)
    SNP_count = SNP_counts_chromosomal_domains(SNP_file, chr_domains_file, chromo)
    
    # performs ch_square test
    PTC = [PTC_count[1], PTC_count[2]]
    SNP = [SNP_count[1], SNP_count[2]]
    SNP_counts = [PTC, SNP]
    
    chi2 = stats.chi2_contingency(SNP_counts)
    print('chisquare:\t%6.4f' % chi2[0])
    print('p_value:\t{0}'.format(chi2[1]))
    print('arm: PTC / SYN:\t{0}'.format(round(PTC_count[1] / SNP_count[1], 4)))
    print('center: PTC / SYN:\t{0}'.format(round(PTC_count[2] / SNP_count[2], 4)))


