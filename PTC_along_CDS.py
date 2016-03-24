def PTC_positions_along_CDS(PTC_file, window):
    '''
    (file, int) -> list
    Return a list of PTC counts along the CDS in windows
    '''

    # make a list to store the relative position
    # only the 5' most PTC is taken into consideration
    stops = open(PTC_file, 'r')
    header = stops.readline()

    PTC = []
    total_snps = 0
    for line in stops:
        line = line.rstrip()
        if line != '':
            line = line.split()
            position = int(line[8]) / int(line[3]) * 100
            PTC.append(position)
            total_snps += 1

    # make a list of size window that contains only 0s:
    # each value in the list is the count of position for the range [0 - window[ etc
    range_counts = [0] * (100 // window)
        
    # determine the index in the list range_count where the truncation should be added
    # count the number of times a position appear within the window
    for position in PTC:
        which_range = int(position) // window
        range_counts[which_range] += 1

    stops.close()
    return range_counts

def is_PTC_positions_uniform(PTC_file, window):
    '''
    (file, int) -> None
    Performs a chi-square test of uniform distribution of the most 5' PTC along the CDS
    '''

    from scipy import stats
    range_counts = PTC_positions_along_CDS(PTC_file, window)
        
    # performs chi square test
    # note that the degree of freedom = k -1 - ddof
    # with k = number of observed frequencies
    # need to specify the parameter ddof to get the appropriate degree of freedom

    chisquare_test = stats.chisquare(range_counts)
    print('chisquare: %6.4f' % chisquare_test[0])
    print('p_value: {0}'.format(chisquare_test[1]))


def histogram_PTC_positions(PTC_file, window):
    '''
    (file, int) -> list
    Return a list of porportions of the 5' most PTC along the CDS in bins of size window
    to be use to make a histogram plot
    '''
    # count the number of PTC in each bin, taling into consideration only the most 5' PTC
    range_counts = PTC_positions_along_CDS(PTC_file, window)

    # count the number of such PTC
    stops = open(PTC_file, 'r')
    header = stops.readline()
    genes = set()
    for line in stops:
        line = line.rstrip()
        if line != '':
            line  = line.split()
            gene = line[0]
            genes.add(gene)
    total = len(genes)

    # calculate proportions:
    for i in range(len(range_counts)):
        range_counts[i] = range_counts[i] / total

    stops.close()
    return range_counts



def partition_CDS_length(PTC_file):
    '''
    (file) -> list
    Sort the genes in PTC_file according to their length
    and return a list of dictionnaries for each length quartile
    and the header of the file as last item
    '''

    # stote the PTC info in dictionnary
    stops = open(PTC_file, 'r')
    PTC = {}
    header = stops.readline()
    for line in stops:
        line = line.rstrip()
        if line != '':
            line = line.split()
            PTC[line[0]] = line

    # compute the length quartiles
    CDS_length = []
    for gene in PTC:
        CDS_length.append(int(PTC[gene][3]))

    quartiles = compute_quartiles(CDS_length)

    # sort genes into dictionnary of length quartile
    Q1, Q2, Q3, Q4 = {}, {}, {}, {}
    for gene in PTC:
        size = int(PTC[gene][3])
        if size < quartiles[0]:
            Q1[gene] = PTC[gene]
        elif size >= quartiles[0] and size < quartiles[1]:
            Q2[gene] = PTC[gene]
        elif size >= quartiles[1] and size < quartiles[2]:
            Q3[gene] = PTC[gene]
        else:
            Q4[gene] = PTC[gene]

    stops.close()
    return [Q1, Q2, Q3, Q4, header]


def save_quartiles_to_file(PTC_file):
    '''
    Sort the genes according to their CDS length, store them in dictionnaries
    of quartile length and save each dictionnary to a separate file
    '''

    quartiles = partition_CDS_length(PTC_file)
    header = quartiles[-1]

    for i in range(len(quartiles)-1):
        Q_file = open('PTC_polym_SNPS_only_Q' + str(i+1) +'.txt', 'w')
        Q_file.write(header)
        for gene in quartiles[i]:
            for item in quartiles[i][gene][:-1]:
                Q_file.write(item + '\t')
            Q_file.write(quartiles[i][gene][-1] + '\n')
        Q_file.close()


def compute_quartiles(L):
    '''
    (list) -> tuple
    Return the three quartile points of values contained in list L
    
    >>> compute_quartiles([1, 3, 3, 4, 5, 6, 6, 7, 8, 8,])
    (3, 5.5, 7)
    >>> compute_quartiles([3, 4, 4, 5, 6, 8, 8])
    (4, 5, 8)
    '''
  
    L.sort()
    if len(L) % 2 == 1:
        pos_median = int(len(L) / 2)
        median = L[pos_median]
        first_half = L[:pos_median]
        second_half = L[pos_median + 1:]

    elif len(L) % 2 == 0:
        pos1 = int(len(L) / 2)
        pos2 = pos1 - 1
        median = (L[pos1] + L[pos2]) / 2
        first_half = L[:pos1]
        second_half = L[pos1:]

    pos_Q1 = int(len(first_half) / 2)
    Q1 = first_half[pos_Q1]
    pos_Q3 = int(len(second_half) / 2)
    Q3 = second_half[pos_Q3] 

    return (Q1, median, Q3)



