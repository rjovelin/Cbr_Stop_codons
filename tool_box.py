def convert_fasta(fasta_file):
    '''
    (file) -> dict
    Take a file with fasta sequences and return a dictionnary with
    sequence ID as key and single string sequence as value
    '''
    
    sequences = {}
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line != '':
                if line.startswith('>'):
                    sequences[line[1: line.index(' ')]] = ""
                    gene = line[1: line.index(' ')]
                else:
                    sequences[gene] += line
    file.close()
    return sequences

def exclude_cds(excluded_file):
    '''
    (file) -> list
    Returns a set of gene that have excluded in the caling of PTC SNPs
    '''

    errors = open(excluded_file, 'r')

    excluded = set()

    for line in errors:
        gene = line[0: line.index('\t')]
        excluded.add(gene)

    return excluded


def make_set_PTC_genes(PTC_file):
    '''
    (file) -> set
    Return a set of genes that carry a PTC mutation
    '''
    PTC = open(PTC_file, 'r')
    
    # make a set of PTC genes
    PTC_genes = set()
    for line in PTC:
        if line.startswith('CBG'):
            line = line.rstrip().split()
            PTC_genes.add(line[0])

    PTC.close()
    return PTC_genes    


def compute_GC_content(S):
    '''
    (str) -> float
    Compute the GC content of a sequence S, not taking into account ambiguous nucleotides

    >>> compute_GC_content('atcgatcgatg')
    45.4545
    '''

    Sup = S.upper()
    valid = {'A', 'C', 'G', 'T'}
    GC = 0
    invalid = 0
    
    for base in Sup:
        if base in valid:
            if base == 'G':
                GC += 1
            elif base == 'C':
                GC += 1
        else:
            invalid += 1

    GC_content = (GC / (len(Sup) - invalid)) * 100

    return round(GC_content, 6) 

def cds_translate(cds):
    '''
    (str) -> str
    Translate a coding sequence into a protein sequence according to the standard genetic code

    >>> cds_translate('ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA')
    MAMAPRTEINSTRING*
    >>> cds_translate('ATGTACTAA')
    MY*
    '''

    genetic_code = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
                   'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
                   'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
                   'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
                   'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
                   'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                   'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                   'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                   'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
                   'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                   'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                   'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                   'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
                   'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                   'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                   'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
   

    CDS = cds.upper()
    protein = ''

    for i in range(0, len(CDS), 3):
        codon = CDS[i:i+3]
        if codon not in genetic_code:
            protein += 'X'
        else:
            protein += genetic_code[codon]

    return protein

def compute_mean_std_error(L):
    '''
    (list) -> tuple
    Returns a tuple containing the mean and the standard error of a collection of values in the list L
    Pre-condition: the values in L are floats and/or integers
    '''

    # verify the pre-condition
    for item in L:
        try:
            item + 1
        except:
            print('values in L need to be intergers and/or floats')

    import math

    # compute the mean
    total = 0
    for item in L:
        total += item
    mean = total/ len(L)

    # compute the stand error of the mean
    total_diff = 0
    for item in L:
        total_diff += (item - mean)**2
    std_dev = math.sqrt(total_diff / len(L))
    std_error = std_dev / math.sqrt(len(L))

    return (mean, std_error)

    
    
