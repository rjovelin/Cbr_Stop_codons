def make_final_PTC_table(stops_file, PTC_file):
    '''
    (file, file) -> file
    Copy the content of the inputfile into the outputfile without the information about GO, Interpro and frequencies,
    add the number of strains with coverage   
    '''

    inputfile = open(stops_file, 'r')
    outputfile = open(PTC_file, 'w')

    header1 = 'geneID' + '\t' + 'start' + '\t' + 'end' + '\t' + 'length_cds' + '\t' + 'Chr' + '\t'
    header2 = 'Chr_domain' + '\t' + 'N_PTC' + '\t' + 'PTC_pos' + '\t' + 'new_length' + '\t' + 'pos_mutation_in_codon' + '\t' +  'all_strains' + '\t' + 'all_strains_PTC' + '\t'
    header3 = 'temp_strains' + '\t' + 'temp_strains_PTC' + '\t' + 'trop_strains' + '\t' + 'trop_strains_PTC' + '\n'

    outputfile.write(header1)
    outputfile.write(header2)
    outputfile.write(header3)

    # make dictionnary without GO and Interpro
    data = {}
    for line in inputfile:
        if not line.startswith('gene'):
            line = line.rstrip()
            if line != '':
                line = line.split()
                data[line[0]] = []
                for i in range(0, 6):
                    data[line[0]].append(line[i])
                for i in range(-29, -16):
                    data[line[0]].append(line[i])

    # remove fields that are not included in headers
    positions = [16, 13, 10]
    for gene in data:
        for i in positions:
            del data[gene][i]
    
    # replace nb_ref with nb_strains_coverage (nb_strains_coverage = nb_ref + nb_snp)
    for gene in data:
        for i in range(10, 15, 2):
            if ';' in data[gene][i]:
                data[gene][i] = data[gene][i].split(';')
                data[gene][i+1] = data[gene][i+1].split(';')
                for y in range(len(data[gene][i])):
                    data[gene][i][y] = int(data[gene][i][y])
                    data[gene][i+1][y] = int(data[gene][i+1][y])
                    data[gene][i][y] += data[gene][i+1][y]
                    data[gene][i][y] = str(data[gene][i][y])
                    data[gene][i+1][y] = str(data[gene][i+1][y])
            else:
                data[gene][i] = str(int(data[gene][i]) + int(data[gene][i+1]))
    for gene in data:
        for i in range(10, 16):
            if type(data[gene][i]) == list:
                data[gene][i] = ';'.join(data[gene][i])

    # write the content of data to outputfile
    for gene in data:
        for i in range(0, len(data[gene]) -1):
            outputfile.write(data[gene][i] + '\t')
        outputfile.write(data[gene][-1] + '\n')
                
    outputfile.close()
    inputfile.close()

def make_syn_nonsyn_table(SNP_file, var_file):
    '''
    (file) -> file
    Save the content of the SNP_file, removing uncessary fiels
    Copy the content of the SNP_file into the var_file without the information about GO, Interpro and frequencies,
    add the number of strains with coverage 
    '''

    inputfile = open(SNP_file, 'r')
    outputfile = open(var_file, 'w')

    header1 = 'geneID' + '\t' + 'start' + '\t' + 'end' + '\t' + 'length_cds' + '\t' + 'Chr' + '\t'
    header2 = 'Chr_domain' + '\t' + 'N_SNP' + '\t' + 'SNP_pos' + '\t' + 'new_length' + '\t' + 'pos_mutation_in_codon' + '\t' +  'all_strains' + '\t' + 'all_strains_SNP' + '\t'
    header3 = 'temp_strains' + '\t' + 'temp_strains_SNP' + '\t' + 'trop_strains' + '\t' + 'trop_strains_SNP' + '\n'

    outputfile.write(header1)
    outputfile.write(header2)
    outputfile.write(header3)

    # make dictionnary without GO and Interpro
    data = {}
    for line in inputfile:
        if not line.startswith('gene'):
            line = line.rstrip()
            if line != '':
                line = line.split()
                data[line[0]] = []
                for i in range(0, 6):
                    data[line[0]].append(line[i])
                for i in range(-16, -3):
                    data[line[0]].append(line[i])

    # replace the new_length field
    for gene in data:
        data[gene][8] = '--'

    # remove fields that are not included in headers
    positions = [16, 13, 10]
    for gene in data:
        for i in positions:
            del data[gene][i]
    
    # replace nb_ref with nb_strains_coverage (nb_strains_coverage = nb_ref + nb_snp)
    for gene in data:
        for i in range(10, 15, 2):
            if ';' in data[gene][i]:
                data[gene][i] = data[gene][i].split(';')
                data[gene][i+1] = data[gene][i+1].split(';')
                for y in range(len(data[gene][i])):
                    data[gene][i][y] = int(data[gene][i][y])
                    data[gene][i+1][y] = int(data[gene][i+1][y])
                    data[gene][i][y] += data[gene][i+1][y]
                    data[gene][i][y] = str(data[gene][i][y])
                    data[gene][i+1][y] = str(data[gene][i+1][y])
            else:
                data[gene][i] = str(int(data[gene][i]) + int(data[gene][i+1]))
    for gene in data:
        for i in range(10, 16):
            if type(data[gene][i]) == list:
                data[gene][i] = ';'.join(data[gene][i])

    # write the content of data to outputfile
    for gene in data:
        for i in range(0, len(data[gene]) -1):
            outputfile.write(data[gene][i] + '\t')
        outputfile.write(data[gene][-1] + '\n')
                
    outputfile.close()
    inputfile.close()


def make_final_SCL_table(no_stops_file, stops_file, SCL_file):
    '''
    (file, file, file) -> file
    Copy the content of the inputfile into the outputfile without the information about GO, Interpro and frequencies,
    add the number of strains with coverage   
    '''

    inputfile = open(no_stops_file, 'r')
    outputfile = open(SCL_file, 'w')

    header1 = 'geneID' + '\t' + 'start' + '\t' + 'end' + '\t' + 'length_cds' + '\t' + 'Chr' + '\t'
    header2 = 'Chr_domain' + '\t' + 'N_PTC' + '\t' + 'PTC_pos' + '\t' + 'new_length' + '\t' + 'pos_mutation_in_codon' + '\t' +  'all_strains' + '\t' + 'all_strains_PTC' + '\t'
    header3 = 'temp_strains' + '\t' + 'temp_strains_PTC' + '\t' + 'trop_strains' + '\t' + 'trop_strains_PTC' + '\n'

    outputfile.write(header1)
    outputfile.write(header2)
    outputfile.write(header3)

    # make dictionnary without GO and Interpro
    data = {}
    for line in inputfile:
        if not line.startswith('gene'):
            line = line.rstrip()
            if line != '':
                line = line.split()
                data[line[0]] = []
                for i in range(0, 6):
                    data[line[0]].append(line[i])
                for i in range(-28, -15):
                    data[line[0]].append(line[i])

    # replace the new_length field, nb_stop, pos_stop
    for gene in data:
        data[gene][6] = '--'
        data[gene][7] = '--'
        data[gene][8] = '--'
        
    # remove fields that are not included in headers
    positions = [16, 13, 10]
    for gene in data:
        for i in positions:
            del data[gene][i]
    
    # replace nb_ref with nb_strains_coverage (nb_strains_coverage = nb_ref + nb_snp)
    for gene in data:
        for i in range(10, 15, 2):
            if ';' in data[gene][i]:
                data[gene][i] = data[gene][i].split(';')
                data[gene][i+1] = data[gene][i+1].split(';')
                for y in range(len(data[gene][i])):
                    data[gene][i][y] = int(data[gene][i][y])
                    data[gene][i+1][y] = int(data[gene][i+1][y])
                    data[gene][i][y] += data[gene][i+1][y]
                    data[gene][i][y] = str(data[gene][i][y])
                    data[gene][i+1][y] = str(data[gene][i+1][y])
            else:
                data[gene][i] = str(int(data[gene][i]) + int(data[gene][i+1]))
    for gene in data:
        for i in range(10, 16):
            if type(data[gene][i]) == list:
                data[gene][i] = ';'.join(data[gene][i])


    # remove genes that are PTC to keep only SCL genes
    from tool_box import make_set_PTC_genes
    PTC_genes = make_set_PTC_genes(stops_file)
    for gene in PTC_genes:
        if gene in data:
            del data[gene]

    # write the content of data to outputfile
    for gene in data:
        for i in range(0, len(data[gene]) -1):
            outputfile.write(data[gene][i] + '\t')
        outputfile.write(data[gene][-1] + '\n')
                
    outputfile.close()
    inputfile.close()




def remove_genes_with_PSC_frameshift(PTC_file, indel_file, new_PTC_file):
    '''
    (file, file, file) -> file
    Remove the genes from file PTC for which indels are causing the PTC and save the output to a new file
    '''

    # make a dictionnary of the PTC file
    PTC_snps = {}
    PTC = open(PTC_file, 'r')
    header = PTC.readline()
    for line in PTC:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            PTC_snps[gene] = line
    
    # make a dictionnary of the indel file
    indel_snps = {}
    indel = open(indel_file, 'r')
    indel.readline()
    for line in indel:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            indel_snps[gene] = line


    # remove genes for which the PTC is caused by an indel and genes that have frameshifts and PTC
    to_remove = []
    for gene in indel_snps:
        if indel_snps[gene][11] == 'yes':
            if indel_snps[gene][9] == 'indel':
                to_remove.append(gene)
            elif indel_snps[gene][9] == 'snp' and int(indel_snps[gene][4]) > 0:
                to_remove.append(gene)

    for gene in to_remove:
        if gene in PTC_snps:
            del PTC_snps[gene]

    # save the output to a new file
    new_PTC = open(new_PTC_file, 'w')
    new_PTC.write(header) 
    for gene in PTC_snps:
        for item in PTC_snps[gene][:-1]:
            new_PTC.write(item + '\t')
        new_PTC.write(PTC_snps[gene][-1] + '\n')
     
    PTC.close()
    indel.close()
    new_PTC.close()
    

def make_set_PTC_upstream_indel(PTC_file, indel_file, new_file):
    '''
    (file, file, file) -> file
    Make a table for the PTC genes that have frameshifts but that have a PTC SNP upstream the indel, keep only the 5' most SNP is multiple PTC SNPs are present
    and save the table to a new file
    '''

    # make a dictionnary of the PTC file
    PTC_snps = {}
    PTC = open(PTC_file, 'r')
    header = PTC.readline()
    for line in PTC:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            PTC_snps[gene] = line

    # make a dictionnary of the indel file
    indel_snps = {}
    indel = open(indel_file, 'r')
    indel.readline()
    for line in indel:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene = line[0]
            indel_snps[gene] = line

    # make a set of genes to keep
    to_keep = set()
    for gene in indel_snps:
        if indel_snps[gene][11] == 'yes' and indel_snps[gene][9] == 'snp':
            if int(indel_snps[gene][4]) > 0:
                to_keep.add(gene)

    # make a set of gene to remove:
    to_remove = set()
    for gene in PTC_snps:
        if gene not in to_keep:
            to_remove.add(gene)

    # remove genes from the PTC dict
    for gene in to_remove:
        del PTC_snps[gene]

    # keep only the most upstream SNP
    for gene in PTC_snps:
        if int(PTC_snps[gene][6]) > 1:
            PTC_snps[gene][7].split(';')
            PTC_snps[gene][7] = PTC_snps[gene][7][0]
            for i in range(8, 16):
                PTC_snps[gene][i].split(';')
                PTC_snps[gene][i] = PTC_snps[gene][i][0]

    # save the output to a new file
    new_PTC = open(new_file, 'w')
    new_PTC.write(header) 
    for gene in PTC_snps:
        for item in PTC_snps[gene][:-1]:
            new_PTC.write(item + '\t')
        new_PTC.write(PTC_snps[gene][-1] + '\n')
     
    PTC.close()
    indel.close()
    new_PTC.close()

    
    





























