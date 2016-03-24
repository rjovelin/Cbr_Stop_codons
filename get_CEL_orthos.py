def get_PTC_true_CEL_orthologs(CB_gene_IDs_file, CEL_CB_orthos_file, PTC_file, outputfile):
    '''
    (file, file, file, file) -> file
    Save the C. elegans 1 to 1 orthologs of the C. briggsae PTC genes into a new outputfile
    '''

    CBIDs = open(CB_gene_IDs_file, 'r')
    orthos = open(CEL_CB_orthos_file, 'r')
    CELorthoPTC = open(outputfile, 'w')

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

    # make a dictionnary of 1:1 cel-cbr orthologs with cbr-WBG as key and cel-WBG as value
    orthologs = {}
    orthos.readline()
    for line in orthos:
        line = line.rstrip()
        if line != '':
            line  = line.split()
            if int(line[-1]) == 1:
                orthologs[line[-4]] = line[0]

    # check if the PTC gene has ortholog
    # if so, write the cel-WBG to outputfile
    for gene in ID:
        if ID[gene] in orthologs and gene in PTC:
            CELorthoPTC.write(orthologs[ID[gene]] + '\n')

    CBIDs.close()
    orthos.close()
    CELorthoPTC.close()

    

            
            
        
