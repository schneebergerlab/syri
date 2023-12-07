def align(coords, fix=False):
    '''
    Read alignments between chromosomes and create ancestral chromosome states.
    This is done by fusing the chromosomes to get larger 'ancestral chromosomes',
    while maximising syntenic regions.
    Optional: consider expanding it to multiple genomes
    :param coords: Alignments between the genomes (output of readCoords)
    :param fix: Fix the reference chromosomes. Only query chromosomes would be
    considered for fusing
    :return: Order of reference and query chromosome
    '''
    

    return
# END