from syri.scripts.func import revcomp, readfasta

def test_ins_coord():
    '''
    Test that the matching bases for INS (before the inserted sequence) are same
    :return:
    '''
    refseq = readfasta("col.TAIR.fasta.filtered")
    qryseq = readfasta("Ler.chr.all.v1.0.fasta.filtered")
    with open("sv.txt", 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[0] != 'INS': continue
            try:
                if int(line[3]) < int(line[4]):
                    assert refseq[line[5]][int(line[1])-1] == qryseq[line[6]][int(line[3])-2]
                else:
                    assert refseq[line[5]][int(line[1])-1] == revcomp(qryseq[line[6]][int(line[3])+1-1])
            except AssertionError:
                print(line)
                return line
    return 0


def test_del_coord():
    '''
    Test that the matching bases for DEL (before the deleted sequence) are same
    :return:
    '''
    refseq = readfasta("col.TAIR.fasta.filtered")
    qryseq = readfasta("Ler.chr.all.v1.0.fasta.filtered")
    with open("sv.txt", 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[0] != 'DEL': continue
            try:
                # if int(line[3]) <= int(line[4]):
                assert refseq[line[5]][int(line[1])-2] == qryseq[line[6]][int(line[3])-1] or refseq[line[5]][int(line[1])-2] == revcomp(qryseq[line[6]][int(line[3])-1])
                # else:
                #     assert refseq[line[5]][int(line[1])-2] == revcomp(qryseq[line[6]][int(line[3])-1])
            except AssertionError:
                print(line)
                return line
    return 0



#TODO: test to check that the SV sequences are correct
