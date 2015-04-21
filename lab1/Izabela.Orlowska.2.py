from collections import Counter
import re
import sys
# how to: python Izabela.Orlowska.2.py hemoglobin_amino_acid.txt hemoglobin.txt amino_acid_sequences.txt 

#czesc1:
#gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
#    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
#    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
#    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
#    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
#    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
#    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
#    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
#    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
#    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
#    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
#    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
#    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
#    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
#    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
#    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }
REGEX_ORG = "G{1,3}.[^P]{2}[RGI]*L+V"

# G=GGA, GGC, GGG, GGT = GG.
# P=CCA, CCC, CCG, CCT = CC.
# R=AGA, AGG, CGA, CGC, CGG, CGT = AG[AG] | CG.
# I=ATA, ATC, ATT = AT[^G]
# L=CTA, CTC, CTG, CTT, TTA, TTG = CT. | TT[AG]
# V=GTA, GTC, GTG, GTT = GT.
REGEX_TRAN = '(GG.){1,3}...([^C][^C].){2}(AG[AG]|CG.|GG.|AT[^G])*(CT.|TT[AG])+GT.'

#czesc2:
#group(0) - calosc, group(1) - fragment pasujacy, group(2) - Q, group(3) - ILM
AMINO_ACID_REGEX = "(([^.{0,3}]| Q).[^FHWY]([ILM])[^P][^FHILVWYP][DHFM][FMY])"

#try matching
def try_match(filename, regex, print_most_common=False):
    counter = Counter()
    i = 0
    for line in open(filename).readlines():
        match = re.search(regex, line)
        if match:
            print 'line ', i, 'matches pattern'
            if print_most_common:
                ilm = match.group(3)
                counter[ilm] += 1
        else:
            print 'line ', i, 'does not match pattern'
        i += 1
    if print_most_common:
        print 'Most common: ', counter.most_common()

#main
if __name__ == "__main__":
    print 'Czesc 1:'
    print '\t1.a'
    print 'regex=', REGEX_TRAN
    try_match(sys.argv[1], REGEX_TRAN) # hemoglobin_amino_acid.txt
    print '\t1.b'
    print 'regex=', REGEX_ORG
    try_match(sys.argv[2], REGEX_ORG) # hemoglobin.txt
    print
    print 'Czesc 2:'
    filename3 = sys.argv[3]
    try_match(filename3, AMINO_ACID_REGEX, print_most_common=True) # amino_acid_sequences.txt
