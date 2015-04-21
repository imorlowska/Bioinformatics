# Sekwencje DNA - bioinformatyka lab 1 zadanie domowe
# Izabela Orlowska, python 3.4
# wywolanie: $ python Izabela.Orlowska.1.py <input file>
# przyklad: $ python Izabela.Orlowska.1.py hemoglobin.txt
def count(sequence):
        return {'A': sequence.count('A'), 'T': sequence.count('T'), 'C': sequence.count('C'), 'G': sequence.count('G')}

compl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
def complement(sequence):
	com = ""
	for x in sequence:
		com += compl[x]
	return com

gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }
	
def translate(sequence):
	return translate2(sequence, "", [])
	
def translate2(sequence, currstr="", all=[]):
	if (len(sequence) < 3):
		return all
	else:
		curr = gencode.get(sequence[0:3])
		if curr == '*':
			if currstr != "":
				all.append(currstr)
			return translate2(sequence[3:], "", all)
		else: 
			return translate2(sequence[3:], (currstr + curr), all)

def readFromFile(fileName):
	f = open(fileName, 'r')
	str = ""
	for line in f:
		str = str + line.rstrip()
	print ("\ncount:")
	print (count(str))
	print ("\ncomplement:")
	print (complement(str))
	print ("\ntranslate:")
	print (translate(str))
	print ("\ntranslate from 2nd:")
	print (translate(str[1:]))
	print ("\ntranslate from 3rd:")
	print (translate(str[2:]))
	print ("\nfinished")


if __name__ == "__main__":
    import sys
    readFromFile(sys.argv[1])
