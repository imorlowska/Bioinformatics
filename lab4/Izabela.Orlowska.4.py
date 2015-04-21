from Bio import Entrez, SeqIO
from Bio.SubsMat import MatrixInfo as mi

#macierz substytucji
values = mi.blosum62
#kara za przerwanie
gap_penalty = -7

#zwraca wartosc dopasowania z blosum62, np A i B daja -2
def lookup(a, b):
	if (a, b) in values:
		return values[(a, b)]
	else:
		return values[(b, a)]

#zwraca tablice wypelniana zerami oraz i*d w 0 kolumnie i 0 wierszu
#przygotowanie tablicy do pozniejszych obliczen
def init(str1, str2):
	tab = [[0] * (len(str1) + 1) for i in range(len(str2) + 1)]
	for i in range(1, len(str1) + 1):
		tab[0][i] = i * gap_penalty
	for i in range(1, len(str2) + 1):
		tab[i][0] = i * gap_penalty
	return tab

#generuje pusta tablice (zawierajaca puste listy), tu beda trzymane kierunki
def gen_array(n, m):
    arr = []
    for x in range(0, m):
        tmp = []
        for y in range(0, n):
            tmp.append([])
        arr.append(tmp)
    return arr

#glowna czesc algorytmu, dopasowuje wartosc wg schematu:
# tab[i][j] = max(
#				tab[i-1][j-1] + values(a,b), // pasuja
#				tab[i-1][j] + gap_penalty) // delecja
#				tab[i][j-1] + gap_penalty)) // insercja
# przy okazji w tablicy move usupelnia informacje, ktory ruch byl najkorzystniejszy.
# poniewaz chcemy wszystkie optymalne wyniki, dlatego musimy sprawdzic czy przypadkiem
# wszystkie przypadki nie sa rowne max
def fill(str1, str2, tab, move):
	for j in range(1, len(str1) + 1):
		for i in range(1, len(str2) + 1):
			match = tab[i-1][j-1] + lookup(str1[j-1], str2[i-1])
			delete = tab[i-1][j] + gap_penalty
			insert = tab[i][j-1] + gap_penalty
			maks = max(match, delete, insert)
			if maks == match:
				move[i][j].append('diag')
			if maks == delete:
				move[i][j].append('up')
			if maks == insert:
				move[i][j].append('left')
			tab[i][j] = maks

# backtracking na podstawie tablicy move, idziemy od konca wyrazow, zbieramy je w akumulatorach
# acc1 i acc2, gdy dochodzimy do 'ramki' od gory lub lewej - pustych miejsc w move - to wypisujemy
# wynik. tutaj tez nieladne if-y, ale mozemy miec wiele optymalnych rozwiazan (np. dla "AA" i "AAAA")
def backtrack2(str1, str2, move, i, j, acc1, acc2):
	x = move[i][j]
	if x.count('diag') > 0: #match
		backtrack2(str1, str2, move, i-1, j-1, str1[i-1]+acc1, str2[j-1]+acc2)
	if x.count('up') > 0: #delete
		backtrack2(str1, str2, move, i-1, j, str1[i-1]+acc1, "-"+acc2)
	if x.count('left') > 0: #insert
		backtrack2(str1, str2, move, i, j-1, "-"+acc1, str2[j-1]+acc2)
	if len(x) == 0:
		print 'best match:'
		print "1: " + acc1
		print "2: " + acc2

#opakowuje backtrack2 dla wygody testowania
def backtrack(str1, str2, move):
	backtrack2(str2, str1, move, len(str2), len(str1), "", "")
	
#glowna metoda, inicjalizuje tablice, wykonuje algorytm, wypisuje wynik(i)
def nw(str1, str2):
	tab = init(str1, str2)
	move = gen_array(len(str1)+1, len(str2)+1)
	fill(str1, str2, tab, move)
	print 'best score: ', tab[len(str2)][len(str1)]
	backtrack(str1, str2, move)

#pobiera przez entrez sekwencje aminokwasow o podanym id
def getSequence(gi):
	handle = Entrez.efetch(db="protein", id=gi, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	return record.seq

#main
def main():
	print 'starting...'
	Entrez.email = "i.orlowska@uj.edu.pl"
	print 'fetching data...'
	str1 = getSequence("40886941") # hemoglobina czlowieka
	str2 = getSequence("34849618") # szczura
	print 'matching...'
	nw(str1, str2)
	print 'end'

if __name__ == "__main__":
	main()


