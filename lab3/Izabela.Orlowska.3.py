from Bio import Entrez, SeqIO

def part1(term):
	handle = Entrez.esearch(db="gene", term=term)
	result = Entrez.read(handle)
	ids = result['IdList']
	foundid = 0
	print "{:>12} {:>10} {:>11} {:>10}".format('Name', 'Chromosome', 'MapLocation', 'Description')
	for id in ids:
		handle = Entrez.esummary(db="gene", id=id)
		r = Entrez.read(handle)[0]
		print "{:>12} {:>10} {:>11} {:>10}".format(r['Name'], r['Chromosome'], (r['MapLocation'] or 'none'), r['Description'])
		if r['Name'] == "BRCA1":
			foundid = id
	return foundid

def part2(id):
	if id == 0:
		print "\nCouldn't find BRCA1\n"
	else:
		print "\nFound BRCA1 gene id = ", id, "\n"

def part3(term):
	handle = Entrez.esearch(db="omim", term=term)
	result = Entrez.read(handle)
	ids = result['IdList']
	print 'OMIM results\' titles:'
	for id in ids:
		handle = Entrez.esummary(db="omim", id=id)
		r = Entrez.read(handle)[0]
		print r['Title']

def part4(id):
	handle = Entrez.elink(dbfrom="gene", db="protein", id=id)
	result = Entrez.read(handle)
	ids = result[0]['LinkSetDb'][0]['Link']
	print '\nProtein Gi\'s:'
	for id in ids:
		handle = Entrez.esummary(db="protein", id=id['Id'])
		r = Entrez.read(handle)
		print r[0]['Gi']

def part5():
	handle = Entrez.efetch(db="protein", id="121949022", rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	print "\nFound protein in genbank, details:"
	print "name: " + record.name
	print "seq: " + record.seq

def part6(term):
	handle = Entrez.esearch(db="snp", term="Homo sapiens AND " + term)
	result = Entrez.read(handle)
	ids = result['IdList'][:50]
	print "\nSNP results:"
	print "{:>10} {:>10} {:>10} {:>30}".format('SNP ID', 'SNP CLASS', 'GENE', 'CONTIGPOS')
	for id in ids:
		handle = Entrez.esummary(db="snp", id=id)
		r = Entrez.read(handle)[0]
		print "{:>10} {:>10} {:>10} {:>30}".format(r['SNP_ID'], r['SNP_CLASS'], (r['GENE'] or 'empty'), r['CONTIGPOS'])

def main():
	term = "BRCA1"
	Entrez.email = "i.orlowska@uj.edu.pl"
	BRCA1_gene_id = part1(term)
	part2(BRCA1_gene_id)
	part3(term)
	part4(BRCA1_gene_id)
	part5()
	part6(term)

if __name__ == "__main__":
	main()

# part7
# Entrez jest bardzo przydatnym narzedziem, przydalaby sie jednak mozliwosc latwego cachowania
# czesci danych (gdy brakuje dostepu do sieci). 
# Ze wzgledu na prostote kodu i wygode stosowanie rettype genbank i retmode text wydaje sie
# ladniejsza alternatywa niz XML (record.name vs result['Name'] i duze prawdopodobienstwo
# wystapienia literowek). 
# Wygodnym narzedziem przy uzywaniu Entrez wydaje sie interpreter ipython, poniewaz bledy
# zwracane same z siebie z blednego uzywania Entrez sa malo intuicyjne.
#