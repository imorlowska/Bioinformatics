from Bio.Blast import NCBIWWW, NCBIXML

seq = """ >Lookup Sequence (in FASTA format)
MKSILDGLADTTFRTITTDLLGSPFQEKMTAGDNPQLVPADQVNITEFYNKSLSSFKENEENIQCGENFMDIECFMVLNPSQQLAIAVLSLTLGTFTVLENLLVLCVILHSRSLRCRPSYHFIGSLAVADLLGSVIFVYSFIDFHVFHRKDSRNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLAYKRIVTRPKAVVAFCLMWTIAIVIAVLPLLGWNCEKLQSVCSDIFPHIDETYLMFWIGVTSVLLLFIVYAYMYILWKAHSHAVRMIQRGTQKSIIIHTSEDGKVQVTRPDQARMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNKLIKTVFAFCSMLCLLNSTVNPIIYALRSKDLRHAFRSMFPSCEGTAQPLDNSMGDSDCLHKHANNAASVHRAAESCIKSTVKIAKVTMSVSTDTSAEAL

"""

def main():
	print 'querying...'
	qblast_output = NCBIWWW.qblast("blastp", "nr", seq, hitlist_size=50) #zapytanie
	blast_records = NCBIXML.parse(qblast_output) #parsowanie
	blast_record = blast_records.next() #bylo tylko jedno zapytanie
	print 'results: '
	ile = 0
	len_sum = 0
	e_ile = 0
	e_sum = 0

	for a in blast_record.alignments: #dla kazdego z wynikow
		ile += 1
		print 'Title: ', a.title
		print 'Length:   ', a.length
		len_sum += a.length
		for hsp in a.hsps:
			e_sum += hsp.expect
			e_ile += 1
	if ile > 0:
		print 'Average length: ', len_sum/ile # srednia dlugosc
	else:
		print 'Average length: ', 0
	if e_sum > 0:
		print 'Average E: ', e_sum/e_ile #srednie e
	else:
		print 'Average E: ', 0

if __name__ == "__main__":
	main()
