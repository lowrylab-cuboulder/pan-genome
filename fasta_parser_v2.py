#########################################################################
# NOTE: NCBI is phasing out GI numbers in September 2016				#
# An HTML parser that will find fasta format genes from a gene name		#
# ALGORITHM:															#
#	1. find fasta file													#
#	2. blast against 29 genome database									#
#	3. evaluate BLAST results to determine where the hit came from and	#
#		and determine if it's worth using								#
#	4. compile results													#
#########################################################################


import sys
import urllib.request
import html.parser
from bs4 import BeautifulSoup
from os import system
from os import listdir
from os.path import isfile, join
import csv

#get list of genes
listname = input('enter list file name: ')
glist = []
with open('./gene_lists/'+listname, 'r') as f:
	for line in f:
		if '\n' in line:
			glist.append(line[:-1])
		else:
			glist.append(line)
f.close()			

#get fasta format file from a gene search
#NEED TO INCLUDE SOME ERROR HANDLING AND WARNINGS IF CANT FIND GENE
def uniprotUrl(name):
	f = open('./tmp/seq.fa','w')
	baseURL = "http://www.uniprot.org/uniprot/"
	#the search variable may need to be edited upon search request
	search = "http://www.uniprot.org/uniprot/?query="+name+'+taxonomy%3A"mycobacterium"&sort=score' #added the taxonomy filter 092816
	with urllib.request.urlopen(search) as search_file:
		soup = BeautifulSoup(search_file, 'html.parser')
		num = 0
		if soup.find(id='noResultsMessage'):
			print('cant find '+name)
			f.close()
			progress = 0
			return progress
		for things in soup.find_all(class_='entry'):
			if num != 1:
				for lines in things.children:
					#if 'Mycobacterium' in str(lines.a): #commented because overlaps with "Actinobacteria"
					geneID = things.find(class_='entryID').string
					num = 1
		newURL = baseURL+geneID+'.fasta'
		with urllib.request.urlopen(newURL) as faFile:
			soup2 = BeautifulSoup(faFile, 'html.parser')
			fasta = soup2.string
			f.write(fasta)
	progress = 1
	f.close()
	return progress
	

#BLAST a fasta file against 29 genome database
def BLAST(eval):
	system('.\\blastp -query ' + '.\\tmp\\seq.fa -db .\\database -evalue ' + str(eval) + ' -outfmt \"6 sseqid pident qcovhsp gapopen evalue bitscore\" -out .\\tmp\\res.txt')
	f = open('./tmp/res.txt', 'r')
	num = 0
	if len(f.readline()) > 0:
		num = 1
	f.close()
	return num

#make a dictionary to associate the trilogons with genome name
def makeTriD():
	triD = {}
	for f in listdir('../flower_plot+bsr/Mycobacterial_Pangenome/flowerplot.v2/genomes/'):
		if isfile(join('../flower_plot+bsr/Mycobacterial_Pangenome/flowerplot.v2/genomes',f)) and '.fa' in f and not '.gz' in f:
			name = f.split('.')[1]
			fo = open('../flower_plot+bsr/Mycobacterial_Pangenome/flowerplot.v2/genomes/'+f,'r')
			triSet = set([])
			for line in fo:
				if line[0] == '>':
					triSet.add(line[1:4])
			fo.close()
			for i in triSet:		#this should probably be warning of some kind
				if i in triD.keys():
					print(i,' already exists as ',triD[i])
				else:
					triD[i] = name
	return triD

#make a dictionary to have a permanent order to the genomes
def orderGenomes(dict):
	counter = 0
	newDict = {}
	for key, value in dict.items():
		if value not in newDict.keys():
			counter += 1
			newDict[value] = counter
	return newDict

#get the max bitscore
def maxBit():
	f = open('./tmp/res.txt','r')
	max = float(f.readline().split('\t')[5])
	return max

#evaluate BLAST results and order the data based on genome order. return a list of whether or
#not a homolog is present.
def writeData(TriD, posiD, gname):
	f = open('./tmp/res.txt','r')
	maxbit = maxBit()
	nameList = []
	valueD = {}
	theDATA = [gname]
	for line in f:
		res = line.split('\t')
		bitscore = float(res[5])/maxbit
		name = TriD[res[0][0:3]]
		if name not in nameList:
			nameList.append(name)
			valueD[posiD[name]] = int(res[2])*bitscore/100
	for key, value in posiD.items():
		if key not in nameList:
			valueD[posiD[key]] = 0
	f.close()
	for i in range(29):
		theDATA.append(valueD[i+1])
	return theDATA
				
eval = input('enter evalue: ')
cover = input('enter coverage minimum: ')
minbit = input('enter minimum bit score ratio: ')
TriD = makeTriD()
posiD = orderGenomes(TriD)
csvfile = open(listname[:-4]+'_results.csv', 'w', newline='')
reswriter = csv.writer(csvfile, delimiter = ',', quotechar = '|', quoting = csv.QUOTE_MINIMAL)
header = [0]*30
header[0] = 'gene_name'
for key, value in posiD.items():
	header[value] = key
reswriter.writerow(header)
for gene in glist:
	print(gene)
	progress = uniprotUrl(gene)
	if progress == 1:
		BLASTres = BLAST(eval)
		if BLASTres == 1:
			data = writeData(TriD, posiD, gene)
			reswriter.writerow(data)
		else:
			print(gene+' does not exisit in genome database...')
csvfile.close()	