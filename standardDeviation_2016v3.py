##################################################################################
# Algorithm to calculate pan-genome growth. Calculate means and standard errors  #
# for a determined number of genome sets. For example, for 29 genomes, calculate #
# the averages and standard deviations for each combination of 1-[determined     #
# number of genomes]. Program reads in pickled dictionary with keys that are     #
# genome intersections and values that are number of genes and gene identifiers. #
# Algorithm:																	 #
#	1. Make list of genome trilogons	                                         #
#	2. Make all genome combinations up to specified value						 #
#	2. Find trilogon combination in dictionary keys                              #
#	3. Sum gene counts for combination                                           #
#	4. Create dictionary of sets with overlapping genes                          #
#	5. dictionary of set lengths built from sums of set dictionary               #
#																				 #
# In version 2, the algorithm has modified to use all 29 genomes and when it     #
# does the gene couting it only counts the first 1000 permutations.              #
# 																				 #
##################################################################################

import sys
from os import listdir
from os.path import isfile, join
import pickle
from itertools import combinations, chain
import math
import random

#makes a dictionary to associate integers with genome trilogons
def makeGenomeDict():
	triD = {}
	counter = 1
	for f in listdir('../flowerplot.v2/genomes/'):
		if isfile(join('../flowerplot.v2/genomes',f)) and '.fa' in f and not '.gz' in f:
			fo = open('../flowerplot.v2/genomes/'+f,'r')
			triValue = fo.readline()[1:4]
			fo.close()
			triD[counter] = triValue
			counter += 1
	return triD
	
#make a dictionary to associate all genome keys in "gIntersect" with a single genome
def makeGeneSets(dict, intersection):
	gSetKeys = {}
	for genome in dict:
		gSetKeys[genome] = set([])
		for keys in intersection:
			if dict[genome] in keys:
				gSetKeys[genome].add(keys)
	return gSetKeys

#make all combinations of the genomes (from genome dictionary) and return an iterator
def makePermutations(dict, ngenomes):
	list = random.sample(dict.keys(),ngenomes)
	return chain.from_iterable(combinations(list,r) for r in range(ngenomes+1))

#function to calculate the gene counts for each permutation. giter is the iterator, dict is the gene set dictionary
#and intersection is the pickled results file
def countGenes(subL, dict, intersection, l):
	permD={}
	i = 0 #variable to count the number of permutations evaluated
	binom = [] #list to store the binomial coefficient for each combination. 
	for j in range(ngenomes+1):
		binom.append(math.factorial(ngenomes)/(math.factorial(j)*math.factorial(ngenomes-j)))
	if binom[l] > 1000000:
		ii = round(binom[l]/100000)
	else:
		ii = 1
	
	for element in subL: 
		if i == ii: 
			permD[element] = 0
			gSet = set([])
			for k in range(len(element)):
				gSet = gSet.union(dict[element[k]])
			for combo in gSet:
				permD[element] += intersection[combo][2]
			i = 0
		i += 1
	return permD

def main():
	deviation_file = open('deviations.txt', 'w')
	avgs_file = open('averages.txt', 'w')
	triD = makeGenomeDict()
	gSetKeys = makeGeneSets(triD, gIntersect)
	giter = makePermutations(triD, ngenomes)
	giter.__next__()
	
	sublist = []
	avg=[]
	sdev=[]
	c1=0
	c2=2
	for element in giter:
		c1 = len(element)
		if c1 == c2:
			print('working on ',len(element))
			avgN = 0
			s_devN = 0
			permD = countGenes(sublist, gSetKeys, gIntersect, len(element)-1) #make list of counts
			for key, value in permD.items():
				avgN += value
			avgN = avgN/len(permD)
			for key, value in permD.items():
				s_devN += (value - avgN)**2
			s_devN = math.sqrt(s_devN/(len(permD)-1))
			avg.append(avgN)
			sdev.append(s_devN)
			sublist  = []
			pass
		if len(element) > 0:
			sublist.append(element)
		c2 = c1 + 1
					
	for i in avg:
		avgs_file.write(i+'\n')
	for i in sdev:
		deviation_file.write(i+'\n')

	deviation_file.close()
	avgs_file.close()
	
gIntersect = pickle.load(open('../flowerplot.v2/result_all/result','rb'))
ngenomes = 29
main()
