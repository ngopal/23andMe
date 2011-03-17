###################################################################################
# Module written to help me explore 23andMe raw data sets 
# Written by Nikhil Gopal 2-6-2011                        
# http://blog.nikhilgopal.com                              
#
# Copyright (c) 2011 Nikhil Gopal
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###################################################################################

import os, sys
import sqlite3
import json
from json import JSONEncoder

class ParseToDB:
	'''
	This class loads the input files into a sqlite database
	'''
	def __init__(self,list_of_files):
		'''
		This function parses the 23andMe raw data file to a sqlite database
		'''
		dbName = 'AllData.db'
		try:
			self.createDB(dbName,list_of_files)
		except sqlite3.OperationalError:
			print "Database", dbName, "already exists. Moving on..."
			return None
		for i in list_of_files:
				self.readInFilesToDB(i,dbName)
		
	def createDB(self,name,list_of_files):
		'''
		This function creates a database
		'''
		DB = sqlite3.connect(name)
		cursor = DB.cursor()
		for i in list_of_files:
			TableName = str(i.split('.')[0])
			cursor.execute('create table '+TableName+' (rsid TEXT, chr TEXT, pos TEXT, geno TEXT)')
		DB.close()
	
	def readInFilesToDB(self,infile,dbName):
		'''
		This function reads in the 23andMe raw data file into a dictionary
		'''
		f = open(infile,'r')
		DB = sqlite3.connect(dbName)
		cursor = DB.cursor()
		for i in f.readlines():
			if '#' not in i:
				line = i.strip().split('\t')
				TableName = str(infile.split('.')[0])
				cursor.execute("insert into "+TableName+" values(?,?,?,?)", (line[0],line[1:][0],line[1:][1],line[1:][2]))
		DB.commit()
		DB.close()
		f.close()


class ParseToJSON:
	'''
	This function parses input data into json format
	'''
	def __init__(self,list_of_files):
		'''
		This function reads the input data into a dictionary, converts it to json, and prints it to an outfile
		'''
		outfile = open('23andMe.json', 'w')
		self.Data = {}
		for i in list_of_files:
				self.Data[i] = self.readInFile(i)
		self.Data = self.convertToJSON(self.Data)
		outfile.write(str(self.Data))
		outfile.close()

	def readInFile(self,infile):
		'''
		This function reads in the data from a file into a dictionary
		'''
		f = open(infile, 'r')
		rsid = {}
		for i in f.readlines():
			if '#' not in i:
				line = i.strip().split('\t')
				rsid[line[0]] = line[1:]
		f.close()
		return rsid

	def convertToJSON(self,data):
		'''
		This function converts the data structure (dictionary) into a json string
		'''
		out = JSONEncoder().encode(data)
		return out


class ParseToDict:
	'''
	This class reads in the raw data files into a dictionary
	'''
	def __init__(self,list_of_files):
		'''
		This function initializes the data contained in the class
		'''
		# The main data dictionary
		self.Data = {}
		for i in list_of_files:
			self.Data[i] = self.readInFile(i)
		
		# A data dictionary with a list of rsids encountered in each file
		self.RSids = self.getRSids()
		
		# Orders files by their number of RSids
		self.FileRSidsTuple = sorted([(i,len(self.RSids[i])) for i in self.RSids], key = lambda a: -a[1])
		
		# File with the maximum number of RSids
		self.MaxSize = self.FileRSidsTuple[0]
		
		# File with the minimum number of RSids
		self.MinSize = self.FileRSidsTuple[-1]
		
		# Find intersection between the input datasets
		self.Intersection = self.calcIntersectionAll()
		
		# A data dictionary which include only the snps in the intersection list
		self.intersectionData = self.commonDict()
		
		# A variable that contains the files provided as input
		self.files = [i for i in set(self.Data.keys())]
		
		
	def readInFile(self,infile):
		'''
		This function reads in the data from a file into a dictionary
		'''
		f = open(infile, 'r')
		rsid = {}
		for i in f.readlines():
			if '#' not in i:
				line = i.strip().split('\t')
				rsid[line[0]] = line[1:]
		f.close()
		return rsid
		
	def getRSids(self):
		'''
		This function returns a dictionary identical to self.Data, but each file (key) returns a list of corresponding RSids
		'''
		self.Sizes = {}
		for i in self.Data:
			self.Sizes[i] = self.Data[i].keys()
		return self.Sizes

	def getIntersection(self,A,B):
		'''
		This function returns the intersection of two lists
		'''
		self.inter = set(A).intersection( set(B) )
		return self.inter
		
	def searchSNP(self,rsid):
		'''
		This function looks for SNPS across the input files
		'''
		list = []
		for i in self.Data:
			try:
				tuple = (i, self.Data[i][rsid][-1])
				list.append(tuple)
			except:
				tuple = (i, 'N/A')
				list.append(tuple)
		return list
		
	def searchChromosomes(self,chrom):
		'''
		This function outputs all of the data across all of the files for a specified chromosome
		'''
		self.biglist = []
		for i in self.intersectionData:
			for k in self.intersectionData[i]:
				if str(chrom) == str(self.intersectionData[i][k][0]) and str(min) < str(self.intersectionData[i][k][1]) and str(max) > str(self.intersectionData[i][k][1]):
					li = []
					for j in self.intersectionData:
						tup = (j,self.intersectionData[i][k])
						li.append(tup)
					self.biglist.append(li)
		return self.biglist
		
	def calcIntersectionAll(self):
		'''
		This function calculates the intersection between all of the input files
		'''
		self.li = set(self.Data[self.MinSize[0]].keys())
		for i in range(len(self.files)):
			self.li = (set(self.li) & set(self.Data[self.files[i]].keys()))
		return self.li
		
	def commonDict(self):
		'''
		This function outputs a dictionary (identical to the Data dictionary) which only contains the snps from the intersection list
		'''
		self.newDict = {}
		for i in self.Data:
			self.rsiddict = {}
			for k in self.Data[i]:
				if k in self.Intersection:
					self.rsiddict[k] = self.Data[i][k]
				else:
					pass
			self.newDict[i] = self.rsiddict
		return self.newDict

	def identity(self, dict1, dict2):
		'''
		This function compares two dictionaries and returns a dictionary which contains only the identical elements between
		'''
	    return dict( (key, dict1[key]) for key in (set(dict1) & set(dict2)) if dict1[key] == dict2[key] )

	def halfIdentity(self, dict1, dict2):
		'''
		This function compares two dictionaries and returns a dictionary which contains only the half identical elements between
		'''
	    return dict( (key, dict1[key]) for key in (set(dict1) & set(dict2)) if dict1[key][-1][0] == dict2[key][-1][0] or dict1[key][-1][0] == dict2[key][-1][1] or dict1[key][-1][1] == dict2[key][-1][0] or dict1[key][-1][1] == dict2[key][-1][1])

	def phylogeny(self,rsid=None,metric='identity'):
		'''
		This function generates a phylogeny from the raw data files. If rsid is lookup is requested, each raw data file in the phylogeny can be coupled with its genotype
		'''
		def level(x):
			'''
			This function essentially removes the tuples and sub-lists contained in a list and returns a single flattened list
			'''
		    result = []
		    for n in x:
		        if hasattr(n, "__iter__") and not isinstance(n, basestring):
		            result.extend(level(n))
		        else:
		            result.append(n)
		    return result

		def makeTree(li):
			'''
			This function creates a visual clustering effect on the final output list for the phylogeny function
			'''
			for i in range(len(li)):
				if i == 0 or i ==1:
					tup = [li[0],li[1]]
				else:
					tup = [tup,li[i]]
			return tup


		self.threshold = {}
		self.people = [i for i in self.files]
		self.combo = {}

		# Prepare the ModDict dictionary by removing the SNPS that are not in common
		self.ModDict = self.commonDict()

		# Generate a list of combinations and scores
		for i in range(len(self.people)):
			for k in range(len(self.people)):
				if self.people[i] == self.people[k]:
					pass
				else:
					# a different metric is used if specified, otherwise identity is the default
					if metric == 'identity':
						self.common = self.identity(self.ModDict[self.people[i]],self.ModDict[self.people[k]])
					elif metric == 'halfidentity':
						self.common = self.halfIdentity(self.ModDict[self.people[i]],self.ModDict[self.people[k]])
					else:
						self.common = self.identity(self.ModDict[self.people[i]],self.ModDict[self.people[k]])
					self.combo[(self.people[i], self.people[k])] = len(self.common)
					# record the highest score to identify which two files are first clustered
					if len(self.common) > len(self.threshold):
						self.threshold == self.common.copy()
						self.coords = (i,k)
					else:
						pass

		# Create the phylogeny
		self.phylo = []
		self.used = []
		for i in sorted(self.combo.items(), key=lambda x: x[1], reverse=True):
			if len(self.phylo) == 0:
				self.phylo.append(i[0])
				self.used.append(i[0][0])
				self.used.append(i[0][1])
			else:
				if i[0][0] in self.used:
					pass
				else:
					self.phylo.append(i[0][0])
					self.used.append(i[0][0])

		self.flatphy = level(self.phylo)
		
		# Return the genotype tupled with the raw file if rsid provided
		if rsid != None:
			tree = zip(self.flatphy,[Datasets.Data[i][rsid][-1] for i in self.flatphy if isinstance(i,str)])
			return makeTree(tree)
		else:
			return makeTree(self.flatphy)

if __name__ == "__main__":
	
	# Program takes 23andMe raw data files as input and they are stored in a list
	files_to_process = sys.argv[1:]
	
	# Program takes 23andMe raw data files as input and store them in a dictionary
	Datasets = ParseToDict(sys.argv[1:])