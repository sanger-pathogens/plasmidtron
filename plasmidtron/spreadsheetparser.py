import sys
import os
import csv
import plasmidtron

class SpreadsheetParser:
	def __init__(self,filename):
		self.filename = filename
	
	def extract_samples():
		samples = []
		with open(self.filename) as csvfile:
			spreadsheetreader = csv.reader(csvfile, delimiter = ',')
			for row in spreadsheetreader:
				forward_file = row[0]
				reverse_file = row[1]

				if not os.path.exists(forward_file) or not os.path.exists(reverse_file):
					 sys.exit( "The input files do not exist: "+forward_file+ " "+reverse_file)
				samples.append( SampleData(forward_file,reverse_file) )
		return samples
		
	def is_valid_file():
		return true
		