import sys
import os
import csv
from plasmidtron.SampleData import SampleData

class SpreadsheetParser:
	def __init__(self,filename):
		self.filename = filename
	
	def extract_samples(self):
		samples = []
		with open(self.filename) as csvfile:
			spreadsheetreader = csv.reader(csvfile, delimiter = ',')
			for row in spreadsheetreader:
				forward_file = row[0]
				reverse_file = row[1]

				for filename in [forward_file, reverse_file]:
					if not os.path.exists(filename):
						raise Exception('Input file in spreadsheet doesnt exit: '+ filename)
						  
				samples.append( SampleData(forward_file,reverse_file) )
		return samples
		