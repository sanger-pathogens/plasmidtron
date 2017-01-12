import sys
import os
import csv
import logging
from plasmidtron.SampleData import SampleData

class SpreadsheetParser:
	def __init__(self,filename):
		self.filename = filename
		self.logger = logging.getLogger(__name__)
	
	def extract_samples(self):
		samples = []
		self.logger.info("Reading input spreadsheet")
		with open(self.filename) as csvfile:
			spreadsheetreader = csv.reader(csvfile, delimiter = ',')
			for row in spreadsheetreader:
				forward_file = row[0]
				reverse_file = row[1]

				for filename in [forward_file, reverse_file]:
					if not os.path.exists(filename):
						raise Exception('Input file in spreadsheet doesnt exit: '+ filename)
						
				self.logger.info("Found input files")
				samples.append( SampleData(forward_file,reverse_file) )
		return samples
		