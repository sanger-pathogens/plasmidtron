import unittest
import os
import dill
from plasmidtron.SpreadsheetParser import SpreadsheetParser

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','spreadsheetparser')

class TestSpreadsheetParser(unittest.TestCase):
	
	def test_invalid_spreadsheet_doesnt_exist(self):
		'''test_invalid_spreadsheet_doesnt_exist'''
		i = SpreadsheetParser(os.path.join(data_dir, 'file_which_doesnt_exist'), False)
		
		with self.assertRaises(Exception):
			i.extract_samples()
		
	def test_spreadsheet_with_one_set_of_files(self):
		'''A spreadsheet with one set of files should give a single object back'''
		i = SpreadsheetParser(os.path.join(data_dir, 'spreadsheet_with_one_set_of_files.csv'), False)
		samples = i.extract_samples()
		self.assertEqual(len(samples),1)
		self.assertEqual(samples[0].forward_file,'plasmidtron/tests/data/spreadsheetparser/sampleA_1.fastq.gz')
		self.assertEqual(samples[0].reverse_file,'plasmidtron/tests/data/spreadsheetparser/sampleA_2.fastq.gz')
		
	def test_spreadsheet_with_a_non_paired_file(self):
		'''A spreadsheet with a single file should be okay'''
		i = SpreadsheetParser(os.path.join(data_dir, 'spreadsheet_with_a_non_paired_file.csv'), False)
		samples = i.extract_samples()

		self.assertEqual(len(samples),1)
		self.assertEqual(samples[0].forward_file,'plasmidtron/tests/data/spreadsheetparser/sampleA.fasta')
		self.assertEqual(samples[0].basename, 'sampleA')
		
	def test_spreadsheet_with_dash_in_filenames(self):
		'''A spreadsheet with one set of files should give a single object back'''
		i = SpreadsheetParser(os.path.join(data_dir, 'spreadsheet_with_dash_in_filenames.csv'), False)
		samples = i.extract_samples()
		self.assertEqual(len(samples),1)
		self.assertEqual(samples[0].forward_file,'plasmidtron/tests/data/spreadsheetparser/sample-A_1.fastq.gz')
		self.assertEqual(samples[0].reverse_file,'plasmidtron/tests/data/spreadsheetparser/sample-A_2.fastq.gz')
		self.assertEqual(samples[0].basename, 'sample_A')