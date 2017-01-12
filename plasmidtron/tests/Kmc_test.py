import unittest
import os
from plasmidtron.Kmc import Kmc
from plasmidtron.SampleData import SampleData

class TestKmc(unittest.TestCase):
	
	def test_basename_set_correctly(self):
		'''test basename extracted correctly from input file with gz'''
		i = SampleData('/path/to/sample_1.fastq.gz','/path/to/sample_2.fastq.gz' )
		self.assertEqual(i.basename, 'sample')

	def test_populate_fofn_name(self):
		sample = SampleData('/path/to/sample_1.fastq.gz','/path/to/sample_2.fastq.gz' )
		i = Kmc(os.getcwd(), sample, 1, 51, 30)
		self.assertEqual(sample.file_of_fastq_files, i.temp_working_dir+"/fofn")
	
	def test_create_file_of_file_names(self):
		sample = SampleData('/path/to/sample_1.fastq.gz','/path/to/sample_2.fastq.gz' )
		i = Kmc(os.getcwd(), sample, 1, 51, 30)
		i.create_file_of_file_names(sample.file_of_fastq_files)
		
		with open(sample.file_of_fastq_files, 'r') as actual_file:
			actual_fofn_content = actual_file.read()
			expected_fofn_content = """\
/path/to/sample_1.fastq.gz
/path/to/sample_2.fastq.gz
"""
			self.assertEqual(actual_fofn_content, expected_fofn_content)
			