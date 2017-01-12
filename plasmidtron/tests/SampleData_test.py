import unittest
import os
from plasmidtron.SampleData import SampleData


class TestSampleData(unittest.TestCase):
	
	def test_basename_set_correctly(self):
		'''test basename extracted correctly from input file with gz'''
		i = SampleData('/path/to/sample_1.fastq.gz','/path/to/sample_2.fastq.gz' )
		self.assertEqual(i.basename, 'sample')
		
	def test_basename_set_correctly_non_gzipped(self):
		'''test basename extracted correctly from input file non gzipped'''
		i = SampleData('/path/to/sample_1.fastq','/path/to/sample_2.fastq' )
		self.assertEqual(i.basename, 'sample')
	
