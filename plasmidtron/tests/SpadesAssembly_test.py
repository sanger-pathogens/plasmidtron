import unittest
import os
from plasmidtron.SpadesAssembly import SpadesAssembly
from plasmidtron.SampleData import SampleData

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','spadesassembly')

class TestSpadesAssembly(unittest.TestCase):
	
	def test_filtering_small_contigs(self):
		sample = SampleData('/path/to/sample_1.fastq.gz','/path/to/sample_2.fastq.gz' )
		
		s = SpadesAssembly(sample, 'abc', 1, 1, '', 20)
		s.remove_small_contigs(os.path.join(data_dir, 'assembly_with_small_contigs.fa'), os.path.join(data_dir, 'actual_assembly_without_small_contigs.fa'))
		
		with open(os.path.join(data_dir, 'actual_assembly_without_small_contigs.fa'), 'r') as actual_file, open(os.path.join(data_dir, 'expected_assembly_without_small_contigs.fa'), 'r') as expected_file:
			actual_config_content = actual_file.read()
			expected_config_content = expected_file.read()
			
			self.assertEqual(actual_config_content,expected_config_content)
		os.remove(os.path.join(data_dir, 'actual_assembly_without_small_contigs.fa'))
		