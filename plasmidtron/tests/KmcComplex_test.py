import unittest
import os
from plasmidtron.KmcComplex import KmcComplex
from plasmidtron.SampleData import SampleData

class TestKmcComplex(unittest.TestCase):
	
	def test_generate_complex_config_filename(self):
		k = KmcComplex(os.getcwd(),1,1,[],[], 'union')
		# a default should be set.
		self.assertTrue(k.complex_config_filename != '')
		
		k.temp_working_dir = '/path/to'
		self.assertEqual(k.generate_complex_config_filename(),'/path/to/complex_config_file')

	def test_sample_definition_line(self):
		k = KmcComplex(os.getcwd(),1,1,[],[], 'union')
		s = SampleData('/path/to/sample#ABC_1.fastq','/path/to/sample#ABC_2.fastq' )

		self.assertEqual(k.sample_definition_line(s),'sample_ABC = kmc_sample#ABC')
		k.cleanup()

	def test_sample_definitions_str(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'union')
		expected_output = """\
a = kmc_a
b = kmc_b
c = kmc_c
"""
		self.assertEqual(k.sample_definitions_str(),expected_output)
		k.cleanup()
		
	def test_samples_to_set_operation_str(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'union')
		self.assertEqual(k.samples_to_set_operation_str(),'result = (a+b)-(c)')
		k.cleanup()
		
	def test_create_config_file_union(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'union')
		k.create_config_file()
		
		with open(k.complex_config_filename, 'r') as actual_file:
			actual_config_content = actual_file.read()
		self.assertEqual(actual_config_content, """\
INPUT:
a = kmc_a
b = kmc_b
c = kmc_c
OUTPUT:
result = (a+b)-(c)
OUTPUT_PARAMS:
-ci1
""")
		
		k.cleanup()
		
		
	def test_create_config_file_intersection(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'intersection')
		k.create_config_file()

		with open(k.complex_config_filename, 'r') as actual_file:
			actual_config_content = actual_file.read()
		self.assertEqual(actual_config_content, """\
INPUT:
a = kmc_a
b = kmc_b
c = kmc_c
OUTPUT:
result = (a*b)-(c)
OUTPUT_PARAMS:
-ci1
""")
		
		k.cleanup()

	def test_kmc_complex_command(self):
		k = KmcComplex(os.getcwd(), 1, 1, [], [], 'union')
		k.complex_config_filename = '/path/to/config'
		self.assertEqual(k.kmc_complex_command(), 'kmc_tools -t1 complex complex_config_file')
		k.cleanup()
	

