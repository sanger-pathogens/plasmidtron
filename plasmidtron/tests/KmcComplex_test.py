import unittest
import os
from plasmidtron.KmcComplex import KmcComplex
from plasmidtron.SampleData import SampleData

class TestKmcComplex(unittest.TestCase):

	def test_sample_definition_line(self):
		k = KmcComplex(os.getcwd(),1,1,[],[], 'union', False)
		s = SampleData('/path/to/sample#ABC_1.fastq','/path/to/sample#ABC_2.fastq' )

		self.assertEqual(k.sample_definition_line(s),'sample_ABC = kmc_sample#ABC')
		k.cleanup()

	def test_sample_definitions_str(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'union', False)
		expected_output = """\
a = kmc_a
b = kmc_b
c = kmc_c
"""
		self.assertEqual(k.sample_definitions_str(),expected_output)
		k.cleanup()
		
	def test_samples_to_set_operation_str(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'union', False)
		self.assertEqual(k.trait_samples_to_set_operation_str(),'traits = a+b')
		self.assertEqual(k.nontrait_samples_to_set_operation_str(),'nontraits = c')
		k.cleanup()
		
	def test_create_config_file_union(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'union', False)
		k.create_config_files()
		
		with open(os.path.join(k.temp_working_dir, 'traits_config_file'), 'r') as actual_file:
			actual_config_content = actual_file.read()
			self.assertEqual(actual_config_content, """\
INPUT:
a = kmc_a
b = kmc_b
c = kmc_c
OUTPUT:
traits = a+b
OUTPUT_PARAMS:
-ci1
""")
		with open(os.path.join(k.temp_working_dir, 'nontraits_config_file'), 'r') as actual_file:
			actual_config_content = actual_file.read()
			self.assertEqual(actual_config_content, """\
INPUT:
a = kmc_a
b = kmc_b
c = kmc_c
OUTPUT:
nontraits = c
OUTPUT_PARAMS:
-ci1
""")
		with open(os.path.join(k.temp_working_dir, 'combined_config_file'), 'r') as actual_file:
			actual_config_content = actual_file.read()
			self.assertEqual(actual_config_content, """\
INPUT:
set1 = traits
set2 = nontraits
OUTPUT:
result = set1-set2
OUTPUT_PARAMS:
-ci1
""")
		
		k.cleanup()
		
		
	def test_create_config_file_intersection(self):
		k = KmcComplex(os.getcwd(), 1, 1, [SampleData('a_1.fastq', 'a_2.fastq'), SampleData('b_1.fastq', 'b_2.fastq')], [SampleData('c_1.fastq', 'c_2.fastq')], 'intersection', False)
		k.create_config_files()

		with open(os.path.join(k.temp_working_dir, 'traits_config_file'), 'r') as actual_file:
			actual_config_content = actual_file.read()
			self.assertEqual(actual_config_content, """\
INPUT:
a = kmc_a
b = kmc_b
c = kmc_c
OUTPUT:
traits = a*b
OUTPUT_PARAMS:
-ci1
""")
		with open(os.path.join(k.temp_working_dir, 'nontraits_config_file'), 'r') as actual_file:
			actual_config_content = actual_file.read()
			self.assertEqual(actual_config_content, """\
INPUT:
a = kmc_a
b = kmc_b
c = kmc_c
OUTPUT:
nontraits = c
OUTPUT_PARAMS:
-ci1
""")
		with open(os.path.join(k.temp_working_dir, 'combined_config_file'), 'r') as actual_file:
			actual_config_content = actual_file.read()
			self.assertEqual(actual_config_content, """\
INPUT:
set1 = traits
set2 = nontraits
OUTPUT:
result = set1-set2
OUTPUT_PARAMS:
-ci1
""")
		
		k.cleanup()

	def test_kmc_complex_command(self):
		k = KmcComplex(os.getcwd(), 1, 1, [], [], 'union', False)
		self.assertEqual(k.kmc_complex_command('complex_config_file'), 'kmc_tools -t1 complex complex_config_file > /dev/null 2>&1')
		k.cleanup()
