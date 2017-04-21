import os
import re 
 
class SampleData:
	def __init__(self,forward_file, reverse_file = None):
		self.forward_file = forward_file
		self.reverse_file = reverse_file
		self.file_of_fastq_files = ''
		self.basename = self.calculate_basename(forward_file)
		self.filtered_forward_file=''
		self.filtered_reverse_file=''
		self.database_name = 'kmc_'+self.basename
	     
	def calculate_basename(self,filename):
		basename = os.path.basename(filename)
		basename = basename.replace('.gz','')
		basename = basename.replace('_1.fastq','')
		return basename
		
	def is_a_fasta(self):
		m = re.search('\.(fasta|fa|fsa|fna)(.gz)?$', self.forward_file)
		if m and m.group(0):
			return True
		else:
			return False
	
	def cleanup(self):
		os.remove(self.filtered_forward_file)
		os.remove(self.filtered_reverse_file)
		os.remove(self.file_of_fastq_files)
