'''Run KMC over a single sample (2 FASTQs) to get kmers database'''

class Kmc:
	def __init__(self,output_directory, sample_data, verbose, threads, kmer, min_kmers_threshold):
		self.output_directory = output_directory
		self.sample_data = sample_data
		self.verbose = verbose
		self.threads = threads
		self.kmer = kmer
		self.min_kmers_threshold = min_kmers_threshold
		
	def run(self):	
		temp_working_dir = tempfile.mkdtemp(dir=self.output_directory)
		# TODO make sure to strip off path
		basename = sample.forward_file.replace('_1.fastq.gz','')
		sample.basename = basename
	
		with open(temp_working_dir+'/fofn', 'w') as file_of_sample_fastqs:
			file_of_sample_fastqs.write(sample.forward_file + "\n")
			file_of_sample_fastqs.write(sample.reverse_file + "\n")
	
		database_name =  temp_working_dir+"/kmc_"+basename
		file_of_fastq_files = temp_working_dir+"/fofn"
	
		sample.database_name = database_name
		sample.file_of_fastq_files =file_of_fastq_files
	
		kmc_command = "kmc -t"+str(self.threads)+" -ci"+str(self.min_kmers_threshold)+" -k"+ str(self.kmer) +" @"+file_of_fastq_files+" "+ temp_working_dir+"/kmc_"+basename +" "+ temp_working_dir
		print('DEBUG: '+ kmc_command)
		subprocess.call(kmc_command,shell=True)
