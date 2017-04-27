import unittest
import os
import re
from plasmidtron.KmcVersionDetect import KmcVersionDetect

'''KMC command syntax changed for kmc_tools between version 2 and 3 so detect the version here'''
class TestKmcVersionDetect(unittest.TestCase):
	
	'''Can the version of KMC be extracted'''
	def test_kmc_command(self):
		k = KmcVersionDetect(False)
		self.assertTrue(k.major_version())
		
		full_version = k.full_version()
		version_search_results = re.search("([\d]+\.[\d]+\.[\d]+)", full_version)
		self.assertTrue(version_search_results.group(1))
		self.assertTrue(k.major_version())

