import unittest
import os
import tempfile
from plasmidtron.PlotKmers import PlotKmers

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plotkmers')

class TestPlotKmers(unittest.TestCase):
	
	def test_get_kmers_from_db(self):
		'''Given a database extract the kmers into an array'''
		temp_working_dir = tempfile.mkdtemp(dir=os.path.abspath(data_dir))
		p = PlotKmers([],temp_working_dir,1,31,255, False)
		kmers  = p.get_kmers_from_db(os.path.join(data_dir,'simple_kmers'))
		self.assertEqual(kmers, ['AAAAAAAAAAAGAAAAAAAAAAAAA', 'AAAAAAAAAAGAAAAAAAAAAAAAA'])