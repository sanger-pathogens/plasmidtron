import unittest
import filecmp
import os
import re
import plasmidtron

modules_dir = os.path.dirname(os.path.abspath(spreadsheetparser.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'spreadsheet')

class TestCheckSpreadsheet(unittest.TestCase):
    
    def test_invalid_spreadsheet_doesnt_exist(self):
        '''test_invalid_spreadsheet_doesnt_exist'''
        i = SpreadsheetParser(os.path.join(data_dir, 'file_which_doesnt_exist'))
        self.assertFalse(i.is_valid_file())
        self.assertEqual(i.extract_samples(),[])
		