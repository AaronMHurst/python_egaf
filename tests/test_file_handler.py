import pytest
import unittest
import numpy as np
import pandas as pd
from collections.abc import Iterable
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

class FileHandlerTests(unittest.TestCase):

    __doc__ = """Unit tests for the `get_ripl` method of the RIPL class, the 
    `get_json` method of the JSONFile class, and the `get_ensdf` method of the
    ENSDF class.
    """

    # RIPL file tests:
    def test_get_ripl_returns_list_for_all_EGAF_residuals(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_ripl_files = []
        for r in res:
            ripl = e.get_ripl(edata, r)
            list_ripl_files.append(ripl)
            self.assertIsInstance(ripl, list)
            assert isinstance(ripl, Iterable)
        self.assertEqual(len(res), len(list_ripl_files))
            
    def test_get_ripl_returns_None_if_not_in_EGAF_residuals(self):
        ripl = e.get_ripl(edata, "Se70")
        self.assertIsNone(ripl)

    # JSON file tests:
    def test_get_json_returns_dict_for_all_EGAF_residuals(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_json_files = []
        for r in res:
            js = e.get_json(edata, r)
            list_json_files.append(js)
            self.assertIsInstance(js, dict)
            assert isinstance(js, Iterable)
        self.assertEqual(len(res), len(list_json_files))            

    def test_get_json_returns_None_if_not_in_EGAF_residuals(self):
        js = e.get_json(edata, "Se70")
        self.assertIsNone(js)

    # ENSDF file tests:
    def test_get_ensdf_returns_list_for_all_EGAF_residuals(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_ensdf_files = []
        for r in res:
            ensdf = e.get_ensdf(edata, r)
            list_ensdf_files.append(ensdf)
            self.assertIsInstance(ensdf, list)
            assert isinstance(ensdf, Iterable)
            assert len(ensdf) > 0
        self.assertEqual(len(res), len(list_ensdf_files))            

    def test_get_ensdf_returns_None_if_not_in_EGAF_residuals(self):
        ensdf = e.get_ensdf(edata, "Se70")
        self.assertIsNone(ensdf)        
