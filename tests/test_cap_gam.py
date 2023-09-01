import pytest
import unittest
import numpy as np
import pandas as pd
from collections.abc import Iterable
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

def test_capgam_returns_DataFrame():
    assert type(e.capgam(edata,"Si29")) is pd.core.frame.DataFrame
    assert type(e.capgam(edata,"La140","more")) is pd.core.frame.DataFrame

class CapGamTests(unittest.TestCase):

    __doc__ = """Unit tests for the `capgam` method of the CapGam class."""
    
    def test_capgam_returns_NoneType_if_no_gammaspe_data(self):
        return_value_capgam = e.capgam(edata,"Si28")
        self.assertIsNone(return_value_capgam)
        return_value_capgam = e.capgam(edata,"Se70","more")
        self.assertIsNone(return_value_capgam)
        return_value_capgam = e.capgam(edata,"Kr88")
        self.assertIsNone(return_value_capgam)
        return_value_capgam = e.capgam(edata,"Cs119","more")
        self.assertIsNone(return_value_capgam)

    def test_contents_capgam_df(self):
        df = e.capgam(edata,"C13")
        d = df.to_dict('list')
        assert isinstance(d, Iterable)
        for (key,value) in d.items():
            if key == 'Type':
                testValue = False
                for v in value:
                    if (v == 'primary') or (v == 'secondary'):
                        testValue = True
                    self.assertTrue(testValue,"")
                    self.assertIsInstance(v, str)
            else:
                for v in value:
                    self.assertIsInstance(v, float)

    def test_contents_capgam_df_for_every_EGAF_residual(self):
        """This function resulted in the addition of a TypeError handling 
        exception in the capgam method to account for compound-nucleus 
        residuals containing only one type of gamma: Either primaries or 
        secondaries.  The following list of residuals contain only one type of
        gamma:
        
        Cd112, Ce139, Er163, H2, H3, He4, Hg197, La139, Pd103, Pd105, Ru97, 
        Ru99, Sm155, Sn120, Sn123, Sn125, Sr85, Te126, W181, Xe125, Xe129, Xe131
        """
        res = e.egaf_residual_list(edata)
        for r in res:
            df = e.capgam(edata,r)
            d = df.to_dict('list')
            assert isinstance(d, Iterable)
            for (key,value) in d.items():
                if key == 'Type':
                    testValue = False
                    for v in value:
                        if (v == 'primary') or (v == 'secondary'):
                            testValue = True
                        self.assertTrue(testValue,"")
                        self.assertIsInstance(v, str)
                else:
                    for v in value:
                        self.assertIsInstance(v, float)

    def test_more_contents_capgam_df(self):
        df = e.capgam(edata,"C13","more")
        d = df.to_dict('list')
        assert isinstance(d, Iterable)
        for (key,value) in d.items():
            if key == 'Type':
                testValue = False
                for v in value:
                    if (v == 'primary') or (v == 'secondary'):
                        testValue = True
                    self.assertTrue(testValue,"")
                    self.assertIsInstance(v, str)
            elif (key == 'i' or key == 'f'):
                for v in value:
                    self.assertIsInstance(v, int)
            else:
                for v in value:
                    self.assertIsInstance(v, float)                        

    def test_more_contents_capgam_df_for_every_EGAF_residual(self):
        res = e.egaf_residual_list(edata)
        for r in res:
            df = e.capgam(edata,r,"more")
            d = df.to_dict('list')
            assert isinstance(d, Iterable)
            for (key,value) in d.items():
                if key == 'Type':
                    testValue = False
                    for v in value:
                        if (v == 'primary') or (v == 'secondary'):
                            testValue = True
                        self.assertTrue(testValue,"")
                        self.assertIsInstance(v, str)
                elif (key == 'i' or key == 'f'):
                    for v in value:
                        self.assertIsInstance(v, int)
                else:
                    for v in value:
                        self.assertIsInstance(v, float)
