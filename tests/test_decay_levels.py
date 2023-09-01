import pytest
import unittest
import numpy as np
import pandas as pd
from collections.abc import Iterable
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

class LevelsTests(unittest.TestCase):

    __doc__ = """Unit tests for the following methods of the Levels class:

    `get_residual_levels`;
    `find_multiple_jpi`;
    `find_unique_jpi`;
    `find_isomers`.
    """

    # Level scheme tests:
    def test_get_residual_levels_for_29Si_returns_list(self):
        levels_pass_str = e.get_residual_levels(edata,"Si29")
        self.assertIsInstance(levels_pass_str, list)
        assert isinstance(levels_pass_str, Iterable)
        levels_pass_two_ints = e.get_residual_levels(edata,14,29)
        self.assertIsInstance(levels_pass_two_ints, list)
        assert isinstance(levels_pass_two_ints, Iterable)

    def test_get_residual_levels_for_nucleus_not_in_EGAF_returns_None(self):
        levels_pass_str = e.get_residual_levels(edata,"Se70")
        self.assertIsNone(levels_pass_str)
        levels_pass_two_ints = e.get_residual_levels(edata,34,70)
        self.assertIsNone(levels_pass_two_ints)

    def test_get_residual_levels_for_illegal_string_returns_None(self):
        levels_pass_str = e.get_residual_levels(edata,"%^&#@Jnsjksd125ASB")
        self.assertIsNone(levels_pass_str)        

    def test_get_residual_levels_for_wrong_arguments_returns_None(self):
        levels_wrong_args = e.get_residual_levels(edata,6,"C13")
        self.assertIsNone(levels_wrong_args)
        levels_wrong_args = e.get_residual_levels(edata,73.8,93.2)
        self.assertIsNone(levels_wrong_args)
        levels_wrong_args = e.get_residual_levels(edata,"C13","Si29")
        self.assertIsNone(levels_wrong_args)

    def test_get_residual_levels_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.get_residual_levels(6, 13, edata)
        with self.assertRaises(TypeError):
            e.get_residual_levels("C13", edata)
        bad_dict = {'a':0, 'b':1, 'c':2}
        with self.assertRaises(TypeError):
            e.get_residual_levels(bad_dict, edata)

    def test_get_residual_levels_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.get_residual_levels(XXXdata, "Si29")            
            
    def test_get_residual_levels_raises_KeyError_if_bad_dict_items_in_list(self):   
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.get_residual_levels(bad_dict_items_in_list, "Si29")

    def test_contents_of_returned_list_get_residuals_all_EGAF(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_level_schemes = []
        for r in res:
            levels = e.get_residual_levels(edata, r)
            assert isinstance(levels, Iterable)
            list_level_schemes.append(levels)
            for level in levels:
                self.assertIsInstance(levels, list)
                self.assertIsInstance(level, list)
                self.assertIsInstance(level[0], int)
                self.assertIsInstance(level[1], float)
                self.assertIsInstance(level[2], float)
                self.assertIsInstance(level[3], int)
                self.assertIsInstance(level[4], int)
                try:
                    self.assertIsInstance(level[5], float)
                except AssertionError:
                    self.assertIsNone(level[5])
                try:
                    self.assertIsInstance(level[6], int)
                except AssertionError:
                    print(r,level)
                    raise
                self.assertIsInstance(level[7], int)
                self.assertIsInstance(level[8], int)
        self.assertEqual(len(list_level_schemes), len(res))
        assert isinstance(list_level_schemes, Iterable)
        

    # Multiple JPi tests
    def test_find_multiple_jpi_in_24Na_returns_list(self):
        mult_jpi_str = e.find_multiple_jpi(edata,"Na24")
        self.assertIsInstance(mult_jpi_str, list)
        assert isinstance(mult_jpi_str, Iterable)
        mult_jpi_ints = e.find_multiple_jpi(edata,11,24)
        self.assertIsInstance(mult_jpi_ints, list)
        assert isinstance(mult_jpi_ints, Iterable)

    def test_find_multiple_jpi_not_in_29Si_returns_None(self):        
        mult_jpi_str = e.find_multiple_jpi(edata,"Si29")
        self.assertIsNone(mult_jpi_str)

    def test_find_multiple_jpi_not_in_EGAF_returns_None(self):        
        mult_jpi_str = e.find_multiple_jpi(edata,"Kr88")
        self.assertIsNone(mult_jpi_str)

    def test_find_multiple_jpi_for_illegal_string_returns_None(self):
        mult_jpi = e.find_multiple_jpi(edata,"%^&#@Jnsjksd125ASB")
        self.assertIsNone(mult_jpi)        

    def test_find_multiple_jpi_for_wrong_arguments_returns_None(self):
        mult_jpi = e.get_residual_levels(edata,6,"C13")
        self.assertIsNone(mult_jpi)
        mult_jpi = e.get_residual_levels(edata,300.564,1.25673)
        self.assertIsNone(mult_jpi)
        mult_jpi = e.get_residual_levels(edata,"O17","La140",8)
        self.assertIsNone(mult_jpi)

    def test_find_multiple_jpi_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.find_multiple_jpi(6, 13, edata)
        with self.assertRaises(TypeError):
            e.find_multiple_jpi("C13", edata)
        bad_dict = {'a':0, 'b':1, 'c':2}
        with self.assertRaises(TypeError):
            e.find_multiple_jpi(bad_dict, edata)

    def test_find_multiple_jpi_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.find_multiple_jpi(XXXdata, "Si29")            
            
    def test_find_multiple_jpi_raises_KeyError_if_bad_dict_items_in_list(self):   
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.find_multiple_jpi(bad_dict_items_in_list, "Si29")
        
    def test_contents_of_returned_list_find_multiple_jpi_all_EGAF(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_level_schemes = []
        for r in res:
            levels = e.find_multiple_jpi(edata, r)
            try:
                assert isinstance(levels, Iterable)
            except AssertionError:
                assert isinstance(levels, type(None))
            list_level_schemes.append(levels)
            try:
                for level in levels:
                    self.assertIsInstance(levels, list)
                    self.assertIsInstance(level, list)
                    self.assertIsInstance(level[0], int)
                    self.assertIsInstance(level[1], float)
                    self.assertIsInstance(level[2], float)
                    self.assertIsInstance(level[3], int)
                    self.assertIsInstance(level[4], int)
                    try:
                        self.assertIsInstance(level[5], float)
                    except AssertionError:
                        self.assertIsNone(level[5])
                    self.assertIsInstance(level[6], int)
            except AssertionError:
                self.assertIsNone(levels)
            except TypeError:
                self.assertIsNone(levels)
        self.assertEqual(len(list_level_schemes), len(res))
        assert isinstance(list_level_schemes, Iterable)

    # Unique JPi tests
    def test_find_unique_jpi_in_24Na_returns_list(self):
        unique_jpi_str = e.find_unique_jpi(edata,"Na24")
        self.assertIsInstance(unique_jpi_str, list)
        unique_jpi_ints = e.find_unique_jpi(edata,11,24)
        self.assertIsInstance(unique_jpi_ints, list)

    def test_find_unique_jpi_not_in_EGAF_returns_None(self):
        unique_jpi_str = e.find_unique_jpi(edata,"Kr88")
        self.assertIsNone(unique_jpi_str)

    def test_find_unique_jpi_for_illegal_string_returns_None(self):
        unique_jpi = e.find_unique_jpi(edata,"%^&#@Jnsjksd125ASB")
        self.assertIsNone(unique_jpi)        

    def test_find_unique_jpi_for_wrong_arguments_returns_None(self):
        unique_jpi = e.find_unique_jpi(edata,6,"C13")
        self.assertIsNone(unique_jpi)
        unique_jpi = e.find_unique_jpi(edata,300.564,1.25673)
        self.assertIsNone(unique_jpi)
        unique_jpi = e.find_unique_jpi(edata,"O17","La140",8)
        self.assertIsNone(unique_jpi)

    def test_find_unique_jpi_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.find_unique_jpi(6, 13, edata)
        with self.assertRaises(TypeError):
            e.find_unique_jpi("C13", edata)
        bad_dict = {'a':0, 'b':1, 'c':2}
        with self.assertRaises(TypeError):
            e.find_unique_jpi(bad_dict, edata)

    def test_find_unique_jpi_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.find_unique_jpi(XXXdata, "Si29")            
            
    def test_find_unique_jpi_raises_KeyError_if_bad_dict_items_in_list(self):   
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.find_unique_jpi(bad_dict_items_in_list, "Si29")
        
    def test_contents_of_returned_list_find_unique_jpi_all_EGAF(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_level_schemes = []
        for r in res:
            levels = e.find_unique_jpi(edata, r)
            try:
                assert isinstance(levels, Iterable)
            except AssertionError:
                assert isinstance(levels, type(None))
            list_level_schemes.append(levels)
            try:
                for level in levels:
                    self.assertIsInstance(levels, list)
                    self.assertIsInstance(level, list)
                    self.assertIsInstance(level[0], int)
                    self.assertIsInstance(level[1], float)
                    self.assertIsInstance(level[2], float)
                    self.assertIsInstance(level[3], int)
                    try:
                        self.assertIsInstance(level[4], float)
                    except AssertionError:
                        self.assertIsNone(level[4])
                    self.assertIsInstance(level[5], int)
            except AssertionError:
                self.assertIsNone(levels)
            except TypeError:
                self.assertIsNone(levels)
        self.assertEqual(len(list_level_schemes), len(res))

    # Find isomers tests:
    def test_find_isomers_in_24Na_returns_list(self):
        isomers_str = e.find_isomers(edata,"Na24",units='best')
        self.assertIsInstance(isomers_str, list)
        isomers_ints = e.find_isomers(edata,11,24,units='best')
        self.assertIsInstance(isomers_ints, list)
        isomers_str = e.find_isomers(edata,"Na24",units='seconds')
        self.assertIsInstance(isomers_str, list)
        isomers_ints = e.find_isomers(edata,11,24,units='seconds')
        self.assertIsInstance(isomers_ints, list)
        isomers_str = e.find_isomers(edata,"Na24",units='s')
        self.assertIsInstance(isomers_str, list)
        isomers_ints = e.find_isomers(edata,11,24,units='s')
        self.assertIsInstance(isomers_ints, list)

    def test_find_isomers_in_24Na_with_wrong_or_no_units_returns_None(self):
        isomers_str = e.find_isomers(edata,"Na24",units='worst')
        self.assertIsNone(isomers_str)
        isomers_ints = e.find_isomers(edata,11,24,units='worst')
        self.assertIsNone(isomers_ints)
        isomers_str = e.find_isomers(edata,"Na24")
        self.assertIsNone(isomers_str)
        isomers_ints = e.find_isomers(edata,11,24)
        self.assertIsNone(isomers_ints)
        isomers_str = e.find_isomers(edata,"Na24",units='ms')
        self.assertIsInstance(isomers_str, type(None))
        isomers_ints = e.find_isomers(edata,11,24,units='ps')
        self.assertIsInstance(isomers_ints, type(None))            

    def test_find_isomers_not_in_EGAF_returns_None(self):
        isomers_str = e.find_isomers(edata,"Kr88")
        self.assertIsNone(isomers_str)

    def test_find_isomers_for_illegal_string_returns_None(self):
        isomers_str = e.find_isomers(edata,"%^&#@Jnsjksd125ASB")
        self.assertIsNone(isomers_str)

    def test_find_isomers_for_wrong_arguments_returns_None(self):
        isomers = e.find_isomers(edata,6,"C13",units='best')
        self.assertIsNone(isomers)
        isomers = e.find_isomers(edata,300.564,1.25673,units='s')
        self.assertIsNone(isomers)
        isomers = e.find_isomers(edata,"O17","La140",8,units='seconds')
        self.assertIsNone(isomers)
        
    def test_find_isomers_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.find_isomers(6, 13, edata, units='best')
        with self.assertRaises(TypeError):
            e.find_isomers("C13", edata, units='best')
        bad_dict = {'a':0, 'b':1, 'c':2}
        with self.assertRaises(TypeError):
            e.find_isomers(bad_dict, edata, units='s')

    def test_find_isomers_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.find_isomers(XXXdata, "Si29", units='seconds')            
            
    def test_find_isomers_raises_KeyError_if_bad_dict_items_in_list(self):   
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.find_isomers(bad_dict_items_in_list, "Si29", units='best')        

    def test_contents_of_returned_list_find_isomers_all_EGAF(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_level_schemes = []
        for r in res:
            levels = e.find_isomers(edata, r, units='best')
            try:
                assert isinstance(levels, Iterable)
            except AssertionError:
                assert isinstance(levels, type(None))
            list_level_schemes.append(levels)
            try:
                for level in levels:
                    self.assertIsInstance(levels, list)
                    self.assertIsInstance(level, list)
                    self.assertIsInstance(level[0], int)
                    self.assertIsInstance(level[1], float)
                    self.assertIsInstance(level[2], float)
                    self.assertIsInstance(level[3], float)
                    self.assertIsInstance(level[4], float)
                    self.assertIsInstance(level[5], str)
            except TypeError:
                self.assertIsNone(levels)
        self.assertEqual(len(list_level_schemes), len(res))
