import pytest
import unittest
import numpy as np
import pandas as pd
from collections.abc import Iterable
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

class GammasTests(unittest.TestCase):

    __doc__ = """Unit tests for the following methods of the Gammas class:

    `get_gammas`;
    `find_all_gammas_feeding_gs`;
    `get_gamma_types`;
    `find_gamma`;
    `get_strongest_gammas`.
    """

    # Gamma array tests
    def test_get_gammas_for_29Si_returns_array(self):
        gammas_str = e.get_gammas(edata, "Si29", intensity="elemental")
        self.assertIsInstance(gammas_str, np.ndarray)
        gammas_str = e.get_gammas(edata, "Si29", intensity="isotopic")
        self.assertIsInstance(gammas_str, np.ndarray)
        gammas_str = e.get_gammas(edata, "Si29", intensity="population")
        self.assertIsInstance(gammas_str, np.ndarray)
        assert isinstance(gammas_str, Iterable)

        gammas_ints = e.get_gammas(edata, 14, 29, intensity="elemental")
        self.assertIsInstance(gammas_ints, np.ndarray)
        gammas_ints = e.get_gammas(edata, 14, 29, intensity="isotopic")
        self.assertIsInstance(gammas_ints, np.ndarray)
        gammas_ints = e.get_gammas(edata, 14, 29, intensity="population")
        self.assertIsInstance(gammas_ints, np.ndarray)
        assert isinstance(gammas_ints, Iterable)

    def test_get_gammas_for_nucleus_not_in_EGAF_returns_None(self):
        gammas_str = e.get_gammas(edata, "Se70", intensity="elemental")
        self.assertIsNone(gammas_str)
        gammas_ints = e.get_gammas(edata, 34, 88, intensity="isotopic")
        self.assertIsNone(gammas_ints)
        gammas_ints = e.get_gammas(edata, 55, 119, intensity="population")
        self.assertIsNone(gammas_ints)

    def test_get_gammas_for_illegal_string_returns_None(self):
        gammas_str = e.get_gammas(edata, "THisIsB@LL@CK$", intensity="isotopic")
        self.assertIsNone(gammas_str)
        gammas_str = e.get_gammas(edata, "Si29", intensity="XXXelemental")
        self.assertIsNone(gammas_str)
        gammas_str = e.get_gammas(edata, "Si29.", intensity="population")
        self.assertIsNone(gammas_str)
        
    def test_get_gammas_for_wrong_arguments_returns_None(self):
        gammas_str = e.get_gammas(edata,14.9,28.9,intensity="isotopic")
        self.assertIsNone(gammas_str)
        gammas_str = e.get_gammas(edata,"14","29","Si29",intensity="isotopic")
        self.assertIsNone(gammas_str)

    def test_get_gammas_for_wrong_or_missing_kwargs_returns_None(self):
        gammas_str = e.get_gammas(edata, 14, 29, intensity="XXXpopulationXX")
        self.assertIsNone(gammas_str)
        gammas_str = e.get_gammas(edata, "Si29")
        self.assertIsNone(gammas_str)

    def test_get_gammas_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.get_gammas("Si29", edata, intensity="elemental")
        with self.assertRaises(TypeError):
            e.get_gammas(14, 29, edata, intensity="elemental")

    def test_get_gammas_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.get_gammas(XXXedata, "Si29", intensity="population")
            
    def test_get_gammas_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.get_gammas(bad_dict_items_in_list, "Si29", intensity="isotopic")

    def test_get_gammas_in_all_EGAF_returned_contents_of_array(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_gamma_arrays = []
        for r in res:
            gammas = e.get_gammas(edata, r, intensity="elemental")
            self.assertIsInstance(gammas, np.ndarray)
            self.assertIsInstance(gammas, Iterable)
            list_gamma_arrays.append(gammas)
            for gamma in gammas:
                self.assertIsInstance(gamma, np.ndarray)
                self.assertIsInstance(gamma, Iterable)
                self.assertIsInstance(gamma[0], float)
                self.assertIsInstance(gamma[1], float)
                self.assertIsInstance(gamma[2], float)
                self.assertIsInstance(gamma[3], float)
                self.assertIsInstance(gamma[4], float)
                self.assertIsInstance(gamma[5], float)
                self.assertIsInstance(gamma[6], float)
                self.assertIsInstance(gamma[7], float)
                self.assertIsInstance(gamma[8], float)
                self.assertIsInstance(gamma[9], float)
        self.assertEqual(len(list_gamma_arrays), len(res))
        self.assertIsInstance(list_gamma_arrays, Iterable)


    # Gammas feeding GS and associated properties
    def test_find_all_gammas_feeding_gs_for_29Si_returns_array(self):
        all_gammas_str = e.find_all_gammas_feeding_gs(edata, "Si29", intensity="isotopic")
        self.assertIsInstance(all_gammas_str, list)
        all_gammas_str = e.find_all_gammas_feeding_gs(edata, "Si29", intensity="population")
        self.assertIsInstance(all_gammas_str, list)
        assert isinstance(all_gammas_str, Iterable)

        all_gammas_ints = e.find_all_gammas_feeding_gs(edata, 14, 29, intensity="isotopic")
        self.assertIsInstance(all_gammas_ints, list)
        all_gammas_ints = e.find_all_gammas_feeding_gs(edata, 14, 29, intensity="population")
        self.assertIsInstance(all_gammas_ints, list)
        assert isinstance(all_gammas_ints, Iterable)

    def test_find_all_gammas_feeding_gs_for_illegal_string_returns_None(self):
        all_gammas_str = e.get_gammas(edata, "THisIsB@LL@CK$", intensity="isotopic")
        self.assertIsNone(all_gammas_str)

    def test_find_all_gammas_feeding_gs_for_wrong_kwargs_returns_None(self):
        all_gammas_str = e.find_all_gammas_feeding_gs(edata, "Si29", intensity="elemental")
        self.assertIsNone(all_gammas_str)
        all_gammas_ints = e.find_all_gammas_feeding_gs(edata, 14, 29, intensity="elemental")
        self.assertIsNone(all_gammas_ints)

    def test_find_all_gammas_feeding_gs_for_missing_kwargs_returns_None(self):
        all_gammas_str = e.find_all_gammas_feeding_gs(edata, "Si29")
        self.assertIsNone(all_gammas_str)
        all_gammas_ints = e.find_all_gammas_feeding_gs(edata, 14, 29)
        self.assertIsNone(all_gammas_ints)

    def test_find_all_gammas_feeding_gs_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.find_all_gammas_feeding_gs(14, 29, edata, intensity="population")
        with self.assertRaises(TypeError):
            e.find_all_gammas_feeding_gs("Si29", edata, intensity="population")

    def test_find_all_gammas_feeding_gs_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.find_all_gammas_feeding_gs(XXXedata, "Si29", intensity="population")
        with self.assertRaises(NameError):
            e.find_all_gammas_feeding_gs(XXXedata, 14, 29, intensity="population")

    def test_find_all_gammas_feeding_gs_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.find_all_gammas_feeding_gs(bad_dict_items_in_list, "Si29", intensity="isotopic")
                
    def test_find_all_gammas_feeding_gs_in_all_EGAF_returned_contents_of_list(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_gammas = []
        for r in res:
            gammas = e.find_all_gammas_feeding_gs(edata, r, intensity="isotopic")
            list_gammas.append(gammas)
            try:
                self.assertIsInstance(gammas, list)
                self.assertIsInstance(gammas, Iterable)
                for gamma in gammas:
                    self.assertIsInstance(gamma, list)
                    self.assertIsInstance(gamma, Iterable)
                    self.assertIsInstance(gamma[0], int)
                    self.assertIsInstance(gamma[1], int)
                    self.assertIsInstance(gamma[2], float)
                    self.assertIsInstance(gamma[3], float)
                    self.assertIsInstance(gamma[4], float)
                    self.assertIsInstance(gamma[5], float)
                    self.assertIsInstance(gamma[6], float)
                    self.assertIsInstance(gamma[7], float)
                    self.assertIsInstance(gamma[8], float)
                    self.assertIsInstance(gamma[9], float)
            except AssertionError:
                self.assertIsNone(gammas)
        self.assertEqual(len(list_gammas), len(res))
        self.assertIsInstance(list_gammas, Iterable)            


    # Gamma types
    def test_get_gamma_types_for_29Si_returns_array(self):
        all_gammas_str = e.get_gamma_types(edata, "Si29",gammas="primary", intensity="isotopic")
        self.assertIsInstance(all_gammas_str, list)
        all_gammas_str = e.get_gamma_types(edata, "Si29",gammas="primary", intensity="population")
        self.assertIsInstance(all_gammas_str, list)
        assert isinstance(all_gammas_str, Iterable)
        all_gammas_str = e.get_gamma_types(edata, "Si29",gammas="secondary", intensity="isotopic")
        self.assertIsInstance(all_gammas_str, list)
        all_gammas_str = e.get_gamma_types(edata, "Si29",gammas="secondary", intensity="population")
        self.assertIsInstance(all_gammas_str, list)
        all_gammas_str = e.get_gamma_types(edata, "Si29",gammas="secondary", intensity="relative")
        self.assertIsInstance(all_gammas_str, list)
        assert isinstance(all_gammas_str, Iterable)

        all_gammas_ints = e.get_gamma_types(edata, 14, 29,gammas="primary", intensity="isotopic")
        self.assertIsInstance(all_gammas_ints, list)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29,gammas="primary", intensity="population")
        self.assertIsInstance(all_gammas_ints, list)
        assert isinstance(all_gammas_ints, Iterable)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29,gammas="secondary", intensity="isotopic")
        self.assertIsInstance(all_gammas_ints, list)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29,gammas="secondary", intensity="population")
        self.assertIsInstance(all_gammas_ints, list)
        assert isinstance(all_gammas_ints, Iterable)

    def test_get_gamma_types_for_illegal_string_returns_None(self):
        all_gammas_str = e.get_gamma_types(edata, "THisIsB@LL@CK$",gammas="primary", intensity="isotopic")
        self.assertIsNone(all_gammas_str)
        all_gammas_str = e.get_gamma_types(edata, "THisIsB@LL@CK$",gammas="secondary", intensity="isotopic")
        self.assertIsNone(all_gammas_str)

    def test_get_gamma_types_for_wrong_kwargs_returns_None(self):
        all_gammas_str = e.get_gamma_types(edata, "Si29",gammas="primary", intensity="isotopes")
        self.assertIsNone(all_gammas_str)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29,gammas="secondary", intensity="elements")
        self.assertIsNone(all_gammas_ints)
        all_gammas_str = e.get_gamma_types(edata, "Si29",gammas="tertiary", intensity="isotopic")
        self.assertIsNone(all_gammas_str)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29,gammas="observed", intensity="population")
        self.assertIsNone(all_gammas_ints)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29,gammas="nogammas", intensity="nointensities")
        self.assertIsNone(all_gammas_ints)

    def test_get_gamma_types_for_missing_kwargs_returns_None(self):
        all_gammas_str = e.get_gamma_types(edata, "Si29")
        self.assertIsNone(all_gammas_str)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29)
        self.assertIsNone(all_gammas_ints)
        all_gammas_str = e.get_gamma_types(edata, "Si29", gammas="primary")
        self.assertIsNone(all_gammas_str)
        all_gammas_ints = e.get_gamma_types(edata, 14, 29, intensity="population")
        self.assertIsNone(all_gammas_ints)

    def test_get_gamma_types_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.get_gamma_types(14, 29, edata,gammas="primary", intensity="population")
        with self.assertRaises(TypeError):
            e.get_gamma_types("Si29", edata,gammas="primary", intensity="elemental")
        with self.assertRaises(TypeError):
            e.get_gamma_types(14, 29, edata,gammas="secondary", intensity="isotopic")
        with self.assertRaises(TypeError):
            e.get_gamma_types("Si29", edata,gammas="secondary", intensity="isotopic")            

    def test_get_gamma_types_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.get_gamma_types(XXXedata, "Si29",gammas="primary", intensity="population")
        with self.assertRaises(NameError):
            e.get_gamma_types(XXXedata, 14, 29,gammas="primary", intensity="population")
        with self.assertRaises(NameError):
            e.get_gamma_types(XXXedata, "Si29",gammas="secondary", intensity="population")
        with self.assertRaises(NameError):
            e.get_gamma_types(XXXedata, 14, 29,gammas="secondary", intensity="population")

    def test_get_gamma_types_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.get_gamma_types(bad_dict_items_in_list, "Si29",gammas="primary", intensity="isotopic")
        with self.assertRaises(KeyError):
            e.get_gamma_types(bad_dict_items_in_list, "Si29",gammas="secondary", intensity="isotopic")

        

    def test_get_gamma_types_in_all_EGAF_returned_contents_of_list(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_gammas = []
        for r in res:
            gammas = e.get_gamma_types(edata, r, gammas='primary', intensity="isotopic")
            #gammas = e.get_gamma_types(edata, r, gammas='secondary', intensity="isotopic")
            list_gammas.append(gammas)
            try:
                self.assertIsInstance(gammas, list)
                self.assertIsInstance(gammas, Iterable)
                for gamma in gammas:
                    self.assertIsInstance(gamma, list)
                    self.assertIsInstance(gamma, Iterable)
                    self.assertIsInstance(gamma[0], int)
                    self.assertIsInstance(gamma[1], int)
                    self.assertIsInstance(gamma[2], float)
                    self.assertIsInstance(gamma[3], float)
                    self.assertIsInstance(gamma[4], float)
                    self.assertIsInstance(gamma[5], float)
                    self.assertIsInstance(gamma[6], float)
                    self.assertIsInstance(gamma[7], float)
                    self.assertIsInstance(gamma[8], float)
                    self.assertIsInstance(gamma[9], float)
                    self.assertIsInstance(gamma[10], str)
            except AssertionError:
                self.assertIsNone(gammas)
        self.assertEqual(len(list_gammas), len(res))
        self.assertIsInstance(list_gammas, Iterable)


    # Find gammas
    def test_find_gamma_at_1273keV_returns_DataFrame(self):
        gamma_str = e.find_gamma(edata, 1273)
        self.assertIsInstance(gamma_str, pd.core.frame.DataFrame)
        gamma_str = e.find_gamma(edata, 1273, intensity="isotopic")
        self.assertIsInstance(gamma_str, pd.core.frame.DataFrame)
        gamma_str = e.find_gamma(edata, 1273, intensity="population")
        self.assertIsInstance(gamma_str, pd.core.frame.DataFrame)
        gamma_str = e.find_gamma(edata, 1273, 2.5, intensity="elemental")
        self.assertIsInstance(gamma_str, pd.core.frame.DataFrame)
        gamma_str = e.find_gamma(edata, 1273, 2.5, intensity="population")
        self.assertIsInstance(gamma_str, pd.core.frame.DataFrame)
        gamma_str = e.find_gamma(edata, 1273, 2.5, intensity="Gobbledygook")
        self.assertIsInstance(gamma_str, pd.core.frame.DataFrame)
        gamma_str = e.find_gamma(edata, 1273, intensity="Gobbledygook")
        self.assertIsInstance(gamma_str, pd.core.frame.DataFrame)
        assert isinstance(gamma_str, Iterable)

    def test_find_gamma_returns_None_if_no_gammas_within_search_window(self):
        gamma_str = e.find_gamma(edata, 10000)
        self.assertIsNone(gamma_str)
        gamma_str = e.find_gamma(edata, 10000, intensity="elemental")
        self.assertIsNone(gamma_str)
        
    def test_find_gamma_at_1273keV_as_string_raises_TypeError(self):        
        with self.assertRaises(TypeError):
            e.find_gamma(edata, "1273")
        with self.assertRaises(TypeError):
            e.find_gamma(edata, "1273", intensity="population")
            
    def test_find_gamma_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.find_gamma("Si29", edata, intensity="elemental")
        with self.assertRaises(TypeError):
            e.find_gamma(1273, edata)            

    def test_find_gamma_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.find_gamma(XXXedata, 1273, intensity="population")
        with self.assertRaises(NameError):
            e.find_gamma(XXXedata, 1273)

    def test_find_gamma_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.find_gamma(bad_dict_items_in_list, 1273, intensity="isotopic")
        with self.assertRaises(KeyError):
            e.find_gamma(bad_dict_items_in_list, 1273, 0.75)
        with self.assertRaises(KeyError):
            e.find_gamma(bad_dict_items_in_list, 1273, 0.75, intensity="elemental")            


    # Strongest gammas
    def test_get_strongest_gammas_for_29Si_returns_array(self):
        all_gammas_str = e.get_strongest_gammas(edata, "Si29", intensity="isotopic")
        self.assertIsInstance(all_gammas_str, pd.core.frame.DataFrame)
        all_gammas_str = e.get_strongest_gammas(edata, "Si29", intensity="elemental")
        self.assertIsInstance(all_gammas_str, pd.core.frame.DataFrame)
        all_gammas_str = e.get_strongest_gammas(edata, 14, 29, intensity="population")
        self.assertIsInstance(all_gammas_str, pd.core.frame.DataFrame)
        all_gammas_str = e.get_strongest_gammas(edata, 14, 29, intensity="relative")
        self.assertIsInstance(all_gammas_str, pd.core.frame.DataFrame)
        all_gammas_str = e.get_strongest_gammas(edata, 14, 29)
        self.assertIsInstance(all_gammas_str, pd.core.frame.DataFrame)
        all_gammas_str = e.get_strongest_gammas(edata, "Si29")
        self.assertIsInstance(all_gammas_str, pd.core.frame.DataFrame)
        self.assertIsInstance(all_gammas_str, Iterable)

    def test_get_strongest_gammas_for_illegal_string_returns_None(self):
        all_gammas_str = e.get_strongest_gammas(edata, "THisIsB@LL@CK$", intensity="isotopic")
        self.assertIsNone(all_gammas_str)

    def test_get_strongest_gammas_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.get_strongest_gammas(14, 29, edata, intensity="population")
        with self.assertRaises(TypeError):
            e.get_strongest_gammas("Si29", edata, intensity="elemental")
        with self.assertRaises(TypeError):
            e.get_strongest_gammas(14, 29, edata, intensity="isotopic")
        with self.assertRaises(TypeError):
            e.get_strongest_gammas("Si29", edata, intensity="relative")

    def test_get_strongest_gammas_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.get_strongest_gammas(XXXedata, "Si29", intensity="population")
        with self.assertRaises(NameError):
            e.get_strongest_gammas(XXXedata, 14, 29, intensity="relative")

    def test_get_strongest_gammas_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.get_strongest_gammas(bad_dict_items_in_list, "Si29", intensity="isotopic")
        with self.assertRaises(KeyError):
            e.get_strongest_gammas(bad_dict_items_in_list, 14, 29, intensity="elemental")

