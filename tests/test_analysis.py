import pytest
import unittest
import numpy as np
import pandas as pd
from collections.abc import Iterable
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

class AnalysisTests(unittest.TestCase):

    # Internal-conversion intensity tests
    def test_intensity_conversion_returns_tuple(self):
        total_trans = e.intensity_conversion(0.00075, 0.00012, 0.000335, 0.000005)
        self.assertIsInstance(total_trans, tuple)
        self.assertIsInstance(total_trans, Iterable)

    def test_intensity_conversion_raises_TypeError_if_too_many_or_too_few_args(self):
        # Too few args
        with self.assertRaises(TypeError):
            e.intensity_conversion(0.00075, 0.00012)
        # Too many args
        with self.assertRaises(TypeError):
            e.intensity_conversion(0.00075, 0.00012, 0.000335, 0.000005, 9.999999)

    def test_intensity_conversion_returned_floats_contents_of_tuple(self):
        total_trans = e.modeled_sigma0(0.00075, 0.00012, 0.000335, 0.000005)
        self.assertIsInstance(total_trans[0], float)
        self.assertIsInstance(total_trans[1], float)

    # Tests for sigma0 determined using experimental cs and P0 from model
    def test_modeled_sigma0_returns_tuple(self):
        s0 = e.modeled_sigma0(0.187, 0.0033, 0.02314, 0.00080)
        self.assertIsInstance(s0, tuple)
        self.assertIsInstance(s0, Iterable)

    def test_modeled_sigma0_raises_TypeError_if_too_many_or_too_few_args(self):
        # Too few args
        with self.assertRaises(TypeError):
            e.modeled_sigma0(0.187, 0.0033)
        # Too many args
        with self.assertRaises(TypeError):
            e.modeled_sigma0(0.187, 0.0033, 0.02314, 0.00080, 0.9)

    def test_modeled_sigma0_returned_floats_contents_of_tuple(self):
        s0 = e.modeled_sigma0(0.187, 0.0033, 0.02314, 0.00080)
        self.assertIsInstance(s0[0], float)
        self.assertIsInstance(s0[1], float)

    # Tests for sigma0 as a function of Ecrit
    def test_modeled_sigma0_ecrit_returns_list(self):
        s0_Ec = e.modeled_sigma0_ecrit(edata, 11, 0.02256, 0.00064, "Si29")
        self.assertIsInstance(s0_Ec, list)
        s0_Ec = e.modeled_sigma0_ecrit(edata, 11, 0.02256, 0.00064, 14, 29)
        self.assertIsInstance(s0_Ec, list)
        
    def test_modeled_sigma0_ecrit_returns_None_if_Ecrit_too_high_or_too_low(self):
        # Too low
        s0_Ec = e.modeled_sigma0_ecrit(edata, 0, 0.02256, 0.00064, "Si29")
        self.assertIsNone(s0_Ec)
        s0_Ec = e.modeled_sigma0_ecrit(edata, 0, 0.02256, 0.00064, 14, 29)
        self.assertIsNone(s0_Ec)

        # Too high
        s0_Ec = e.modeled_sigma0_ecrit(edata, 20, 0.02256, 0.00064, "Si29")
        self.assertIsNone(s0_Ec)
        s0_Ec = e.modeled_sigma0_ecrit(edata, 20, 0.02256, 0.00064, 14, 29)
        self.assertIsNone(s0_Ec)

    def test_modeled_sigma0_ecrit_returns_None_if_residual_not_in_EGAF(self):
        s0_Ec = e.modeled_sigma0_ecrit(edata, 0, 0.02256, 0.00064, "Si42")
        self.assertIsNone(s0_Ec)
        s0_Ec = e.modeled_sigma0_ecrit(edata, 0, 0.02256, 0.00064, 14, 42)
        self.assertIsNone(s0_Ec)

    def test_modeled_sigma0_ecrit_returns_None_for_illegal_string(self):
        s0_Ec = e.modeled_sigma0_ecrit(edata, 0, 0.02256, 0.00064, "THisIsB@LL@CK$")
        self.assertIsNone(s0_Ec)

    def test_modeled_sigma0_ecrit_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.modeled_sigma0_ecrit(11, edata, 0.02256, 0.00064, "Si29")
        with self.assertRaises(TypeError):
            e.modeled_sigma0_ecrit(11, edata, 0.02256, 0.00064, 14, 29)            

    def test_modeled_sigma0_ecrit_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.modeled_sigma0_ecrit(XXXedata, 11, 0.02256, 0.00064, "Si29")
        with self.assertRaises(NameError):
            e.modeled_sigma0_ecrit(XXXedata, 11, 0.02256, 0.00064, 14, 29)                        

    def test_modeled_sigma0_ecrit_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.modeled_sigma0_ecrit(bad_dict_items_in_list, 11, 0.02256, 0.00064, "Si29")
        with self.assertRaises(KeyError):
            e.modeled_sigma0_ecrit(bad_dict_items_in_list, 11, 0.02256, 0.00064, 14, 29)            

    def test_modeled_sigma0_returned_floats_contents_of_tuple(self):
        s0_Ec = e.modeled_sigma0_ecrit(edata, 11, 0.02256, 0.00064, "Si29")
        self.assertEqual(len(s0_Ec), 5)
        self.assertIsInstance(s0_Ec[0], float)
        self.assertIsInstance(s0_Ec[1], float)
        self.assertIsInstance(s0_Ec[2], float)
        self.assertIsInstance(s0_Ec[3], float)
        self.assertIsInstance(s0_Ec[4], float)

        s0_Ec = e.modeled_sigma0_ecrit(edata, 11, 0.02256, 0.00064, 14, 29)
        self.assertEqual(len(s0_Ec), 5)
        self.assertIsInstance(s0_Ec[0], float)
        self.assertIsInstance(s0_Ec[1], float)
        self.assertIsInstance(s0_Ec[2], float)
        self.assertIsInstance(s0_Ec[3], float)
        self.assertIsInstance(s0_Ec[4], float)
            
            
    # GS-feeding summation tests
    def test_sum_feeding_gs_for_29Si_returns_array(self):
        sf_str = e.sum_feeding_gs(edata, "Si29", intensity="isotopic")
        self.assertIsInstance(sf_str, tuple)
        sf_str = e.sum_feeding_gs(edata, "Si29", intensity="population")
        self.assertIsInstance(sf_str, tuple)
        assert isinstance(sf_str, Iterable)

        sf_ints = e.sum_feeding_gs(edata, 14, 29, intensity="isotopic")
        self.assertIsInstance(sf_ints, tuple)
        sf_ints = e.sum_feeding_gs(edata, 14, 29, intensity="population")
        self.assertIsInstance(sf_ints, tuple)
        assert isinstance(sf_ints, Iterable)

    def test_sum_feeding_gs_for_illegal_string_returns_None(self):
        sf_str = e.get_gammas(edata, "THisIsB@LL@CK$", intensity="isotopic")
        self.assertIsNone(sf_str)

    def test_sum_feeding_gs_for_wrong_kwargs_returns_None(self):
        sf_str = e.sum_feeding_gs(edata, "Si29", intensity="elemental")
        self.assertIsNone(sf_str)
        sf_ints = e.sum_feeding_gs(edata, 14, 29, intensity="elemental")
        self.assertIsNone(sf_ints)

    def test_sum_feeding_gs_for_missing_kwargs_returns_None(self):
        sf_str = e.sum_feeding_gs(edata, "Si29")
        self.assertIsNone(sf_str)
        sf_ints = e.sum_feeding_gs(edata, 14, 29)
        self.assertIsNone(sf_ints)

    def test_sum_feeding_gs_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.sum_feeding_gs(14, 29, edata, intensity="elemental")
        with self.assertRaises(TypeError):
            e.sum_feeding_gs("Si29", edata, intensity="elemental")

    def test_sum_feeding_gs_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.sum_feeding_gs(XXXedata, "Si29", intensity="population")
        with self.assertRaises(NameError):
            e.sum_feeding_gs(XXXedata, 14, 29, intensity="population")

    def test_sum_feeding_gs_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.sum_feeding_gs(bad_dict_items_in_list, "Si29", intensity="isotopic")
            
    def test_sum_feeding_gs_in_all_EGAF_returned_contents_of_tuple(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_sum_feeding_gs = []
        for r in res:
            sf = e.sum_feeding_gs(edata, r, intensity="isotopic")
            list_sum_feeding_gs.append(sf)
            try:
                self.assertIsInstance(sf, tuple)
                self.assertIsInstance(sf, Iterable)
                self.assertIsInstance(sf[0], float)
                self.assertIsInstance(sf[1], float)
            except AssertionError:
                self.assertIsNone(sf)
        self.assertEqual(len(list_sum_feeding_gs), len(res))
        self.assertIsInstance(list_sum_feeding_gs, Iterable)


    # Primary-gamma summation tests
    def test_sum_primaries_for_29Si_returns_array(self):
        sp_str = e.sum_primaries(edata, "Si29", intensity="isotopic")
        self.assertIsInstance(sp_str, tuple)
        sp_str = e.sum_primaries(edata, "Si29", intensity="population")
        self.assertIsInstance(sp_str, tuple)
        assert isinstance(sp_str, Iterable)

        sp_ints = e.sum_primaries(edata, 14, 29, intensity="isotopic")
        self.assertIsInstance(sp_ints, tuple)
        sp_ints = e.sum_primaries(edata, 14, 29, intensity="population")
        self.assertIsInstance(sp_ints, tuple)
        assert isinstance(sp_ints, Iterable)

    def test_sum_primaries_for_illegal_string_returns_None(self):
        sp_str = e.get_gammas(edata, "THisIsB@LL@CK$", intensity="isotopic")
        self.assertIsNone(sp_str)

    def test_sum_primaries_for_wrong_kwargs_returns_None(self):
        sp_str = e.sum_primaries(edata, "Si29", intensity="elemental")
        self.assertIsNone(sp_str)
        sp_ints = e.sum_primaries(edata, 14, 29, intensity="elemental")
        self.assertIsNone(sp_ints)

    def test_sum_primaries_for_missing_kwargs_returns_None(self):
        sp_str = e.sum_primaries(edata, "Si29")
        self.assertIsNone(sp_str)
        sp_ints = e.sum_primaries(edata, 14, 29)
        self.assertIsNone(sp_ints)

    def test_sum_primaries_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            e.sum_primaries(14, 29, edata, intensity="elemental")
        with self.assertRaises(TypeError):
            e.sum_primaries("Si29", edata, intensity="elemental")

    def test_sum_primaries_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            e.sum_primaries(XXXedata, "Si29", intensity="population")
        with self.assertRaises(NameError):
            e.sum_primaries(XXXedata, 14, 29, intensity="population")

    def test_sum_primaries_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.sum_primaries(bad_dict_items_in_list, "Si29", intensity="isotopic")
            
    def test_sum_primaries_in_all_EGAF_returned_contents_of_tuple(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_sum_primaries = []
        for r in res:
            sp = e.sum_primaries(edata, r, intensity="isotopic")
            list_sum_primaries.append(sp)
            try:
                self.assertIsInstance(sp, tuple)
                self.assertIsInstance(sp, Iterable)
                self.assertIsInstance(sp[0], float)
                self.assertIsInstance(sp[1], float)
            except AssertionError:
                self.assertIsNone(sp)
        self.assertEqual(len(list_sum_primaries), len(res))
        self.assertIsInstance(list_sum_primaries, Iterable)

    # Tests for normalising intensities using P0 from model
    def test_normalise_intensities_returns_list(self):
        norm = e.normalise_intensities(edata, "Si29", 0.02217, 0.00051)
        self.assertIsInstance(norm, list)
        norm = e.normalise_intensities(edata, 14, 29, 0.02217, 0.00051)
        self.assertIsInstance(norm, list)

    def test_normalise_intensities_returns_None_if_residual_not_in_EGAF(self):
        norm = e.normalise_intensities(edata, "Si55", 0.02217, 0.00051)
        self.assertIsNone(norm)
        norm = e.normalise_intensities(edata, 14, 55, 0.02217, 0.00051)
        self.assertIsNone(norm)

    def test_normalise_intensities_returns_None_for_illegal_string(self):
        norm = e.normalise_intensities(edata, "THisIsB@LL@CK$", 0.02217, 0.00051)
        self.assertIsNone(norm)

    def test_normalise_intensities_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            norm = e.normalise_intensities("Si29", edata, 0.02217, 0.00051)
        with self.assertRaises(TypeError):
            norm = e.normalise_intensities(14, 29, edata, 0.02217, 0.00051)
        
    def test_normalise_intensities_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            norm = e.normalise_intensities(XXXedataX, "Si29", 0.02217, 0.00051)
        with self.assertRaises(NameError):
            norm = e.normalise_intensities(XXXedataX, 14, 29, 0.02217, 0.00051)

    def test_normalise_intensities_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.normalise_intensities(bad_dict_items_in_list, "Si29", 0.02217, 0.00051)
        with self.assertRaises(KeyError):
            e.normalise_intensities(bad_dict_items_in_list, 14, 29, 0.02217, 0.00051)

    def test_normalise_intensities_returned_contents_of_list(self):
        norm = e.normalise_intensities(edata, "Si29", 0.02217, 0.00051)
        for n in norm:
            self.assertEqual(len(n), 7)
            self.assertIsInstance(n[0], int)
            self.assertIsInstance(n[1], float)
            self.assertIsInstance(n[2], float)
            try:
                self.assertIsInstance(n[3], float)
            except AssertionError:
                self.assertIsInstance(n[3], int)
            self.assertIsInstance(n[4], float)
            self.assertIsInstance(n[5], float)
            self.assertIsInstance(n[6], float)


    # Tests for generating level-depopulation data
    def test_level_depopulations_returns_list(self):
        levels = e.level_depopulations(edata, "Si29")
        self.assertIsInstance(levels, list)
        levels = e.level_depopulations(edata, 14, 29)
        self.assertIsInstance(levels, list)

    def test_level_depopulations_returns_None_if_residual_not_in_EGAF(self):
        levels = e.level_depopulations(edata, "Si55")
        self.assertIsNone(levels)
        levels = e.level_depopulations(edata, 14, 55)
        self.assertIsNone(levels)

    def test_level_depopulations_returns_None_for_illegal_string(self):
        levels = e.level_depopulations(edata, "THisIsB@LL@CK$")
        self.assertIsNone(levels)

    def test_level_depopulations_raises_TypeError_if_first_arg_not_list(self):
        with self.assertRaises(TypeError):
            levels = e.level_depopulations("Si29", edata)
        with self.assertRaises(TypeError):
            levels = e.level_depopulations(14, 29, edata)
        
    def test_level_depopulations_raises_NameError_if_wrong_list_name(self):
        with self.assertRaises(NameError):
            levels = e.level_depopulations(XXXedataX, "Si29")
        with self.assertRaises(NameError):
            levels = e.level_depopulations(XXXedataX, 14, 29)

    def test_level_depopulations_raises_KeyError_if_bad_dict_items_in_list(self):
        bad_dict_items_in_list = [{'a':0, 'b':1, 'c':2}]
        with self.assertRaises(KeyError):
            e.level_depopulations(bad_dict_items_in_list, "Si29")
        with self.assertRaises(KeyError):
            e.level_depopulations(bad_dict_items_in_list, 14, 29)

    def test_level_depopulations_returned_contents_of_list(self):
        levels = e.level_depopulations(edata, "Si29")
        for l in levels:
            self.assertEqual(len(l), 5)
            self.assertIsInstance(l[0], int)
            self.assertIsInstance(l[1], float)
            self.assertIsInstance(l[2], float)
            try:
                self.assertIsInstance(l[3], float)
            except AssertionError:
                self.assertIsInstance(l[3], int)
            self.assertIsInstance(l[4], float)

    def test_level_depopulations_in_all_EGAF_returned_contents_of_list(self):
        res = e.egaf_residual_list(edata)
        self.assertEqual(len(res), 245)
        list_levels = []
        for r in res:
            levels = e.level_depopulations(edata, r)
            list_levels.append(levels)
            self.assertIsInstance(levels, list)
            self.assertIsInstance(levels, Iterable)
            for l in levels:
                self.assertEqual(len(l), 5)
                self.assertIsInstance(l[0], int)
                self.assertIsInstance(l[1], float)
                self.assertIsInstance(l[2], float)
                try:
                    self.assertIsInstance(l[3], float)
                except AssertionError:
                    self.assertIsInstance(l[3], int)
                self.assertIsInstance(l[4], float)

        self.assertEqual(len(list_levels), len(res))
        self.assertIsInstance(list_levels, list)
        self.assertIsInstance(list_levels, Iterable)

