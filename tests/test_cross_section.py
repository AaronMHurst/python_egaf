import pytest
import unittest
import numpy as np
from collections.abc import Iterable
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

# Testing methods of the CrossSection class:

# Individual total cross section functions
def test_total_cs_Si28ng_is_tuple():
    assert type(e.get_total_cross_section(edata,"Si28")) is tuple
    assert type(e.get_total_cross_section(edata,14,28)) is tuple

def test_total_cs_for_nuclides_not_in_EGAF_is_None():
    # Se70(n,g)
    assert e.get_total_cross_section(edata,"Se70") == None
    # 88Kr(n,g)
    assert e.get_total_cross_section(edata,36,88) == None
    
def test_total_cs_for_Si28_ng():
    cs = e.get_total_cross_section(edata,"Si28")[0]
    d_cs = e.get_total_cross_section(edata,"Si28")[1]
    unit = e.get_total_cross_section(edata,"Si28")[2]
    ref = e.get_total_cross_section(edata,"Si28")[3]
    
    assert isinstance(cs, float)
    assert isinstance(d_cs, float)
    assert isinstance(unit, str)
    assert isinstance(ref, str)

    assert cs == pytest.approx(0.177, d_cs)
    assert unit == 'b'
    assert ref == '1981MuZQ'

# Individual abundance functions
def test_abundance_Cl35_is_2_element_tuple():
    assert type(e.get_abundance(edata,"Cl35")) is tuple
    assert len(e.get_abundance(edata,"Cl35")) == 2
    assert type(e.get_abundance(edata,17,35)) is tuple
    assert len(e.get_abundance(edata,17,35)) == 2

def test_abundance_for_nuclides_not_in_EGAF_is_None():
    # Se70(n,g)
    assert e.get_abundance(edata,"Se70") == None
    # 88Kr(n,g)
    assert e.get_abundance(edata,36,88) == None    

def test_abundance_for_Cl35():
    abund = e.get_abundance(edata,"Cl35")[0]
    d_abund = e.get_abundance(edata,"Cl35")[1]
    
    assert isinstance(abund, float)
    assert isinstance(d_abund, float)
    assert abund == pytest.approx(75.8, d_abund)

def test_total_abundance_chlorine_Cl35_plus_Cl37_is_100_precent():
    abund = e.get_abundance(edata,"Cl35")[0]
    abund += e.get_abundance(edata,"Cl37")[0]
    d_abund = e.get_abundance(edata,"Cl35")[1]**2
    d_abund += e.get_abundance(edata,"Cl37")[1]**2
    d_abund = np.sqrt(d_abund)
    
    assert isinstance(abund, float)
    assert isinstance(d_abund, float)
    assert abund == pytest.approx(100, d_abund)

# All total cross section functions
def test_all_total_cross_section_dict_data_types():
    assert type(e.get_all_total_cross_sections(edata)) is dict
    assert len(e.get_all_total_cross_sections(edata)) == 245
    assert isinstance(e.get_all_total_cross_sections(edata), Iterable)

    all_tot_cs = e.get_all_total_cross_sections(edata)
    for i,(k,v) in enumerate(all_tot_cs.items()):
        assert isinstance(k, str)
        assert isinstance(v[0], str)
        assert isinstance(v[1], float)
        assert isinstance(v[2], float)
        assert isinstance(v[3], str)
        assert isinstance(v[4], str)
        assert type(v) is tuple

        assert v[3] == 'b'
        assert v[4].upper() == ('1981MUZQ') or v[4].upper() == ('1984MUZY') or v[4].upper() == ('2018MUZY')

        if i == 0:
            assert k == 'Ag107'
        if i == 244:
            assert k == 'Zr96'

# All total cross section functions
def test_all_abundances_dict_data_types():
    assert type(e.get_all_abundances(edata)) is dict
    assert len(e.get_all_abundances(edata)) == 245
    assert isinstance(e.get_all_abundances(edata), Iterable)

    all_abunds = e.get_all_abundances(edata)
    for i,(k,v) in enumerate(all_abunds.items()):
        assert isinstance(k, str)
        assert isinstance(v[0], float)
        assert isinstance(v[1], float)
        assert type(v) is tuple

        if i == 1:
            assert k == 'Ag109'
        if i == 243:
            assert k == 'Zr94'            

# Test exception errors thrown by cross-section methods
class CrossSectionTests(unittest.TestCase):

    __doc__ = """Unit tests for the following methods of the CrossSection class:
    
    `get_total_cross_section`;
    `get_abundance`;
    `get_all_total_cross_sections`;
    `get_all_abundances`.
    """
    
    # Throw TypeError exception if no parameters are passed
    def test_total_cross_section_throws_TypeError_without_parameters(self):
        self.assertRaises(TypeError,e.get_total_cross_section)
    def test_abundance_throws_TypeError_without_parameters(self):
        self.assertRaises(TypeError,e.get_abundance)
        
    # Throw TypeError exception if no list gets passed to get_all methods
    def test_all_total_cross_section_throws_TypeError_without_list(self):
        self.assertRaises(TypeError,e.get_all_total_cross_sections)
    def test_all_abundances_TypeError_without_list(self):
        self.assertRaises(TypeError,e.get_all_abundances)
    # Throw TypeError exception if extra positional arguments are given
    def test_all_total_cross_section_throws_TypeError_with_extra_args(self):
        self.assertRaises(TypeError,e.get_all_total_cross_sections,6,12.0)
    def test_all_abundances_TypeError_with_extra_args(self):
        self.assertRaises(TypeError,e.get_all_abundances,"La139")    
    

    # Throw TypeError if list of wrong types gets passed:
    def test_total_cross_section_throws_TypeError_with_list_of_wrong_types(self):
        list_of_ints = [1,2,3,4,5]
        list_of_strings = ['a','b','c','d','e']
        list_of_floats = [6.0,7.0,8.0,9.0,10.0]
        with self.assertRaises(TypeError):
            e.get_total_cross_section(list_of_ints, "Si28")
            e.get_total_cross_section(list_of_strings, "Al27")
            e.get_total_cross_section(list_of_floats, "C12")

    def test_abundance_throws_TypeError_with_list_of_wrong_types(self):
        list_of_ints = [1,2,3,4,5]
        list_of_strings = ['a','b','c','d','e']
        list_of_floats = [6.0,7.0,8.0,9.0,10.0]
        with self.assertRaises(TypeError):
            e.get_abundance(list_of_ints, "Si28")
            e.get_abundance(list_of_strings, "Al27")
            e.get_abundance(list_of_floats, "C12")

    # Throw KeyError for list of bad dictionary items (i.e. not an expected
    # JSON structure):
    def test_cross_section_throws_KeyError_with_bad_dict_in_list(self):
        bad_dict_items_in_list = [{'x':(3.0,4.0,5.0), 'y':3, 'z':'capture-gamma'}]
        with self.assertRaises(KeyError):
            e.get_total_cross_section(bad_dict_items_in_list, "Si28")

    def test_abundance_throws_KeyError_with_bad_dict_in_list(self):
        bad_dict_items_in_list = [{'x':(3.0,4.0,5.0), 'y':3, 'z':'capture-gamma'}]
        with self.assertRaises(KeyError):
            e.get_abundance(bad_dict_items_in_list, "Si28")

    # Throw a NameError exception if incorrect name is provided for a data list
    def test_cross_section_throws_NameError_for_wrong_data_list(self):
        with self.assertRaises(NameError):
            e.get_total_cross_section(XdataX, "C12")

    def test_abundance_throws_NameError_for_wrong_data_list(self):
        with self.assertRaises(NameError):
            e.get_abundance(this_data_list_does_not_exist, "W186")            
            
    def test_all_total_cross_section_throws_NameError_for_wrong_data_list(self):
        with self.assertRaises(NameError):
            e.get_all_total_cross_sections(XXXdataXXX)

    def test_all_abundances_throws_NameError_for_wrong_data_list(self):
        with self.assertRaises(NameError):
            e.get_all_abundances(this_data_list_does_not_exist)
