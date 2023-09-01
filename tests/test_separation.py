import pytest
import unittest
import numpy as np
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

# Testing methods of the Separation class:

# Sn (AME2020) methods:
def test_residual_Sn_AME_returns_tuple():
    assert type(e.get_residual_Sn_AME(edata,"Al28")) is tuple
def test_residual_Sn_AME_tuples_contain_floats():
    assert isinstance(e.get_residual_Sn_AME(edata,"Al28")[0], float)
    assert isinstance(e.get_residual_Sn_AME(edata,"Al28")[1], float)

def test_residual_Sn_AME_28Al_is_approx_residual_Sn_EGAF():
    assert e.get_residual_Sn_AME(edata,"Al28")[0] == pytest.approx(e.get_residual_Sn_EGAF(edata,"Al28")[0], 0.005)
def test_residual_Sn_AME_28Al_is_consistent_residual_Sn_EGAF():
    delta_Sn = np.sqrt(e.get_residual_Sn_AME(edata,"Al28")[1]**2 + e.get_residual_Sn_EGAF(edata,"Al28")[1]**2)
    assert e.get_residual_Sn_AME(edata,"Al28")[0] == pytest.approx(e.get_residual_Sn_EGAF(edata,"Al28")[0], delta_Sn)

# Sp (AME2020) methods:
def test_residual_Sp_AME_returns_tuple():
    assert type(e.get_residual_Sp_AME(edata,14,29)) is tuple
def test_residual_Sp_AME_tuples_contain_floats():
    assert isinstance(e.get_residual_Sp_AME(edata,14,29)[0], float)
    assert isinstance(e.get_residual_Sp_AME(edata,14,29)[1], float)
def test_residual_Sp_AME_29Si_value():
    delta_Sp = e.get_residual_Sp_AME(edata,14,29)[1]
    assert e.get_residual_Sp_AME(edata,14,29)[0] == pytest.approx(12333.33, delta_Sp)
    
# Sn (EGAF) methods:
def test_residual_Sn_EGAF_returns_tuple():
    assert type(e.get_residual_Sn_EGAF(edata,6,13)) is tuple
def test_residual_Sn_EGAF_tuples_contain_floats():
    assert isinstance(e.get_residual_Sn_EGAF(edata,6,13)[0], float)
    assert isinstance(e.get_residual_Sn_EGAF(edata,6,13)[1], float)
def test_residual_Sn_EGAF_13C_value():
    delta_Sn = e.get_residual_Sn_EGAF(edata,6,13)[1]
    assert e.get_residual_Sn_EGAF(edata,6,13)[0] == pytest.approx(4946.309, delta_Sn)
def test_residual_Sn_not_in_EGAF_returns_None():
    assert e.get_residual_Sn_EGAF(edata,34,70) == None
    assert e.get_residual_Sn_EGAF(edata,"Kr88") == None
def test_residual_Sn_EGAF_wrong_input_sequence_returns_None():
    assert e.get_residual_Sn_EGAF(edata,13,6) == None
    assert e.get_residual_Sn_EGAF(edata,6,13,100) == None
    assert e.get_residual_Sn_EGAF(edata,6,"C13") == None
    assert e.get_residual_Sn_EGAF(edata,"C13",13) == None
    assert e.get_residual_Sn_EGAF(edata) == None
    assert e.get_residual_Sn_EGAF("C13") == None
def test_residual_Sn_EGAF_wrong_list_name_returns_None():
    assert e.get_residual_Sn_EGAF(edata,"mcsnOWwq3&*^#$XCJ1285867Ksadk") == None


# All separation-energy methods
def test_all_EGAF_Sn_are_floats_or_None():
    egaf_Sn = e.get_all_separation_energies(edata,"egaf")
    for (k,v) in egaf_Sn.items():
        if v[0] == None:
            assert v[1] == None            
        else:
            assert isinstance(v[0], float)
            assert isinstance(v[1], float)
        assert isinstance(k, str)
    assert len(egaf_Sn) == 245

def test_all_AME2020_Sn_are_floats():
    AME_Sn = e.get_all_separation_energies(edata,"neutron")
    for (k,v) in AME_Sn.items():
        assert isinstance(v[0], float)
        assert isinstance(v[1], float)
        assert isinstance(k, str)
    assert len(AME_Sn) == 245

def test_all_AME2020_Sp_are_floats():
    AME_Sp = e.get_all_separation_energies(edata,"proton")
    for (k,v) in AME_Sp.items():
        assert isinstance(v[0], float)
        assert isinstance(v[1], float)
        assert isinstance(k, str)
    assert len(AME_Sp) == 245


# Test exceptions thrown by separation-energy functions
class SeparationEnergyTests(unittest.TestCase):

    __doc__ = """Unit tests for the following methods of the Separation class:

    `get_residual_Sn_AME`;
    `get_residual_Sp_AME`;
    `get_residual_Sn_EGAF`.
    """
    
    # Throw TypeError exception if no parameters are passed
    def test_residual_Sn_EGAF_throws_TypeError_without_parameters(self):
        self.assertRaises(TypeError,e.get_residual_Sn_EGAF)
    def test_residual_Sn_AME_throws_TypeError_without_parameters(self):
        self.assertRaises(TypeError,e.get_residual_Sn_AME)
    def test_residual_Sp_AME_throws_TypeError_without_parameters(self):
        self.assertRaises(TypeError,e.get_residual_Sp_AME)                

    # Throw TypeError if the wrong list gets passed:
    def test_residual_Sn_EGAF_throws_TypeError_with_list_of_wrong_types(self):
        list_of_ints = [1,2,3,4,5]
        self.assertRaises(TypeError,e.get_residual_Sn_EGAF,list_of_ints,"C13")

    def test_residual_Sn_AME_throws_TypeError_with_list_of_wrong_types(self):
        list_of_strings = ['a','b','c','d','e']
        with self.assertRaises(TypeError):
            e.get_residual_Sn_AME(list_of_strings,"W187")

    def test_residual_Sp_AME_throws_TypeError_with_list_of_wrong_types(self):
        list_of_floats = [6.0,7.0,8.0,9.0,10.0]
        with self.assertRaises(TypeError):
            e.get_residual_Sn_AME(list_of_floats,"Si29")
            
    # Throw KeyError for list of bad dictionary items (i.e. not an expected
    # JSON structure):
    def test_residual_Sn_EGAF_throws_KeyError_with_bad_dict_in_list(self):
        bad_dict_items_in_list = [{'x':(3.0,4.0,5.0), 'y':3, 'z':'capture-gamma'}]
        with self.assertRaises(KeyError):
            e.get_residual_Sn_EGAF(bad_dict_items_in_list,"C13")
            e.get_residual_Sn_AME(bad_dict_items_in_list,"W187")
            e.get_residual_Sp_AME(bad_dict_items_in_list,"Si29")
            
    # Raise a NameError exception if the wrong data list name gets passed
    def test_residual_separation_funcs_throw_NameError_for_wrong_data_list(self):
        try:
            #self.assertRaises(NameError,e.get_residual_Sn_EGAF,data,"13C")
            e.get_residual_Sn_EGAF(XXXdataXX,"C13")
            e.get_residual_Sn_AME(data,"W187")
            e.get_residual_Sp_AME(my_data_list_doesnt_exist,"Si29")
        except NameError:
            pass
        else:
            raise
        
