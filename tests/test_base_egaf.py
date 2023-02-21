import pytest
from collections.abc import Iterable
import pyEGAF as egaf
e = egaf.EGAF()
edata = e.load_egaf()

# Testing methods of the BaseEGAF class:

def test_egaf_datasets_are_in_list_object():
    assert type(e.load_egaf()) is list
def test_load_245_egaf_datasets_into_list():
    assert len(e.load_egaf()) == 245

def test_targets_are_in_list_object():
    assert type(e.egaf_target_list(edata)) is list    
def test_first_target_element_is_Ag107():
    assert e.egaf_target_list(edata)[0] == str('Ag107')
def test_last_target_element_is_Zr96():
    assert e.egaf_target_list(edata)[-1] == str('Zr96')

def test_residuals_are_in_list_object():
    assert type(e.egaf_residual_list(edata)) is list
def test_first_residual_element_is_Ag108():
    assert e.egaf_residual_list(edata)[0] == str('Ag108')
def test_last_residual_element_is_Zr97():
    assert e.egaf_residual_list(edata)[-1] == str('Zr97')        

def test_targets_and_residuals_are_in_dict_object():
    assert type(e.egaf_target_residual_dict(edata)) is dict
def test_heaviest_CN_in_target_residual_dict_is_U239():
    assert [(k,v) for (k,v) in e.egaf_target_residual_dict(edata).items() if k=='U239'] == [('U239', ('U238', 92, 239))]
def test_lightest_CN_in_target_residual_dict_is_H2():
    assert [(k,v) for (k,v) in e.egaf_target_residual_dict(edata).items() if k=='H2'] == [('H2', ('H1', 1, 2))]

def test_get_stats_for_140La():
    assert e.get_stats(edata,"La140") == [187, 102, 289, 207]
    assert e.get_stats(edata,57, 140) == [187, 102, 289, 207]
    assert len(e.get_stats(edata,"La140")) == 4
    
def test_get_stats_for_70Se_not_in_EGAF_returns_None():
    assert e.get_stats(edata,"Se70") == None
    assert e.get_stats(edata, 34, 70) == None

def test_get_stats_for_illegal_string_returns_None():
    assert e.get_stats(edata,"THisIsB@LL@CK$") == None

def test_get_stats_for_wrong_number_arguments_returns_None():
    assert e.get_stats(edata,"C13",6,13) == None
    assert e.get_stats(edata,6,13,300.46) == None    

def test_get_stats_for_missing_arguments_returns_None():
    assert e.get_stats(edata) == None    
    
def test_get_stats_returns_iterable_object_for_nucleus_in_EGAF():
    assert type(e.get_stats(edata,"La140")) is list
    assert isinstance(e.get_stats(edata,"La140"), Iterable)

def test_get_stats_raises_NameError_for_wrong_name_list():
    try:
        e.get_stats(XXXedata, "C13")
    except NameError:
        pass
    else:
        raise

def test_get_stats_raises_TypeError_if_first_argument_not_list():
    try:
        e.get_stats("C13",edata)
        e.get_stats(6,13,edata)
    except TypeError:
        pass
    else:
        raise
    
def test_get_stats_raises_KeyError_for_bad_dict_items():
    bad_dict_items_in_list = [{'x':(3.0,4.0,5.0), 'y':3, 'z':'capture-gamma'}]
    try:
        e.get_stats(bad_dict_items_in_list, "C13")
    except KeyError:
        pass
    else:
        raise

def test_get_stats_all_EGAF_confirm_returned_contents():
    res = e.egaf_residual_list(edata)
    assert len(res) == 245
    list_stats = []
    for r in res:
        stats = e.get_stats(edata, r)
        assert len(stats) == 4
        assert type(stats) is list
        assert isinstance(stats, Iterable)
        list_stats.append(stats)
        assert isinstance(stats[0], int)
        assert isinstance(stats[1], int)
        assert isinstance(stats[2], int)
        assert isinstance(stats[3], int)
    assert len(list_stats) == len(res)
    assert isinstance(list_stats, Iterable)
    
        
# Testing the methods of the Meta class

def test_Si29_has_13_primaries_pass_string():
    assert e.num_primaries(edata,"Si29") == 13
def test_Si29_has_13_primaries_pass_integers():
    assert e.num_primaries(edata,14,29) == 13

def test_Si29_has_33_secondaries_pass_string():
    assert e.num_secondaries(edata,"Si29") == 33
def test_Si29_has_33_secondaries_pass_integers():
    assert e.num_secondaries(edata,14,29) == 33    

def test_Si29_has_46_total_gammas_pass_string():
    assert e.num_gammas(edata,"Si29") == 46
def test_Si29_has_46_total_gammas_pass_integers():
    assert e.num_gammas(edata,14,29) == 46

def test_Si29_has_14_levels_pass_string():
    assert e.num_levels(edata,"Si29") == 14
def test_Si29_has_14_levels_pass_integers():
    assert e.num_levels(edata,14,29) == 14    

def test_nucleus_string_not_in_EGAF_in_num_primaries():
    assert e.num_primaries(edata,"Si42") == None
def test_nucleus_integers_not_in_EGAF_in_num_primaries():
    assert e.num_primaries(edata,14,42) == None    
def test_illegal_string_in_num_primaries():
    assert e.num_primaries(edata,"XDFMLklasdcf7y6786+++!!@") == None

    
# Testing the methods of the Uncertainties class:

# Calculate error on 3.0(1) * 7.0(4)
# >>>e.quad_error(21,3,0.1,7,0.4)
# >>>1.3892443989449803
def test_quad_error_float_maximum_output_precision():
    assert e.quad_error(21.0,3.0,0.1,7.0,0.4) == pytest.approx(1.3892443989449803)
def test_quad_error_float_least_required_precision():
    assert e.quad_error(21.0,3.0,0.1,7.0,0.4) == pytest.approx(1.389244)    
def test_quad_error_float_3dp():
    assert e.quad_error(21.0,3.0,0.1,7.0,0.4) == pytest.approx(1.389, 0.001)

