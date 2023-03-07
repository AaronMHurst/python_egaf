from .pyEGAF import *
from .base_egaf import *
from .separation import Separation
from .cross_section import CrossSection
from .decay import Levels, Gammas
from .analysis import Analysis
from .cap_gam import CapGam
from .file_handler import RIPL, JSONFile, ENSDF

__version__='0.1.0'
__author__='Aaron M. Hurst'

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    """Function to return absolute path of the data files inside the root of 
    the Python package."""
    return os.path.join(_ROOT, path)
