from .base_egaf import *
from .separation import Separation
from .cross_section import CrossSection
from .decay import Levels, Gammas
from .analysis import Analysis
from .cap_gam import CapGam
from .file_handler import RIPL, JSONFile, ENSDF

class EGAF(ENSDF):
    __doc__="""Class to handle EGAF data sets and methods."""

    def __init__(self):
        BaseEGAF.__init__(self)
        Meta.__init__(self)
        Uncertainties.__init__(self)
        Separation.__init__(self)
        CrossSection.__init__(self)
        Levels.__init__(self)
        Gammas.__init__(self)
        Analysis.__init__(self)
        CapGam.__init__(self)
        RIPL.__init__(self)
        JSONFile.__init__(self)
        ENSDF.__init__(self)

