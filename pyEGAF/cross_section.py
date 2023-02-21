from .base_egaf import *
from .separation import Separation

class CrossSection(Separation):
    __doc__="""Class to handle adopted total thermal-neutron capture 
    cross section data and natural abundances"""

    def __init__(self):
        BaseEGAF.__init__(self)
        Meta.__init__(self)
        Uncertainties.__init__(self)
        Separation.__init__(self)

    def get_total_cross_section(self,list,*args):
        """Adopted total radiative thermal neutron-capture cross sections for a 
        measured target from the EGAF database.  All values are taken from 
        various editions of the "Atlas of Neutron Resonances":
        
        [1981MuZQ]: S.F. Mughabghab, M.Divadeenam, N.E.Holden, "Neutron Cross 
        Sections, Neutron Resonance Parameters and Thermal Cross 
        Sections, Part A, Z = 1-60", Vol.1, Academic Press, New York (1981).

        [1984MuZY]: S.F. Mughabghab, "Neutron Cross Sections, Neutron Resonance 
        Parameters and Thermal Cross Sections, Part B, Z = 61-100", Vol. 2, 
        Academic Press, New York (1984).

        [2018MuZY]: S.F. Mughabghab, "Atlas of Neutron Resonances, Resonance 
        Properties and Thermal Cross Sections Z = 1-60", 6th ed. Vol.1, 
        Elsevier, Amsterdam (2018).

        [2018MuZZ]: S.F. Mughabghab, "Atlas of Neutron Resonances, Resonance 
        Properties and Thermal Cross Sections Z = 61-102", 6th ed. Vol.2, 
        Elsevier, Amsterdam (2018).

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  target: The target ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the target nucleus passed as an integer 
                     argument.

        Returns: 
            A tuple object containing thermal neutron-capture information 
            associated with the target nucleus:

            [0]: Adopted total radiative thermal neutron capture cross section 
                 from the Atlas of the Neutron Resonances (float);
            [1]: Associated neutron-capture cross section uncertainty (float);
            [2]: Cross section units (str);
            [3]: Reference keynumber (str).

        Examples:
            For 12C(n,g) thermal neutron-capture:
            get_total_cross_section(edata,"C12")

            For 28Si(n,g) thermal neutron-capture:
            get_total_cross_section(edata,"Si28")
        """
        self.list = list
        self.args = args
        total_capture_cs = None
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True
            
        for jdict in self.list:
            try:
                if (len(self.args) == 1 and str(args[0]) == jdict["nucleusTargetID"]) or (len(self.args) == 2 and int(args[0]) == jdict["nucleusTargetZ"] and int(args[1]) == jdict["nucleusTargetA"]):
                    for each_n in jdict["neutronCaptureNormalization"]:
                        for each_r in each_n["normalizationRecord"]:
                            total_capture_cs = (each_r["adoptedTotalThermalCaptureCrossSection"],each_r["dAdoptedTotalThermalCaptureCrossSection"],each_r["unitAdoptedCrossSection"],each_r["keyNumber"])
            except ValueError:
                WRONG_INPUTS = True
        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_total_cross_section(edata,\"Si28\")")
            print("or:")
            print(" get_total_cross_section(edata,14,28)")
            return                        
        if total_capture_cs == None:
            print("No target nucleus in EGAF file for defined input.")
            return None
        else:
            return total_capture_cs

    def get_abundance(self,list,*args):
        """Natural isotopic abundances of a measured sample in the EGAF 
        database.

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  target: The target ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the target nucleus passed as an integer 
                     argument.

        Returns: 
            A tuple object containing isotopic abundance associated with the 
            (n,g) target sample:

            [0]: Natural isotopic abundance of measured target sample (float);
            [1]: Natural abundance associated uncertainty (float).

        Examples:
            For 27Al abundance:
            get_abundance(edata,"Al27")

            For 35Cl abundance:
            get_abundance(edata,"Cl35")
        """
        self.list = list
        self.args = args
        abundance = None
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True        
        for jdict in self.list:
            try:
                if (len(self.args) == 1 and str(args[0]) == jdict["nucleusTargetID"]) or (len(self.args) == 2 and int(args[0]) == jdict["nucleusTargetZ"] and int(args[1]) == jdict["nucleusTargetA"]):
                    for each_n in jdict["neutronCaptureNormalization"]:
                        for each_r in each_n["normalizationRecord"]:
                            abundance = (each_r["naturalIsotopicAbundance"],each_r["dNaturalIsotopicAbundance"])
            except ValueError:
                WRONG_INPUTS = True
        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_abundance(edata,\"Al27\")")
            print("or:")
            print(" get_abundance(edata,13,27)")
            return                              
        if abundance == None:
            print("No target nucleus in EGAF file for defined input.")
            return None
        else:
            return abundance
        return abundance

    def get_all_total_cross_sections(self,list):
        """Adopted total radiative thermal neutron-capture cross sections for 
        all measured targets in the EGAF database.  All values are taken from 
        various editions of the "Atlas of Neutron Resonances":
        
        [1981MuZQ]: S.F. Mughabghab, M.Divadeenam, N.E.Holden, "Neutron Cross 
        Sections, Neutron Resonance Parameters and Thermal Cross 
        Sections, Part A, Z = 1-60", Vol.1, Academic Press, New York (1981).

        [1984MuZY]: S.F. Mughabghab, "Neutron Cross Sections, Neutron Resonance 
        Parameters and Thermal Cross Sections, Part B, Z = 61-100", Vol. 2, 
        Academic Press, New York (1984).

        [2018MuZY]: S.F. Mughabghab, "Atlas of Neutron Resonances, Resonance 
        Properties and Thermal Cross Sections Z = 1-60", 6th ed. Vol.1, 
        Elsevier, Amsterdam (2018).

        [2018MuZZ]: S.F. Mughabghab, "Atlas of Neutron Resonances, Resonance 
        Properties and Thermal Cross Sections Z = 61-102", 6th ed. Vol.2, 
        Elsevier, Amsterdam (2018).

        Arguments:
            list: A list of EGAF-data JSON objects.        

        Returns: 
            An alphabetically-sorted dictionary object.  The key is a string 
            object associated with the (n,g) target in the EGAF database.  The 
            corresponding value is a tuple with the following elements

            [0]: Residual compound nucleus associated with the target nucleus 
                 given in the key (str);
            [1]: Adopted total radiative thermal neutron capture cross section 
                 from the Atlas of the Neutron Resonances for the (n,g) target 
                 nucleus(float);
            [2]: Associated neutron-capture cross section uncertainty (float);
            [3]: Cross section units (str);
            [4]: Reference keynumber (str).

        Examples:
            get_all_total_cross_sections(edata)
        """
        self.list = list
        egaf_targets = BaseEGAF.egaf_target_list(self,self.list)
        all_cs_dict = {}
        for target in egaf_targets:
            regex_symbol = re.compile(r'\D+')
            regex_target_mass = re.compile(r'\d+')
            symbol = regex_symbol.findall(target)[0]
            target_mass = regex_target_mass.findall(target)[0]
            residual_mass = str(int(target_mass)+1)
            residual = symbol+residual_mass
            sigma_0 = CrossSection.get_total_cross_section(self,self.list,target)
            all_cs_dict.update({target: (residual, sigma_0[0], sigma_0[1], sigma_0[2], sigma_0[3])})
        return all_cs_dict

    def get_all_abundances(self,list):
        """Natural isotopic abundances of all measured samples in the EGAF 
        database.

        Arguments:
            list: A list of EGAF-data JSON objects.

        Returns: 
            An alphabetically-sorted dictionary object.  The key is a string 
            object associated with the (n,g) target in the EGAF database.  The 
            corresponding value is a tuple with the following elements

            [0]: Natural isotopic abundance of measured target sample (float);
            [1]: Natural abundance associated uncertainty (float).

        Examples:
            get_all_abundances(edata)
        """
        self.list = list
        egaf_targets = BaseEGAF.egaf_target_list(self,self.list)
        all_abundances = {}
        for target in egaf_targets:
            abundance = CrossSection.get_abundance(self,self.list,target)
            all_abundances.update({target: (abundance[0], abundance[1])})
        return all_abundances
