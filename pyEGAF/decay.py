from .base_egaf import *
from .separation import Separation
from .cross_section import CrossSection


class Levels(CrossSection):
    __doc__="""Class to handle levels in EGAF data sets."""

    def __init__(self):
        BaseEGAF.__init__(self)
        Meta.__init__(self)
        Uncertainties.__init__(self)
        Separation.__init__(self)
        CrossSection.__init__(self)

    def get_residual_levels(self,list,*args):
        """Levels and associated properties of the residual compound nucleus 
        populated in radiative thermal neutron capture.

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns: 
            A list object containing decay-scheme information associated with 
            the levels of the residual compound-nucleus populated following 
            radiative thermal neutron capture.

            [0]: Level index corresponding to populated level (int);
            [1]: Associated level energy populated (float);
            [2]: Associated level energy uncertainty (float);           
            [3]: Number of gammas deexciting level (int);
            [4]: Number of spin-parity assignments for level (int);
            [5]: Spin assignment for level (float).  If more than one spin is 
                 permissable, the first permutation is given.
            [6]: Spin flag (int).  Only permitted integers are: 
                 1: Firm spin assignment; 
                 -1: Tentative spin assignment.
            [7]: Parity assignment for level (int).  If more than one parity is 
                 permissable, the first permutation is given.  Parity 
                 assignments take the following integers only:
                 1: Positive parity;
                 -1: Negative parity;
                 0: No parity given.
            [8]: Parity flag (int).  Only permitted integers are: 
                 1: Firm parity assignment; 
                 -1: Tentative parity assignment.

        Examples:
            get_residual_levels(edata, "Si29")
            get_residual_levels(edata, 14, 29)
        """
        self.list = list
        self.args = args
        level_scheme = []
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True
            
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):        
                    for each_l in jdict["levelScheme"]:
                        level_index = each_l["levelIndex"]
                        level_energy = each_l["levelEnergy"]
                        d_level_energy = each_l["dLevelEnergy"]
                        num_gammas = each_l["numberOfGammas"]
                        num_spins = each_l["numberOfSpins"]
                        spin, parity = None, None
                        j_firm, pi_firm = None, None
                        for each_s in each_l["spins"]:
                            if num_spins > 1:
                                if each_s["spinIndex"] == 0:
                                    spin = each_s["spinReal"]
                                    parity = each_s["parity"]
                                    if each_s["spinIsTentative"] == False:
                                        #j_firm = str("J_FIRM")
                                        j_firm = 1
                                    else:
                                        #j_firm = str("J_TENTATIVE")
                                        j_firm = -1
                                    if each_s["parityIsTentative"] == False:
                                        #pi_firm = str("PI_FIRM")
                                        pi_firm = 1
                                    else:
                                        #pi_firm = str("PI_TENTATIVE")
                                        pi_firm = -1
                            else:
                                spin = each_s["spinReal"]
                                parity = each_s["parity"]
                                if each_s["spinIsTentative"] == False:
                                    #j_firm = str("J_FIRM")
                                    j_firm = 1
                                else:
                                    #j_firm = str("J_TENTATIVE")
                                    j_firm = -1
                                if each_s["parityIsTentative"] == False:
                                    #pi_firm = str("PI_FIRM")
                                    pi_firm = 1
                                else:
                                    #pi_firm = str("PI_TENTATIVE")
                                    pi_firm = -1

                        level_scheme.append([level_index, level_energy, d_level_energy, num_gammas, num_spins, spin, j_firm, parity, pi_firm])

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_residual_levels(edata,\"Na24\")")
            print("or:")
            print(" get_residual_levels(edata,11,24)")
            return
        
        if level_scheme == []:
            print("No decay-scheme data available for defined compound nucleus.")
            return
        else:
            return level_scheme
        

    def find_multiple_jpi(self,list,*args):
        """Finds all levels in decay-scheme of residual nucleus with multiple 
        spin-parity permutations.     

        Arguments:                                       
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.
        
        Returns:
            A list object containing all spin-parity permutations for each 
            level.  The list elements for each level correspond to:

            [0]: Level index (int);
            [1]: Associated level energy (float);
            [2]: Associated level energy uncertainty (float);
            [3]: Number of spin-parity assignments for level (int);
            [4]: Spin-parity sequence (int):
                 0 for the first,
                 1 for the second, etc.
            [5]: Spin assignment for level (float);
            [6]: Parity assignment for level (int):
                 1: Positive parity;
                 -1: Negative parity;
                 0: No parity given.

        Examples:
            find_multiple_jpi(edata,"Na24")
            find_multiple_jpi(edata,11,24)
        """
        self.list = list
        self.args = args
        level_scheme = []
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True
            
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):        
                    for each_l in jdict["levelScheme"]:
                        level_index = each_l["levelIndex"]
                        level_energy = each_l["levelEnergy"]
                        d_level_energy = each_l["dLevelEnergy"]
                        num_spins = each_l["numberOfSpins"]
                        if num_spins > 1:
                            for each_s in each_l["spins"]:
                                spin_index = each_s["spinIndex"]
                                spin = each_s["spinReal"]
                                parity = each_s["parity"]

                                level_scheme.append([level_index, level_energy, d_level_energy, num_spins, spin_index, spin, parity])

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" find_multiple_jpi(edata,\"Na24\")")
            print("or:")
            print(" find_multiple_jpi(edata,11,24)")
            return
                            
        if level_scheme == []:
            print("No levels with multiple spin-parity assignments.")
            return
        else:
            return level_scheme

        
    def find_unique_jpi(self,list,*args):
        """Finds all levels in decay-scheme of residual nucleus with unique 
        spin-parity permutations, where both spin and parity are firmly 
        assigned.

        Arguments:                                       
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.
        
        Returns:
            A list object containing all levels with unique spin-parity 
            assignments.  The list elements for each level correspond to:

            [0]: Level index (int);
            [1]: Associated level energy (float);
            [2]: Associated level energy uncertainty (float);
            [3]: Spin-parity sequence (int):
                 0 for the first,
                 1 for the second, etc.
            [4]: Spin assignment for level (float);
            [5]: Parity assignment for level (int):
                 1: Positive parity;
                 -1: Negative parity;
                 0: No parity given.

        Examples:
            find_unique_jpi(edata,"Na24")
            find_unique_jpi(edata,11,24)
        """
        self.list = list
        self.args = args
        level_scheme = []
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True
            
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):        
                    for each_l in jdict["levelScheme"]:
                        level_index = each_l["levelIndex"]
                        level_energy = each_l["levelEnergy"]
                        d_level_energy = each_l["dLevelEnergy"]
                        num_spins = each_l["numberOfSpins"]
                        if num_spins == 1:
                            for each_s in each_l["spins"]:
                                if each_s["spinIsTentative"] == False and each_s["parityIsTentative"] == False:
                                    spin_index = each_s["spinIndex"]
                                    spin = each_s["spinReal"]
                                    parity = each_s["parity"]

                                    level_scheme.append([level_index, level_energy, d_level_energy, spin_index, spin, parity])

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" find_unique_jpi(edata,\"Na24\")")
            print("or:")
            print(" find_unique_jpi(edata,11,24)")
            return
                            
        if level_scheme == []:
            print("No levels with unique spin-parity assignments.")
            return
        else:
            return level_scheme                                
    
    def find_isomers(self,list,*args,**kwargs):
        """Finds all isomeric levels in decay-scheme of the residual compound 
        nucleus.

        Notes:
            Median-symmetrized values are adopted for the halflife whereupon
            asymmetric quantities are encountered in the source EGAF data set.
            
        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

            kwargs: An additional keyword argument is required for the halflife
                    time units:

                    units='best'     
                    units='seconds'  
                    units='s'       

                    Only the above keyword arguments (case insensitive) are 
                    acceptable.  The 'best' units are those parsed from the 
                    original EGAF file. 

        Returns: 
            A list object with all levels that associated halflife information:

            [0]: Level index (int);
            [1]: Associated level energy (float);
            [2]: Associated level energy uncertainty (float);
            [3]: Halflife in 'best' units or in 'seconds' (float);
            [4]: Halflife uncertainty (float);
            [5]: Halflife time units (str).

        Example:
            find_isomers(edata, "Na24", units="best")
            find_isomers(edata, 11, 24, units="seconds")
        """
        self.list = list
        self.args = args
        UNSPECIFIED_UNIT = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True
        
        level_scheme = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):        
                    for each_l in jdict["levelScheme"]:
                        level_index = each_l["levelIndex"]
                        level_energy = each_l["levelEnergy"]
                        d_level_energy = each_l["dLevelEnergy"]
                        if each_l["levelIsIsomer"] == True:
                            for units in kwargs.values():
                                if units.lower() == str('best'):
                                    for each_i in each_l["isomerDecay"]:
                                        half_life = each_i["halfLifeBest"]
                                        d_half_life = each_i["dHalfLifeBest"]
                                        unit_half_life = each_i["unitHalfLifeBest"]

                                        level_scheme.append([level_index, level_energy, d_level_energy, half_life, d_half_life, unit_half_life])
                                elif (units.lower() == str('seconds')) or (units.lower() == str('s')):
                                    for each_i in each_l["isomerDecay"]:
                                        half_life = each_i["halfLifeConverted"]
                                        d_half_life = each_i["dHalfLifeConverted"]
                                        unit_half_life = each_i["unitHalfLifeConverted"]

                                        level_scheme.append([level_index, level_energy, d_level_energy, half_life, d_half_life, unit_half_life])
                                else:
                                    UNSPECIFIED_UNIT = True
            except ValueError:
                WRONG_INPUTS = True
                                    
        if kwargs == {} or kwargs == None:
            print("A keyword argument is required for the desired halflife units.")
            print("Please pass one of the following arguments:")
            print("units='best'")
            print("units='seconds'")
            print("units='s'")
        if UNSPECIFIED_UNIT == True:
            print("Inavlid units for halflife.")
            print("Only the following keyword arguments are accepted:")
            print("units='best'")
            print("units='seconds'")
            print("units='s'")

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" find_isomers(edata,\"Na24\",units=<str>)")
            print("or:")
            print(" find_isomers(edata,11,24,units=<str>)")
            return
                            
        if level_scheme == []:
            print("No isomeric levels in decay scheme.")
            return
        else:
            return level_scheme 

        
class Gammas(Levels):
    __doc__="""Class to handle gammas in EGAF data sets."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_gammas(self,list,*args,**kwargs):
        """Gamma-ray energies and intensities together with associated 
        transition properties including level energies and internal-conversion 
        coefficients.  The returned data corresponds to the properties of the 
        residual compound nucleus produced in thermal neutron-capture reactions.

        Notes:
            (i) Total internal-conversion coefficients (where given) are 
            calculated values obtained using the BrIcc code:

            [2008Ki07] - T.Kibedi et al., Nucl. Instrum. Methods Phys. Res. 
            Sect. A 589, 202 (2008).

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

            kwargs: An additional keyword arguments is required for the 
                    gamma-ray intensity units:

                    intensity='elemental'  : Elemental partial gamma-ray cross 
                                             sections.
                    intensity='isotopic'   : Isotopic partial gamma-ray cross 
                                             sections.
                    intensity='population' : Populations per neutron capture.
                    intensity='relative'   : Relative intensity (%) to the 
                                             strongest transition in the same
                                             nucleus.

        Returns:
            A numpy array containing the following elements associated with the
            gamma decay of the residual compound nucleus:

            [0]: Level index corresponding to initial level (float);
            [1]: Level index corresponding to final level (float);
            [2]: Associated initial level energy in keV (float);
            [3]: Associated final level energy in keV (float);
            [4]: Deexcitation gamma-ray energy in keV (float);
            [5]: Deexcitation gamma-ray energy uncertainty (float);
            [6]: Gamma-ray <intensity> according to keyword argument provided 
                 (float);
            [7]: Gamma-ray <intensity> uncertainty (float);
            [8]: BrIcc-calculated total internal-conversion coefficient 
                 (float);
            [9]: Total internal-conversion coefficient uncertainty (float).

        Examples:
            For isotopic partial gamma-ray cross sections:
            get_gammas(edata, "Y90", intensity="isotopic") 
            get_gammas(edata, 39, 90, intensity="isotopic")
 
            For elemental partial gamma-ray cross sections:
            get_gammas(edata, "Si29", intensity="elemental") 
            get_gammas(edata, 14, 29, intensity="elemental") 

            For populations per neutron capture:
            get_gammas(edata, "Na24", intensity="population") 
            get_gammas(edata, 11, 24, intensity="population") 

            For intensities relative to the strongest transition (100%):
            get_gammas(edata, "C13", intensity="relative") 
            get_gammas(edata, 6, 13, intensity="relative") 
        """
        self.list = list
        self.args = args
        DECAY_SCHEME_EXISTS = False
        UNSPECIFIED_INTENSITY = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True

        gamma_spec = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):        

                    DECAY_SCHEME_EXISTS = True
                    residual = jdict["nucleusID"]
                    
                    g = Gammas()
                    
                    for each_l in jdict["levelScheme"]:
                        if each_l["numberOfGammas"] > 0:                            
                            for each_g in each_l["gammaDecay"]:
                                initial_index = each_g["levelIndexInitial"]
                                final_index = each_g["levelIndexFinal"]
                                initial_energy = each_g["levelEnergyInitial"]
                                final_energy = each_g["levelEnergyFinal"]
                                gamma_energy = each_g["gammaEnergy"]
                                d_gamma_energy = each_g["dGammaEnergy"]
                                alpha = each_g["calculatedTotalInternalConversionCoefficient"]
                                d_alpha = each_g["dCalculatedTotalInternalConversionCoefficient"]
                                for intensity in kwargs.values():
                                    if intensity.lower() == str("elemental"):
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["partialElementalCrossSection"]
                                            d_gamma_intensity = each_i["dPartialElementalCrossSection"]

                                            gamma_spec.append([initial_index, final_index, initial_energy, final_energy, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity, alpha, d_alpha])
                                    elif intensity.lower() == str("isotopic"):
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["partialIsotopicCrossSection"]
                                            d_gamma_intensity = each_i["dPartialIsotopicCrossSection"]

                                            gamma_spec.append([initial_index, final_index, initial_energy, final_energy, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity, alpha, d_alpha])
                                    elif intensity.lower() == str("population"):
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["populationPerNeutronCapture"]
                                            d_gamma_intensity = each_i["dPopulationPerNeutronCapture"]

                                            gamma_spec.append([initial_index, final_index, initial_energy, final_energy, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity, alpha, d_alpha])
                                            
                                    elif intensity.lower() == str("relative"):
                                        # Call method recursively assuming elemental intensities
                                        spe = g.get_gammas(self.list, residual, intensity = 'elemental')
                                        df = pd.DataFrame({'I':spe[:,6], 'dI':spe[:,7]})
                                        max_I = df['I'].loc[df['I'].idxmax()]
                                        dI_at_max_I = df['dI'].loc[df['I'].idxmax()]
                                        
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["partialElementalCrossSection"]
                                            d_gamma_intensity = each_i["dPartialElementalCrossSection"]
                                            
                                            rel_gamma_intensity = (gamma_intensity/max_I) * 100
                                            if d_gamma_intensity > 0.0:
                                                d_rel_gamma_intensity = rel_gamma_intensity * (d_gamma_intensity/gamma_intensity)
                                            else:
                                                d_rel_gamma_intensity = 0.0
                                                
                                            gamma_intensity = rel_gamma_intensity
                                            d_gamma_intensity = d_rel_gamma_intensity

                                            gamma_spec.append([int(initial_index), int(final_index), initial_energy, final_energy, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity, alpha, d_alpha])
                                    else:
                                        UNSPECIFIED_INTENSITY = True
            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_gammas(edata, \"Y90\", intensity=<str>)")
            print("or:")
            print(" get_gammas(edata, 39, 90, intensity=<str>)")
            return
        
        if kwargs == {} or kwargs == None:
            print("A keyword argument is required for the desired intensity units.")
            print("Please pass one of the following arguments:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            print("intensity='population'")
            print("intensity='relative'")
            return
        
        if UNSPECIFIED_INTENSITY == True:
            print("Incorrect intensity specified.")
            print("Only the following keyword arguments are accepted:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            print("intensity='population'")
            print("intensity='relative'")
            return

        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return            
        
        if gamma_spec == []:
            print("No gammas in decay scheme")
            return
        else:
            return np.array(gamma_spec)

    
    def find_all_gammas_feeding_gs(self,list,*args,**kwargs):
        """Finds all gamma-ray transitions that feed the ground state of the 
        compound nucleus in direct single-step transitions.  For a ground state
        defined with level index 0, the first, second, third, etc. excited 
        states are labeled 1, 2, and 3, respectively, only those transitions 
        with an associated final-level index of 0 are listed, i.e., 1->0, 2->0,
        3->0, etc.  Primary gamma-ray transitions feeding the ground state 
        directly (i.e., capture state -> 0) are also included where observed.
        Associated decay-scheme information including level indices, energies, 
        and internal-conversion coefficients are also returned.

        Notes:
            (i) Total internal-conversion coefficients (where given) are 
            calculated values obtained using the BrIcc code:

            [2008Ki07] - T.Kibedi et al., Nucl. Instrum. Methods Phys. Res. 
            Sect. A 589, 202 (2008).

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

            kwargs: An additional keyword arguments is required for the 
                    gamma-ray intensity units:

                    intensity='isotopic'   : Isotopic partial gamma-ray cross 
                                             sections.
                    intensity='population' : Populations per neutron capture.

        Returns:
            A list object containing  all gamma-ray transitions that feed the 
            ground state in direct single-step transitions.  The list elements 
            correspond to:

            [0]: Level index corresponding to initial level (int);
            [1]: Level index corresponding to final level (int);
            [2]: Associated initial level energy in keV (float);
            [3]: Associated final level energy in keV (float);
            [4]: Deexcitation gamma-ray energy in keV (float);
            [5]: Deexcitation gamma-ray energy uncertainty (float);
            [6]: Gamma-ray <intensity> according to keyword argument provided 
                 (float);
            [7]: Gamma-ray <intensity> uncertainty (float);
            [8]: BrIcc-calculated total internal-conversion coefficient 
                 (float);
            [9]: Total internal-conversion coefficient uncertainty (float).

        Examples:
            For 12C(n,g)13C gamma rays feeding the ground state (isotopic 
            partial gamma-ray cross sections):
            find_all_gammas_feeding_gs(edata, "C13", intensity="isotopic")
            find_all_gammas_feeding_gs(edata, "C13", intensity="isotopic")

            For 28Si(n,g)29Si gamma rays feeding the ground state (populations 
            per neutron capture):
            find_all_gammas_feeding_gs(edata, "Si29", intensity="population")
            find_all_gammas_feeding_gs(edata, 14, 29, intensity="population")
        """
        self.list = list
        self.args = args
        UNSPECIFIED_INTENSITY = False
        DECAY_SCHEME_EXISTS = False
        DIRECT_FEEDING_GS = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True

        feeding_gs = []            
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    DECAY_SCHEME_EXISTS = True
                    for each_l in jdict["levelScheme"]:
                        if each_l["numberOfGammas"] > 0:
                            for each_g in each_l["gammaDecay"]:
                                if each_g["gammaFeedsGroundState"] == True:
                                    DIRECT_FEEDING_GS = True

                                    initial_index = each_g["levelIndexInitial"]
                                    final_index = each_g["levelIndexFinal"]
                                    initial_energy = each_g["levelEnergyInitial"]
                                    final_energy = each_g["levelEnergyFinal"]
                                    gamma_energy = each_g["gammaEnergy"]
                                    d_gamma_energy = each_g["dGammaEnergy"]
                                    alpha = each_g["calculatedTotalInternalConversionCoefficient"]
                                    d_alpha = each_g["dCalculatedTotalInternalConversionCoefficient"]
                                    cs, d_cs = None, None
                                    for each_i in each_g["gammaAbsoluteIntensities"]:
                                        for intensity in kwargs.values():
                                            if intensity.lower() == str("isotopic"):
                                                cs = each_i["partialIsotopicCrossSection"]
                                                d_cs = each_i["dPartialIsotopicCrossSection"]

                                            elif intensity.lower() == str("population"):
                                                cs = each_i["populationPerNeutronCapture"]
                                                d_cs = each_i["dPopulationPerNeutronCapture"]
                                            else:
                                                UNSPECIFIED_INTENSITY = True

                                    feeding_gs.append([initial_index, final_index, initial_energy, final_energy, gamma_energy, d_gamma_energy, cs, d_cs, alpha, d_alpha])
            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" find_all_gammas_feeding_gs(edata, \"C13\", intensity=\"isotopic\")")
            print("or:")
            print(" find_all_gammas_feeding_gs(edata, 6, 13, intensity=\"isotopic\")")
            return

        if kwargs == {} or kwargs == None:
            print("A keyword argument is required for the desired intensity units.")
            print("Please pass one of the following arguments:")
            print("intensity='isotopic'")
            print("intensity='population'")
            return
        
        if UNSPECIFIED_INTENSITY == True:
            print("Incorrect intensity specified.")
            print("Only the following keyword arguments are accepted:")
            print("intensity='isotopic'")
            print("intensity='population'")
            return
                                
        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return
        
        if DIRECT_FEEDING_GS == False:
            print("No transitions feeding ground state directly.")
            return

        if feeding_gs == []:
            print("No transitions feeding ground state directly.")
            return
        else:
            return feeding_gs


    def get_gamma_types(self,list,*args,**kwargs):
        """Lists all gamma-ray transitions of a particular type: Either 
        'primary' or 'secondary' gamma decay from the residual compound nucleus.
        The primary gammas are those that deexcite the capture state, and the 
        secondary gammas are those associated associated with transitions 
        between low-lying levels.  Associated decay-scheme information including
        level indices, energies, and internal-conversion coefficients are also 
        returned.

        Notes:
            (i) Total internal-conversion coefficients (where given) are 
            calculated values obtained using the BrIcc code:

            [2008Ki07] - T.Kibedi et al., Nucl. Instrum. Methods Phys. Res. 
            Sect. A 589, 202 (2008).

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

            kwargs: Two additional keyword arguments are required: 

                    (i) Gamma-ray intensity units:

                    intensity='elemental'  : Elemental partial gamma-ray cross 
                                             sections.
                    intensity='isotopic'   : Isotopic partial gamma-ray cross 
                                             sections.
                    intensity='population' : Populations per neutron capture.
                    intensity='relative'   : Relative intensity (%) to the 
                                             strongest transition in the same
                                             nucleus.

                    (ii) Gamma-ray type:
        
                    gammas='primary'   : Gammas that originate at the 
                                         neutron-capture state.
                    gammas='secondary' : Gammas associated with transitions 
                                         betwenn low-lying levels.

        Returns:
            A list object containing  all gamma-ray transitions of a particular 
            type: Either 'primary' or 'secondary'.  The list elements 
            correspond to:

            [0]: Level index corresponding to initial level (int);
            [1]: Level index corresponding to final level (int);
            [2]: Associated initial level energy in keV (float);
            [3]: Associated final level energy in keV (float);
            [4]: Deexcitation gamma-ray energy in keV (float);
            [5]: Deexcitation gamma-ray energy uncertainty (float);
            [6]: Gamma-ray <intensity> according to keyword argument provided 
                 (float);
            [7]: Gamma-ray <intensity> uncertainty (float);
            [8]: BrIcc-calculated total internal-conversion coefficient 
                 (float);
            [9]: Total internal-conversion coefficient uncertainty (float);
            [10]: Gamma-ray transition type: "primary" or "secondary" (str).

        Examples:
            For 12C(n,g)13C gamma rays (elemental partial gamma-ray cross 
            sections):
            get_gamma_types(edata, "C13", gammas='secondary', intensity='elemental')
            get_gamma_types(edata, 6, 13, gammas='secondary', intensity='elemental')
            get_gamma_types(edata, "C13", gammas='primary', intensity='elemental')
            get_gamma_types(edata, 6, 13, gammas='primary', intensity='elemental')
        """
        self.list = list
        self.args = args
        DECAY_SCHEME_EXISTS = False
        UNSPECIFIED_INTENSITY = None
        PRIMARIES = False
        SECONDARIES = False
        NO_GAMMAS = False
        gamma_type = None
        intensity_type = None
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True        
        
        gammas = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    DECAY_SCHEME_EXISTS = True
                    for kwarg in kwargs.values():
                        if (kwarg.lower() == "primary") or (kwarg.lower() == "secondary"):
                            gamma_type = str(kwarg.lower())
                        if (kwarg.lower() == "elemental") or (kwarg.lower() == "isotopic") or (kwarg.lower() == "population") or (kwarg.lower() == "relative"):
                            intensity_type = str(kwarg.lower())

                    for each_l in jdict["levelScheme"]:
                        if each_l["numberOfGammas"] > 0:
                            for each_g in each_l["gammaDecay"]:
                                if each_g["gammaTransitionType"] == gamma_type:

                                    if gamma_type.lower() == 'primary':
                                        PRIMARIES = True
                                    elif gamma_type.lower() == 'secondary':
                                        SECONDARIES = True
                                    else:
                                        NO_GAMMAS = True

                                    initial_index = each_g["levelIndexInitial"]
                                    final_index = each_g["levelIndexFinal"]
                                    initial_energy = each_g["levelEnergyInitial"]
                                    final_energy = each_g["levelEnergyFinal"]
                                    gamma_energy = each_g["gammaEnergy"]
                                    d_gamma_energy = each_g["dGammaEnergy"]
                                    alpha = each_g["calculatedTotalInternalConversionCoefficient"]
                                    d_alpha = each_g["dCalculatedTotalInternalConversionCoefficient"]
                                    cs, d_cs = None, None

                                    # Find maximum intensity 
                                    g = Gammas()
                                    residual = jdict["nucleusID"]
                                    spe = g.get_gammas(self.list, residual, intensity = 'elemental')
                                    df = pd.DataFrame({'I':spe[:,6], 'dI':spe[:,7]})
                                    max_I = df['I'].loc[df['I'].idxmax()]
                                    dI_at_max_I = df['dI'].loc[df['I'].idxmax()]
                                    
                                    #for intensity in kwargs.values():
                                    #    if intensity == "isotopic":
                                    if intensity_type != None:
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            if intensity_type.lower() == 'elemental':
                                                cs = each_i["partialElementalCrossSection"]
                                                d_cs = each_i["dPartialElementalCrossSection"]
                                            elif intensity_type.lower() == 'isotopic':
                                                cs = each_i["partialIsotopicCrossSection"]
                                                d_cs = each_i["dPartialIsotopicCrossSection"]
                                            elif intensity_type.lower() == 'population':
                                                cs = each_i["populationPerNeutronCapture"]
                                                d_cs = each_i["dPopulationPerNeutronCapture"]
                                            elif intensity_type.lower() == 'relative':
                                                gamma_intensity = each_i["partialElementalCrossSection"]
                                                d_gamma_intensity = each_i["dPartialElementalCrossSection"]
                                            
                                                rel_gamma_intensity = (gamma_intensity/max_I) * 100
                                                if d_gamma_intensity > 0.0:
                                                    d_rel_gamma_intensity = rel_gamma_intensity * (d_gamma_intensity/gamma_intensity)
                                                else:
                                                    d_rel_gamma_intensity = 0.0
                                                
                                                cs = rel_gamma_intensity
                                                d_cs = d_rel_gamma_intensity
                                                
                                            else:
                                                UNSPECIFIED_INTENSITY = True
                                    else:
                                        UNSPECIFIED_INTENSITY = True

                                    gammas.append([initial_index, final_index, initial_energy, final_energy, gamma_energy, d_gamma_energy, cs, d_cs, alpha, d_alpha, gamma_type])
            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_gamma_types(edata, \"C13\", intensity=<str>, gammas=<str>)")
            print("or:")
            print(" get_gamma_types(edata, 6, 13, intensity=<str>, gammas=<str>)")
            return
                                
        if kwargs == {} or kwargs == None:
            print("Two keyword arguments are required.\n")
            print("Please pass one of the following arguments for the required intensity units:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            print("intensity='population'")
            print("intensity='relative'")
            print("\nPlease pass one of the following arguments for the required gamma-ray types:")
            print("gammas='primary'")
            print("gammas='secondary'")
            return

        if len(kwargs)==1 or len(kwargs)>2:
            print("Two keyword arguments are required.\n")
            print("Please pass one of the following arguments for the required intensity units:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            print("intensity='population'")
            print("intensity='relative'")
            print("\nPlease pass one of the following arguments for the required gamma-ray types:")
            print("gammas='primary'")
            print("gammas='secondary'")
            return
        
        if UNSPECIFIED_INTENSITY == True:
            print("Incorrect intensity specified.")
            print("Only the following keyword arguments are accepted:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            print("intensity='population'")
            print("intensity='relative'")
            return
        
        if gamma_type == "primary" and NO_GAMMAS == True:
            print("No primary gamma rays in decay scheme.")
        elif gamma_type == "secondary" and NO_GAMMAS == True:
            print("No secondary gamma rays in decay scheme.")
        elif gamma_type == None:
            print("Incorrect gamma-type specified.")
            print("Only the following keyword arguments are accepted:")
            print("gammas='primary'") 
            print("gammas='secondary'")
            
        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments.")
            return

        if gammas == []:
            print("No gammas of the required type found in decay scheme.")
            return
        else:
            return gammas
    
    def find_gamma(self, list, float, tolerance=0.5, **kwargs):
        """Searches for all gamma rays at a specified energy.  By default the 
        search will find all gamma rays within +/- 0.5 keV of the specified 
        energy passed to the function.

        Arguments:
            list: A list of EGAF-data JSON objects.
            float: Gamma-ray energy in keV (float).
            tolerance: Limit of the energy search; by default +/- 0.5 keV.
            kwargs: The function takes one of the following keyword arguments 
                    according to the desired intensity output:

                    intensity='elemental' : Elemental partial gamma-ray cross 
                                            section.
                    intensity='isotopic'  : Isotopic partial gamma-ray cross 
                                            section.
                    intensity='population': Population per neutron capture.
                    intensity='relative'  : Relative intensity (%) to the 
                                            strongest transition in the same
                                            nucleus.

                    If no keyword argument is provided "relative" intensities 
                    will be adopted by default.

        Returns: 
            A DataFrame object containing the target isotope, residual (n,g) 
            compound-nucleus isotope, gamma-ray energy and its uncertainty, 
            gamma-ray intensity and its uncertainty.
        
        Examples:
            (i) Find all isotopes containing gamma rays within 1273+/-0.5 keV:
            
            For isotopic cross sections:
            find_gamma(edata, 1273) 
            find_gamma(edata, 1273, intensity='isotopic')

            For elemental cross sections:
            find_gamma(edata, 1273, intensity='elemental')

            For populations per neutron capture:
            find_gamma(edata, 1273, intensity='population')

            For relative intensities:
            find_gamma(edata, 1273, intensity='relative')

            (ii) Finds all isotopes containing gamma rays within 1273+/-2.5 keV:

            find_gamma(edata, 1273, 2.5)
            find_gamma(edata, 1273, 2.5, intensity='isotopic')
            find_gamma(edata, 1273, 2.5, intensity='elemental')
            find_gamma(edata, 1273, 2.5, intensity='population')
            find_gamma(edata, 1273, 2.5, intensity='relative')
        """
        self.list = list
        self.float = float

        if kwargs == {} or kwargs == None:
            # Assign default keyword argument for gamma-ray intensities
            kwargs = {'intensity': 'relative'}
            print("No intensity keyword argument provided.")
            print("Default \"relative\" intensities will be adopted.")

        NO_INTENSITY_OPTION = False
        isotope_list = []
        for jdict in self.list:
            if jdict["datasetType"] == "evaluatedGammarayActivationFile":
                for each_l in jdict["levelScheme"]:
                    if each_l["numberOfGammas"] > 0:
                        for each_g in each_l["gammaDecay"]:
                            gamma_energy = each_g["gammaEnergy"]

                            if gamma_energy >= (self.float - tolerance) and gamma_energy <= (self.float + tolerance):

                                target = jdict["nucleusTargetID"]
                                residual = jdict["nucleusID"]

                                residual_Z = jdict["nucleusZ"]
                                residual_A = jdict["nucleusA"]
                                
                                d_gamma_energy = each_g["dGammaEnergy"]
                                gamma_intensity = None
                                d_gamma_intensity = None

                                # Find maximum intensity
                                g = Gammas()
                                spe = g.get_gammas(self.list, residual, intensity = 'elemental')
                                df = pd.DataFrame({'I':spe[:,6], 'dI':spe[:,7]})
                                max_I = df['I'].loc[df['I'].idxmax()]
                                dI_at_max_I = df['dI'].loc[df['I'].idxmax()]
                                
                                for intensity in kwargs.values():
                                    if intensity == str("elemental"):
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["partialElementalCrossSection"]
                                            d_gamma_intensity = each_i["dPartialElementalCrossSection"]

                                            isotope_list.append([target, residual, residual_Z, residual_A, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity])
                                            
                                    elif intensity == str("isotopic"):
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["partialIsotopicCrossSection"]
                                            d_gamma_intensity = each_i["dPartialIsotopicCrossSection"]

                                            isotope_list.append([target, residual, residual_Z, residual_A, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity])
                                            
                                    elif intensity == str("population"):
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["populationPerNeutronCapture"]
                                            d_gamma_intensity = each_i["dPopulationPerNeutronCapture"]

                                            isotope_list.append([target, residual, residual_Z, residual_A, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity])
                                            
                                    elif intensity == str("relative"):
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["partialElementalCrossSection"]
                                            d_gamma_intensity = each_i["dPartialElementalCrossSection"]
                                            
                                            rel_gamma_intensity = (gamma_intensity/max_I) * 100
                                            if d_gamma_intensity > 0.0:
                                                d_rel_gamma_intensity = rel_gamma_intensity * (d_gamma_intensity/gamma_intensity)
                                            else:
                                                d_rel_gamma_intensity = 0.0
                                                
                                            gamma_intensity = rel_gamma_intensity
                                            d_gamma_intensity = d_rel_gamma_intensity
                                            
                                            isotope_list.append([target, residual, residual_Z, residual_A, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity])
                                            
                                    else:
                                        NO_INTENSITY_OPTION = True
                                        # Defaults to "isotopic" intensities
                                        for each_i in each_g["gammaAbsoluteIntensities"]:
                                            gamma_intensity = each_i["partialElementalCrossSection"]
                                            d_gamma_intensity = each_i["dPartialElementalCrossSection"]
                                            
                                            rel_gamma_intensity = (gamma_intensity/max_I) * 100
                                            if d_gamma_intensity > 0.0:
                                                d_rel_gamma_intensity = rel_gamma_intensity * (d_gamma_intensity/gamma_intensity)
                                            else:
                                                d_rel_gamma_intensity = 0.0
                                                
                                            gamma_intensity = rel_gamma_intensity
                                            d_gamma_intensity = d_rel_gamma_intensity
                                
                                            isotope_list.append([target, residual, residual_Z, residual_A, gamma_energy, d_gamma_energy, gamma_intensity, d_gamma_intensity])

        if NO_INTENSITY_OPTION == True:
            print("I didn't understand the intensity keyword-argument provided.")
            print("Default \"relative\" intensities will be adopted.")

        try:
            isotope_list_s = sorted(isotope_list, key=lambda x: (x[2], x[3]))
            isotope_array = np.array(isotope_list_s)
            target = isotope_array[:,0]
            residual = isotope_array[:,1]
            gamma_energy = isotope_array[:,4]
            d_gamma_energy = isotope_array[:,5]
            gamma_intensity = isotope_array[:,6]
            d_gamma_intensity = isotope_array[:,7]

            pd.set_option('display.max_columns', None)
            pd.set_option('display.max_row', None)


            isotope_df = pd.DataFrame({'Target (n,g)': target, 'Residual (CN)': residual, 'Energy (keV)': gamma_energy, 'dE (keV)': d_gamma_energy, 'Intensity': gamma_intensity, 'dI': d_gamma_intensity})

            return isotope_df

        except IndexError:

            if (isotope_list == []) or len(isotope_list) == 0:
                print("No gammas in EGAF database match specified search criterion: {0} \xb1 {1} keV.\nTry a different energy or expand the search window.".format(self.float,tolerance))
                return

            else:
                raise
    
    def get_strongest_gammas(self, list, *args, **kwargs):
        """Finds up to the three strongest gamma-ray transitions in the 
        residual compound nucleus.
        
        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

            kwargs: An additional keyword arguments is required for the 
                    gamma-ray intensity units:

                    intensity='elemental'  : Elemental partial gamma-ray cross 
                                             sections.
                    intensity='isotopic'   : Isotopic partial gamma-ray cross 
                                             sections.
                    intensity='population' : Populations per neutron capture.
                    intensity='relative'  : Relative intensity (%) to the 
                                            strongest transition in the same
                                            nucleus.

                    If no keyword argument is provided "relative" intensities 
                    will be adopted by default.

        Returns:
            A DataFrame object containing the gamma-ray energy and its 
            associated uncertainty, together with the corresponding gamma-ray 
            intensity and its associated uncertainty.

        Examples:
            For relative intensities:
            get_strongest_gammas(edata, "Si29")
            get_strongest_gammas(edata, "Si29", intensity='relative')

            For isotopic cross sections:
            get_strongest_gammas(edata, "Si29", intensity='isotopic')

            For elemental cross sections:
            get_strongest_gammas(edata, "Si29", intensity='elemental')

            For populations per neutron capture:
            get_strongest_gammas(edata, "Si29", intensity='population')
        """
        self.list = list
        self.args = args

        if kwargs == {} or kwargs == None:
            # Assign default keyword argument for gamma-ray intensities
            kwargs = {'intensity': 'relative'}
            print("No intensity keyword argument provided.")
            print("Default \"relative\" intensities will be adopted.")
        
        DECAY_SCHEME_EXISTS = False
        UNSPECIFIED_INTENSITY = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True

        gamma_list = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):

                    DECAY_SCHEME_EXISTS = True
                    residual = jdict["nucleusID"]

                    g = Gammas()

                    try:
                        intensity = [i for i in kwargs.values()][0]
                        spe = g.get_gammas(self.list, residual, intensity=intensity.lower())
                        spe_sorted = sorted(spe, reverse=True, key=lambda x: x[6])
                        gamma_list = spe_sorted[0:3]
                        #gamma_list = spe_sorted.tolist()
                        
                    except IndexError:
                        # Defaults to "relative" intensities
                        UNSPECIFIED_INTENSITY = True
                        spe = g.get_gammas(self.list, residual, intensity="relative")
                        spe_sorted = sorted(spe, reverse=True, key=lambda x: x[6])
                        gamma_list = spe_sorted[0:3]

            except ValueError:
                WRONG_INPUTS = True
                
        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_gammas(edata, \"Y90\", intensity=<str>)")
            print("or:")
            print(" get_gammas(edata, 39, 90, intensity=<str>)")
            return
        
        if UNSPECIFIED_INTENSITY == True:
            print("I did not understand intensity-keyword argument.\n")
            print("Default \"relative\" intensities will be adopted.")
            return

        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return

        if gamma_list == []:
            print("No gammas in decay scheme")
            return
        else:
            gamma_array = np.array(gamma_list)
            df = pd.DataFrame({'E (keV)':gamma_array[:,4], 'dE (keV)':gamma_array[:,5], 'I':gamma_array[:,6], 'dI':gamma_array[:,7]})
            return df
        
