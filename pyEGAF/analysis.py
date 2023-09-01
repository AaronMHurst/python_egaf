from .base_egaf import *
from .separation import Separation
from .cross_section import CrossSection
from .decay import Levels, Gammas

class Analysis(Gammas):
    __doc__="""Class to perform analysis of EGAF observables."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def intensity_conversion(self,float1,float2,float3,float4):
        """Calculates the total transition intensity corrected for internal 
        conversion.

        Arguments:
            float1 = Gamma-ray intensity (float);
            float2 = Associated uncertainty gamma-ray intensity (float);
            float3 = Internal-conversion coefficient (float);
            float4 = Associated uncertainty internal-conversion 
                     coefficient (float).

        Returns:
            A tuple object with the following elements:

            [0]: Total-transition intensity corrected for conversion (float);
            [1]: Associated uncertainty for total-transition intensity (float).

        Example:
            The 2027.98-keV transition in 29Si has a gamma-ray intensity of 
            0.00075(12) b with a corresponding internal-conversion coefficient 
            of 0.000335(5).  The total transition intensity is found from:

            intensity_conversion(0.00075, 0.00012, 0.000335, 0.000005)
        """
        try:
            self.float1 = float(float1)
            self.float2 = float(float2)
            self.float3 = float(float3)
            self.float4 = float(float4)

            total_intensity = self.float1 * (1.0 + self.float3)
            d_total_intensity = Uncertainties.quad_error(self, total_intensity, self.float1, self.float2, (1.0 + self.float3), self.float4)

            return (total_intensity, d_total_intensity)

        except TypeError:
            # Missing number of positional arguments
            print("Missing positional arguments!")
            print("Four arguments are required:")
            print("intensity_conversion(<float>, <float>, <float>, <float>)")
            return

    
    def modeled_sigma0(self, float1, float2, float3, float4):
        """Determination of the total radiative thermal neutron-capture cross 
        section using the sum of experimental partial gamma-ray cross sections, 
        correct for internal conversion, associated with direct-to-ground-state 
        transitions up to the critical energy, and an assumed modeled population
        per neutron capture contribution feeding the ground state from the 
        quasicontinuum.

        Arguments:
            float1 = Sum of experimental cross-sections of direct ground-state 
                     feeding transitions up to the critical energy (float);
            float2 = Associated uncertainty of the gamma-ray cross section 
                     summation (float);
            float3 = Calculated population per neutron capture feeding to the 
                     ground state from the quasicontinuum (float);
            float4 = Associated uncertainty for the popultion per neutron 
                     capture (float).

        Returns:
            A tuple object with the following elements:

            [0]: Total radiative thermal neutron-capture cross section (float);
            [1]: Associated uncertainty for total-capture cross section (float).

        Example:
            In 28Si(n,g) the sum of experimental partial gamma-ray cross 
            sections corrected for conversion is 0.1870(33) b up to a critical 
            energy of 7057.94 keV.  The corresponding continuum contribution 
            feeding the ground state P0=0.02314(80) assuming the Generalized 
            Lorentzian PSF model together with the Back-Shifted Fermi Gas 
            nuclear LD model.  The corresponding value of the total thermal
            neutron capture cross section is determined as:

            modeled_sigma0(0.1827, 0.0027, 0.02314, 0.00080)
        """
        self.float1 = float(float1)
        self.float2 = float(float2)
        self.float3 = float(float3)
        self.float4 = float(float4)

        sigma_0 = self.float1 / (1.0 - self.float3)
        d_sigma_0 = Uncertainties.quad_error(self, sigma_0, self.float1, self.float2, (1.0 - self.float3), self.float4)

        return (sigma_0, d_sigma_0)

    def modeled_sigma0_ecrit(self,list,int,float1,float2,*args):
        """Determination of the total radiative thermal neutron-capture cross 
        section using the sum of experimental partial gamma-ray cross sections, 
        correct for internal conversion, associated with direct-to-ground-state 
        transitions up to the critical energy defined by the corresponding level
        index, and an assumed modeled population per neutron capture 
        contribution feeding the ground state from the quasicontinuum.

        Arguments:
            list: A list of EGAF-data JSON objects.
            int : Integer index corresponding to the level at Ecrit.
            float1 = Calculated population per neutron capture feeding to the 
                     ground state from the quasicontinuum (float).
            float2 = Associated uncertainty for the popultion per neutron 
                     capture (float).
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns:
            A list object containing the following elements:

            [0]: Critical energy corresponding to the input level index (float);
            [1]: Sum of conversion-corrected partial gamma-ray cross sections
                 corresponding to transitions that feed the ground state 
                 directly (float);
            [2]: Associated uncertainty for cross section summation (float).
            [3]: Total radiative thermal neutron capture cross section 
                 determined using experimental cross sections and modeled 
                 population per neutron capture from the quasicontinuum 
                 feeding the ground state (float).
            [4]: Associated uncertainty for total cross section (float).
             

        Examples:
            Total radiative thermal neutron-capture using experimental data 
            with Ecrit set to level index = 5 and modeled population per 
            neutron capture P0 = 0.29201 +/- 0.17387 from quasicontinuum to 
            ground state in 28Si(n,g)29Si:
            
            modeled_sigma0_ecrit(edata, 5, 0.29201, 0.17387, "Si29")
            modeled_sigma0_ecrit(edata, 5, 0.29201, 0.17387, 14, 29)

            For Ecrit level index = 11 and P0 = 0.02256 +/- 0.00064:

            modeled_sigma0_ecrit(edata, 11, 0.02256, 0.00064, "Si29")
            modeled_sigma0_ecrit(edata, 11, 0.02256, 0.00064, 14, 29)
        """
        self.list = list
        self.int = int
        self.float1 = float(float1)
        self.float2 = float(float2)
        self.args = args

        DECAY_SCHEME_EXISTS = False
        DIRECT_FEEDING_GS = False
        LEVEL_ABOVE_ECRIT = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True

        modeled_expt_cs = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and args[0] == jdict["nucleusZ"] and args[1] == jdict["nucleusA"]):
                    DECAY_SCHEME_EXISTS = True

                    compound_nucleus = jdict["nucleusID"]

                    expt_feeding_gs = []
                    d_expt_feeding_gs = []
                    
                    g = Gammas()
                    gammas = g.find_all_gammas_feeding_gs(self.list, compound_nucleus, intensity="isotopic")
                    capture_state_level = max([g[0] for g in gammas])
                    Ecrit = None
                    if self.int >= capture_state_level:
                        LEVEL_ABOVE_ECRIT = True
                    else:
                        for eachl in jdict["levelScheme"]:
                            if eachl["levelIndex"] == self.int:
                                try:
                                    Ecrit = float(eachl["levelEnergy"])
                                except TypeError:
                                    Ecrit = eachl["levelEnergy"]
                                    
                    if gammas != None and len(gammas) > 0:
                        for gdata in gammas:
                            if gdata[0] <= self.int and gdata[1] == 0:
                                DIRECT_FEEDING_GS = True
                                a = Analysis()
                                converted_cs = a.intensity_conversion(gdata[6],gdata[7],gdata[8],gdata[9])[0]
                                d_converted_cs = a.intensity_conversion(gdata[6],gdata[7],gdata[8],gdata[9])[1]
                                expt_feeding_gs.append(converted_cs)
                                d_expt_feeding_gs.append(d_converted_cs**2)

                        # Check for primary gamma to ground state:
                        # Removing the primary otherwise double counting if using experimental populations per neutron capture.
                        # Primary contribution removed by not appending the list.
                        primaries = g.get_gamma_types(self.list, compound_nucleus, intensity="isotopic", gammas="primary")
                        if len(primaries) > 0:
                            if DIRECT_FEEDING_GS == True:
                                #capture_state_level = max([g[0] for g in gammas])
                                for pgamma in primaries:
                                    if pgamma[10] == "primary":
                                        if pgamma[0] == capture_state_level and pgamma[1] == 0:
                                            a = Analysis()
                                            converted_primary = a.intensity_conversion(pgamma[6],pgamma[7],pgamma[8],pgamma[9])[0]
                                            d_converted_primary = a.intensity_conversion(pgamma[6],pgamma[7],pgamma[8],pgamma[9])[1]
                                            #expt_feeding_gs.append(converted_primary)
                                            #d_expt_feeding_gs.append(d_converted_primary**2)

                            
                    if len(expt_feeding_gs) > 0:
                        sum_expt_feeding = sum(expt_feeding_gs)
                        d_sum_expt_feeding = np.sqrt(sum(d_expt_feeding_gs))

                        a = Analysis()
                        sigma_0 = a.modeled_sigma0(sum_expt_feeding, d_sum_expt_feeding, self.float1, self.float2)[0]
                        d_sigma_0 = a.modeled_sigma0(sum_expt_feeding, d_sum_expt_feeding, self.float1, self.float2)[1]
                        modeled_expt_cs.append([Ecrit, sum_expt_feeding, d_sum_expt_feeding, sigma_0, d_sigma_0])

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function using one of the below methods")
            print(" modeled_sigma0_ecrit(edata,<int>,<float>,<float>,\"Si29\")")
            print(" modeled_sigma0_ecrit(edata,<int>,<float>,<float>,14, 29)")
            return

        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return

        if LEVEL_ABOVE_ECRIT == True:
            print("Ecrit set too high!\nTry a lower value.")
            return
        
        if DIRECT_FEEDING_GS == False:
            print("No transitions feeding ground state directly.")
            return

        if len(modeled_expt_cs) > 0:
            #print("Ecrit = {0}".format(Ecrit))
            return modeled_expt_cs[0]
        
    
    def sum_feeding_gs(self,list,bool,*args,**kwargs):
        """Calculates the sum of all internal-conversion-corrected intensities
        (defined as isotopic partial gamma-ray cross sections or populations 
        per neutron capture) corresponding to all transitions that 
        feed the ground state of the compound nucleus in direct single-step 
        transitions.  For a ground state defined with level index 0, the first,
        second, third, etc. excited states are labeled 1, 2, and 3, 
        respectively, only those transitions with an associated final-level 
        index of 0 are included in the summation, i.e., 1->0, 2->0, 3->0, etc.  
        Primary gamma-ray transitions feeding the ground state directly (i.e., 
        capture state -> 0) may be included or removed from the summation 
        according to input arguments.

        Arguments:
            list: A list of EGAF-data JSON objects.
            bool: A boolean value to indicate whether the primary gamma-ray 
                  transition is included (`True`) or not (`False`) in the 
                  summation.
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
            A tuple containing the following elements:

            [0]: Summed internal-conversion-corrected isotopic partial gamma-ray
                 cross sections over all transitions feeding the ground state 
                 directly in single-step transitions (float);
            [1]: Associated uncertainty on the summed cross sections (float).

        Examples:
            For 56Fe(n,g)57Fe summed isotopic parttial gamma-ray cross sections
            including the primary:
            sum_feeding_gs(edata, True, "Fe57", intensity='iostopic')
            sum_feeding_gs(edata, True, 26, 57, intensity='isotopic')

            For 56Fe(n,g)57Fe summed isotopic parttial gamma-ray cross sections
            NOT including the primary:
            sum_feeding_gs(edata, False, "Fe57", intensity='iostopic')
            sum_feeding_gs(edata, False, 26, 57, intensity='isotopic')

            For 56Fe(n,g)57Fe summed populations per neutron capture:
            including the primary:
            sum_feeding_gs(edata, True, "Fe57", intensity='population')
            sum_feeding_gs(edata, True, 26, 57, intensity='population')

            For 56Fe(n,g)57Fe summed populations per neutron capture:
            NOT including the primary:
            sum_feeding_gs(edata, False, "Fe57", intensity='population')
            sum_feeding_gs(edata, False, 26, 57, intensity='population')
        """
        self.list = list
        self.bool = bool
        self.args = args
        UNSPECIFIED_INTENSITY = False
        DECAY_SCHEME_EXISTS = False
        DIRECT_FEEDING_GS = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True

        USE_PRIMARY = True
        if self.bool == False:
            USE_PRIMARY = False
        elif self.bool == True:
            USE_PRIMARY = True
            
        feeding_gs = []
        d_feeding_gs = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    DECAY_SCHEME_EXISTS = True
                    for each_l in jdict["levelScheme"]:
                        if each_l["numberOfGammas"] > 0:
                            for each_g in each_l["gammaDecay"]:
                                if USE_PRIMARY == True:
                                    if each_g["gammaFeedsGroundState"] == True:
                                        DIRECT_FEEDING_GS = True
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

                                        if alpha == None: alpha = 0.0
                                        if d_alpha == None: d_alpha = 0.0
                                        if cs == None: cs = 0.0
                                        if d_cs == None: d_cs = 0.0

                                        converted_cs = float(cs)*(1+float(alpha))
                                        d_converted_cs = Uncertainties.quad_error(self,converted_cs, cs, d_cs, (1+alpha), d_alpha)

                                        feeding_gs.append(converted_cs)
                                        d_feeding_gs.append(d_converted_cs)
                                elif USE_PRIMARY == False:
                                    if each_g["gammaFeedsGroundState"] == True and each_g["gammaTransitionType"] != "primary":
                                        # Don't include the primary
                                        DIRECT_FEEDING_GS = True
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

                                        if alpha == None: alpha = 0.0
                                        if d_alpha == None: d_alpha = 0.0
                                        if cs == None: cs = 0.0
                                        if d_cs == None: d_cs = 0.0

                                        converted_cs = float(cs)*(1+float(alpha))
                                        d_converted_cs = Uncertainties.quad_error(self,converted_cs, cs, d_cs, (1+alpha), d_alpha)

                                        feeding_gs.append(converted_cs)
                                        d_feeding_gs.append(d_converted_cs)

                                    
            except ValueError:
                WRONG_INPUTS = True
                                    
        tot_feeding_gs = sum(feeding_gs)
        d_tot_feeding_gs = np.sqrt(sum([dx**2 for dx in d_feeding_gs]))

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as, for example:")
            print(" sum_feeding_gs(edata, True, \"Fe57\", intensity=\"isotopic\")")
            print("or:")
            print(" sum_feeding_gs(edata, False, 26, 57, intensity=\"population\")")
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
            return (tot_feeding_gs, d_tot_feeding_gs)


    def sum_primaries(self,list,*args,**kwargs):
        """Calculates the sum of all internal-conversion-corrected intensities 
        from the associated primary gamma-ray transitions deexciting the capture
        state.  The intensity is defined according to the keyword argument 
        provided.

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
            A tuple containing the following elements:

            [0]: Summed internal-conversion-corrected gamma-ray intensities 
                 (isotopic partial gamma-ray cross sections or populations per 
                 neutron capture depending on keyword argument passed to the 
                 function) for all primary gamma-ray transitions (float);
            [1]: Associated uncertainty on the summed intensities (float).

        Examples:
            For 12C(n,g)13C summed primary partial gamma-ray isotopic cross 
            sections:
            sum_primaries(edata, "C13", intensity='isotopic')
            sum_primaries(edata, 6, 13, intensity='isotopic')

            For 28Si(n,g)29Si summed populations per neutron capture from the 
            associated primary gamma rays:
            sum_primaries(edata, "Si29", intensity='population')
            sum_primaries(edata, 14, 29, intensity='population')
        """
        self.list = list
        self.args = args
        UNSPECIFIED_INTENSITY = False
        HAS_PRIMARIES = False
        DECAY_SCHEME_EXISTS = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True
        
        primary_cs = []
        d_primary_cs = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    DECAY_SCHEME_EXISTS = True
                    for each_l in jdict["levelScheme"]:
                        if each_l["numberOfGammas"] > 0:
                            for each_g in each_l["gammaDecay"]:
                                if each_g["gammaTransitionType"] == "primary":
                                    HAS_PRIMARIES = True
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

                                    if alpha == None: alpha = 0.0
                                    if d_alpha == None: d_alpha = 0.0
                                    if cs == None: cs = 0.0
                                    if d_cs == None: d_cs = 0.0

                                    converted_cs = float(cs)*(1+float(alpha))
                                    d_converted_cs = Uncertainties.quad_error(self,converted_cs, cs, d_cs, (1+alpha), d_alpha)

                                    primary_cs.append(converted_cs)
                                    d_primary_cs.append(d_converted_cs)
            except ValueError:
                WRONG_INPUTS = True
                            
        tot_primary_cs = sum(primary_cs)
        d_tot_primary_cs = np.sqrt(sum([dx**2 for dx in d_primary_cs]))

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" sum_primaries(edata, \"C13\", intensity=<str>)")
            print("or:")
            print(" sum_primaries(edata, 6, 13, intensity=<str>)")
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
        
        if HAS_PRIMARIES == False:
            print("No primary gamma rays in compound-nucleus decay scheme.")
            return

        if primary_cs == []:
            print("No primary gamma rays in defined compound nucleus.")
            return
        else:
            return (tot_primary_cs, d_tot_primary_cs)
        

    def normalise_intensities(self,list,*args):
        """Calculates the sum of the level depopulation cross sections and 
        normalises the summed-level cross section to the total thermal-neutron 
        capture cross section.  These normalised level intensities are 
        effectively level depopulation intensities per neutron capture and can 
        be compared to the corresponding level population per neutron capture
        from statistical-model calculations.

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1, 2, 4, or 5 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID must be passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

                  (iii) 4 args:
                  residual: The residual ID must be passed as a string argument.
                  float1: Population per neutron capture feeding the ground 
                          state from the quasicontinuum according to a 
                          statistical-model calculation.
                  float2: Associated uncertainty on the calculated population 
                          per neutron capture.
                  int: Level index corresponding to the critical energy (i.e. 
                       cut-off energy).

                  (iv) 5 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.
                  float1: Population per neutron capture feeding the ground 
                          state from the quasicontinuum according to a 
                          statistical-model calculation.
                  float2: Associated uncertainty on the calculated population 
                          per neutron capture.
                  int: Level index corresponding to the critical energy (i.e. 
                       cut-off energy).

        Returns:
            A list containing the following elements:

            [0]: Level index (int);
            [1]: Associated level energy (float);
            [2]: Associated level energy uncertainty (float);
            [3]: Summed conversion-corrected total level-depopulation cross 
                 section (float);
            [4]: Associated uncertainty level-depopulation cross 
                 section (float);
            [5]: Normalised conversion-corrected total level-depopulation 
                 intensity (float);
            [6]: Associated uncertainty normalized level-depopulation 
                 intensity (float).
    
        Examples:
            To calculate normalised level-depopulation intensities using the 
            adopted total neutron-capture cross section:
            normalise_intensities(edata, "Si29")
            normalise_intensities(edata, 14, 29)

            To calculate normalised level-depopulation intensities using the 
            value of P0 obtained from a statistical-model calculation:
            normalise_intensities(edata, "Si29", 0.02217, 0.00051, 12)
            normalise_intensities(edata, 14, 29, 0.02217, 0.00051, 12)
        """
        self.list = list
        self.args = args
        DECAY_SCHEME_EXISTS = False
        UNSPECIFIED_UNIT = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 5:
            WRONG_INPUTS = True

        if type(args[-1]) != int:
            WRONG_INPUTS = True
            
        levels = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):

                    DECAY_SCHEME_EXISTS = True
                    
                    adopted_total_cs = None
                    d_adopted_total_cs = None

                    for each_n in jdict["neutronCaptureNormalization"]:
                        for each_r in each_n["normalizationRecord"]:

                            if each_r["unitAdoptedCrossSection"] == "b":
                                adopted_total_cs = each_r["adoptedTotalThermalCaptureCrossSection"]
                                d_adopted_total_cs = each_r["dAdoptedTotalThermalCaptureCrossSection"]
                            else:
                                UNSPECIFIED_UNIT = True
                                print("Adopted cross section units: {0}".format(each_r["unitAdoptedCrossSection"]))

                    for each_l in jdict["levelScheme"]:
                        level_index = each_l["levelIndex"]
                        level_energy = each_l["levelEnergy"]
                        d_level_energy = each_l["dLevelEnergy"]

                        sum_level_cs = 0
                        d_sum_level_cs = 0
                        if each_l["numberOfGammas"] > 0:                            
                            for each_g in each_l["gammaDecay"]:
                                alpha = each_g["calculatedTotalInternalConversionCoefficient"]
                                d_alpha = each_g["dCalculatedTotalInternalConversionCoefficient"]
                                for each_i in each_g["gammaAbsoluteIntensities"]:
                                    gamma_intensity = each_i["partialIsotopicCrossSection"]
                                    d_gamma_intensity = each_i["dPartialIsotopicCrossSection"]


                                    if alpha == None: alpha = 0.0
                                    if d_alpha == None: d_alpha = 0.0
                                    if gamma_intensity == None: gamma_intensity = 0.0
                                    if d_gamma_intensity == None: d_gamma_intensity = 0.0

                                    converted_cs = float(gamma_intensity)*(1.0+float(alpha))
                                    d_converted_cs = Uncertainties.quad_error(self,converted_cs, gamma_intensity, d_gamma_intensity, (1.0+alpha), d_alpha)

                                    sum_level_cs += converted_cs
                                    d_sum_level_cs += d_converted_cs**2

                        d_sum_level_cs = np.sqrt(d_sum_level_cs)

                        normalised_cs = sum_level_cs/adopted_total_cs
                        d_normalised_cs = Uncertainties.quad_error(self,normalised_cs, sum_level_cs, d_sum_level_cs, adopted_total_cs, d_adopted_total_cs)

                        levels.append([level_index, level_energy, d_level_energy, sum_level_cs, d_sum_level_cs, normalised_cs, d_normalised_cs])


                elif (len(args)==4 and str(args[0]) == jdict["nucleusID"]) or (len(args)==5 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):

                    DECAY_SCHEME_EXISTS = True

                    P0, dP0 = None, None
                    Ecrit = None
                    if len(args) == 4:
                        P0 = float(args[1])
                        dP0 = float(args[2])
                        Ecrit = int(args[3])
                        
                    elif len(args) == 5:
                        P0 = float(args[2])
                        dP0 = float(args[3])
                        Ecrit = int(args[4])

                    a = Analysis()
                    sigma_0 = a.modeled_sigma0_ecrit(self.list, Ecrit, P0, dP0, jdict["nucleusID"])[3]
                    d_sigma_0 = a.modeled_sigma0_ecrit(self.list, Ecrit, P0, dP0, jdict["nucleusID"])[4]

                    for each_l in jdict["levelScheme"]:
                        level_index = each_l["levelIndex"]
                        level_energy = each_l["levelEnergy"]
                        d_level_energy = each_l["dLevelEnergy"]

                        sum_level_cs = 0
                        d_sum_level_cs = 0
                        if each_l["numberOfGammas"] > 0:                            
                            for each_g in each_l["gammaDecay"]:
                                alpha = each_g["calculatedTotalInternalConversionCoefficient"]
                                d_alpha = each_g["dCalculatedTotalInternalConversionCoefficient"]
                                for each_i in each_g["gammaAbsoluteIntensities"]:
                                    gamma_intensity = each_i["partialIsotopicCrossSection"]
                                    d_gamma_intensity = each_i["dPartialIsotopicCrossSection"]

                                    if alpha == None: alpha = 0.0
                                    if d_alpha == None: d_alpha = 0.0
                                    if gamma_intensity == None: gamma_intensity = 0.0
                                    if d_gamma_intensity == None: d_gamma_intensity = 0.0

                                    converted_cs = float(gamma_intensity)*(1.0+float(alpha))
                                    d_converted_cs = Uncertainties.quad_error(self,converted_cs, gamma_intensity, d_gamma_intensity, (1.0+alpha), d_alpha)

                                    sum_level_cs += converted_cs
                                    d_sum_level_cs += d_converted_cs**2

                        d_sum_level_cs = np.sqrt(d_sum_level_cs)

                        #sigma_0 = a.sum_feeding_gs(self.list, False, jdict["nucleusID"],intensity='isotopic')[0]/(1-P0)
                        #d_sigma_0 = Uncertainties.quad_error(self, sigma_0, a.sum_feeding_gs(self.list, False, jdict["nucleusID"],intensity='isotopic')[0], a.sum_feeding_gs(self.list, False, jdict["nucleusID"],intensity='isotopic')[1], P0, dP0)
                        
                        normalised_cs = sum_level_cs/sigma_0
                        d_normalised_cs = Uncertainties.quad_error(self,normalised_cs, sum_level_cs, d_sum_level_cs, sigma_0, d_sigma_0)

                        levels.append([level_index, level_energy, d_level_energy, sum_level_cs, d_sum_level_cs, normalised_cs, d_normalised_cs])
                
                        
            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function using one of the below methods:")
            print(" normalise_intensities(edata, \"C13\")")
            print(" normalise_intensities(edata, 6, 13)")
            print(" normalise_intensities(edata, \"C13\", <float>, <float>, <int>)")
            print(" normalise_intensities(edata, 6, 13, <float>, <float>, <int>)")
            return

        if UNSPECIFIED_UNIT == True:
            print("Not currently handling the adopted cross section units.")
            return

        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return

        if len(levels) > 0:
            return levels
        else:
            print("No depopulation data available.")
            return
        

    def level_depopulations(self,list,*args):
        """Calculates the sum of the level depopulation intensities.  These 
        normalised level intensities are effectively level depopulation 
        intensities per neutron capture and can be compared to the 
        corresponding level population per neutron capture from 
        statistical-model calculations.

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
            A list containing the following elements:

            [0]: Level index (int);
            [1]: Associated level energy (float);
            [2]: Associated level energy uncertainty (float);
            [3]: Total level-depopulation intensity (float);
            [3]: Associated uncertainty level-depopulation intensity (float).

        Examples:
            level_depopulations(edata, "C13")
            level_depopulations(edata, 6, 13)
        """
        self.list = list
        self.args = args
        DECAY_SCHEME_EXISTS = False
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 4:
            WRONG_INPUTS = True

        levels = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):

                    DECAY_SCHEME_EXISTS = True
                    
                    for each_l in jdict["levelScheme"]:
                        level_index = each_l["levelIndex"]
                        level_energy = each_l["levelEnergy"]
                        d_level_energy = each_l["dLevelEnergy"]

                        sum_level_depop = 0
                        d_sum_level_depop = 0
                        if each_l["numberOfGammas"] > 0:                            
                            for each_g in each_l["gammaDecay"]:
                                alpha = each_g["calculatedTotalInternalConversionCoefficient"]
                                d_alpha = each_g["dCalculatedTotalInternalConversionCoefficient"]
                                for each_i in each_g["gammaAbsoluteIntensities"]:
                                    gamma_intensity = each_i["populationPerNeutronCapture"]
                                    d_gamma_intensity = each_i["dPopulationPerNeutronCapture"]

                                    if alpha == None: alpha = 0.0
                                    if d_alpha == None: d_alpha = 0.0
                                    if gamma_intensity == None: gamma_intensity = 0.0
                                    if d_gamma_intensity == None: d_gamma_intensity = 0.0

                                    converted_depop = float(gamma_intensity)*(1.0+float(alpha))
                                    d_converted_depop = Uncertainties.quad_error(self,converted_depop, gamma_intensity, d_gamma_intensity, (1.0+alpha), d_alpha)

                                    sum_level_depop += converted_depop
                                    d_sum_level_depop += d_converted_depop**2
                                    
                        d_sum_level_depop = np.sqrt(d_sum_level_depop)

                        levels.append([level_index, level_energy, d_level_energy, sum_level_depop, d_sum_level_depop])

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" level_depopulations(edata, \"C13\")")
            print(" level_depopulations(edata, 6, 13)")
            return

        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return

        if len(levels) > 0:
            return levels
        else:
            print("No depopulation data available.")
            return


    def intensity_balance(self,list,*args,**kwargs):
        """Gamma-ray intensity balance corrected for internal conversion for 
        all transitions feeding and deexciting each level in the residual 
        compound nucleus.

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

        
        Returns:
            A numpy array containing the following elements associated with the
            intensity balance for each level in the residual compound nucleus:

            [0]: Level index (float);
            [1]: Level energy in keV (float);
            [2]: Total internal-conversion corrected gamma-ray intensity 
                 deexciting the level (float);
            [3]: Associated uncertainty on the total-level deexcitation 
                 intensity (float);
            [4]: Total internal-conversion corrected gamma-ray intensity 
                 feeding the level (float);
            [5]: Associated uncertainty on the total-level feeding 
                 intensity (float);
            [6]: Intensity difference for each level given by  
                 `population - depopulation` (float);
            [7]: Associated uncertainty on the intensity difference (float);
            [8]: Residual intensity difference for each level in units of 
                 standard deviation `sigma` (float).

        Examples:
            Intesnity balance for isotopic partial gamma-ray cross sections:
            intensity_balance(edata, "S35", intensity="isotopic") 
            intensity_balance(edata, 16, 35, intensity="isotopic")
 
            Intensity balance for elemental partial gamma-ray cross sections:
            intensity_balance(edata, "Si29", intensity="elemental") 
            intensity_balance(edata, 14, 29, intensity="elemental")
        """
        self.list = list
        self.args = args
        DECAY_SCHEME_EXISTS = False
        UNSPECIFIED_INTENSITY = True
        intensity_units = None
        WRONG_INPUTS = False
        CAPTURE_STATE = False
        highest_level_if_no_cs = None
        if len(args) == 0 or len(args)> 2:
            WRONG_INPUTS = True
        
        if kwargs == {} or kwargs == None:
            print("A keyword argument is required for the desired intensity units.")
            print("Please pass one of the following arguments:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            return
        else:
            for intensity in kwargs.values():
                if intensity.lower() == str("elemental"):
                    UNSPECIFIED_INTENSITY = False
                    intensity_units = "ELEMENTAL"
                elif intensity.lower() == str("isotopic"):
                    UNSPECIFIED_INTENSITY = False
                    intensity_units = "ISOTOPIC"
                else:
                    UNSPECIFIED_INTENSITY = True

        if UNSPECIFIED_INTENSITY == True:
            print("Incorrect intensity argument.")
            print("Please pass one of the following keyword arguments:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            return
        
        intensity_balance = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):        

                    DECAY_SCHEME_EXISTS = True
                    target = jdict["nucleusTargetID"]
                    residual = jdict["nucleusID"]

                    for eachq in jdict["recordQ"]:
                        if eachq["energyNeutronSeparationEGAF"] == None:
                            CAPTURE_STATE = False
                        else:
                            try:
                                if float(eachq["energyNeutronSeparationEGAF"]) > 0.0:
                                    CAPTURE_STATE = True
                            except ValueError:
                                print("Check Sn data type in {0}(n,g){1} dataset".format(target, residual))
                                raise
                            
                    g = Gammas()
                    gammas_array = None
                    if intensity_units == "ELEMENTAL":
                        gammas_array = g.get_gammas(self.list, residual, intensity="elemental")
                    elif intensity_units == "ISOTOPIC":
                        gammas_array = g.get_gammas(self.list, residual, intensity="isotopic")
                    else:
                        break

                    gammas = gammas_array.tolist()
                    #max_level = max(list(set([int(gamma[0]) for gamma in gammas])))
                    max_level = max(set([int(gamma[0]) for gamma in gammas]))

                    if CAPTURE_STATE == False:
                        highest_level_if_no_cs = set([gamma[2] for gamma in gammas if int(gamma[0]) == max_level])

                    for level in range(0,max_level+1,1):
                        level_energy = None
                        level_depop = 0.0
                        d_level_depop = 0.0
                        LEVEL_IS_DEPOPULATED = False
                        for gamma in gammas:
                            if int(gamma[0]) == level:
                                LEVEL_IS_DEPOPULATED = True
                                level_energy = float(gamma[2])
                                depop = float(gamma[6]) * (1+float(gamma[8]))
                                level_depop += depop
            
                                #f,x,dx,y,dy
                                d_depop = Uncertainties().quad_error(depop, float(gamma[6]), float(gamma[7]), (1+float(gamma[8])), float(gamma[9]))
                                d_level_depop += float(d_depop)**2
                        d_level_depop = np.sqrt(d_level_depop)

                        level_pop = 0.0
                        d_level_pop = 0.0
                        for gamma in gammas:
                            if int(gamma[1]) == level:
                                pop = float(gamma[6]) * (1+float(gamma[8]))
                                level_pop += pop
                                if LEVEL_IS_DEPOPULATED == False:
                                    level_energy = float(gamma[3])
            
                                #f,x,dx,y,dy
                                d_pop = Uncertainties().quad_error(pop, float(gamma[6]), float(gamma[7]), (1+float(gamma[8])), float(gamma[9]))
                                d_level_pop += float(d_pop)**2
                        d_level_pop = np.sqrt(d_level_pop)

                        diff = level_pop - level_depop
                        d_diff = np.sqrt(d_level_depop**2 + d_level_pop**2)
                        if d_diff > 0:
                            res = diff / d_diff
                            if level > 0 and level < max_level:
                                intensity_balance.append([level, level_energy, level_depop, d_level_depop, level_pop, d_level_pop, diff, d_diff, res])

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" intensity_balance(edata, \"S35\", intensity=<str>)")
            print("or:")
            print(" intensity_balance(edata, 16, 35, intensity=<str>)")
            return

        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return
        if len(intensity_balance) == 0:
            print("No gammas in decay scheme")
            return
        else:
            if CAPTURE_STATE == True:
                return np.array(intensity_balance)
            else:
                print("Warning! Neutron-capture state not measured in EGAF.")
                print("Highest level observed in residual nucleus: {0} keV.".format(highest_level_if_no_cs))
                return np.array(intensity_balance)


    def dead_ends(self,list,*args,**kwargs):
        """Finds difference between total internal-conversion corrected 
        intensity deexciting the neutron-capture state and the total internal-
        conversion corrected intensity feeding the ground state.

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

        Returns:
            A list containing the following elements associated with the 
            intensities of the neutron capture state and ground state in the 
            residual compound nucleus:

            [0]: Level index of the ground state (int);
            [1]: Level energy of the ground state in keV (float);
            [2]: Total internal-conversion corrected gamma-ray intensity 
                 feeding the ground state (float);
            [3]: Associated uncertainty on the total-feeding intensity of the 
                 ground state (float);
            [4]: Level index of the neutron-capture state (int);
            [5]: Level energy of the neutron-capture state in keV (float);
            [6]: Total internal-conversion corrected gamma-ray intensity 
                 deexciting the neutron-capture state (float);
            [7]: Associated uncertainty on the total-deexcitation intensity of 
                 the neutron-capture state (float);
            [8]: Intensity difference: 
                 `(ground state intensity) - (capture state intensity)` 
                 This is a signed value depending on missing (negative) or 
                 excess (positive) observed intensity feeding the ground 
                 state (float);
            [9]: Associated uncertainty on the intensity difference (float); 
            [10]: Residual intensity difference between ground and capture 
                  states in units of standard deviation `sigma` (float).
            [11]: Percentage intensity difference relative to the total 
                  deexcitation of the neutron-capture state.  This is a signed 
                  value depending on missing (negative) or excess (positive) 
                  observed intensity feeding the ground state (float).

        Examples:
            Missing intensity to the ground state in 34S(n,g)35S
            dead_ends(edata, "S35", intensity="isotopic")
            dead_ends(edata, 16, 35, intensity="isotopic")

            Excess intensity observed feeding the ground state in 186W(n,g)187W
            dead_ends(edata, "W187", intensity="elemental")
            dead_ends(edata, 74, 187, intensity="elemental")
        """
        self.list = list
        self.args = args
        DECAY_SCHEME_EXISTS = False
        UNSPECIFIED_INTENSITY = True
        intensity_units = None
        WRONG_INPUTS = False
        CAPTURE_STATE = False
        highest_level_if_no_cs = None
        if len(args) == 0 or len(args)> 2:
            WRONG_INPUTS = True
        
        if kwargs == {} or kwargs == None:
            print("A keyword argument is required for the desired intensity units.")
            print("Please pass one of the following arguments:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            return
        else:
            for intensity in kwargs.values():
                if intensity.lower() == str("elemental"):
                    UNSPECIFIED_INTENSITY = False
                    intensity_units = "ELEMENTAL"
                elif intensity.lower() == str("isotopic"):
                    UNSPECIFIED_INTENSITY = False
                    intensity_units = "ISOTOPIC"
                else:
                    UNSPECIFIED_INTENSITY = True

        if UNSPECIFIED_INTENSITY == True:
            print("Incorrect intensity argument.")
            print("Please pass one of the following keyword arguments:")
            print("intensity='elemental'")
            print("intensity='isotopic'")
            return
        
        deadends = []
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):        

                    DECAY_SCHEME_EXISTS = True
                    target = jdict["nucleusTargetID"]
                    residual = jdict["nucleusID"]

                    for eachq in jdict["recordQ"]:
                        if eachq["energyNeutronSeparationEGAF"] == None:
                            CAPTURE_STATE = False
                        else:
                            try:
                                if float(eachq["energyNeutronSeparationEGAF"]) > 0.0:
                                    CAPTURE_STATE = True
                            except ValueError:
                                print("Check Sn data type in {0}(n,g){1} dataset".format(target, residual))
                                raise
                    
                    g = Gammas()
                    gammas_array = None
                    if intensity_units == "ELEMENTAL":
                        gammas_array = g.get_gammas(self.list, residual, intensity="elemental")
                    elif intensity_units == "ISOTOPIC":
                        gammas_array = g.get_gammas(self.list, residual, intensity="isotopic")
                    else:
                        break

                    gammas = gammas_array.tolist()
                    level_index_cs = max(set([int(gamma[0]) for gamma in gammas]))
                    level_energy_cs = None

                    if CAPTURE_STATE == False:
                        highest_level_if_no_cs = set([gamma[2] for gamma in gammas if int(gamma[0]) == level_index_cs])
                    
                    level_depop_cs = 0.0
                    d_level_depop_cs = 0.0                    
                    for gamma in gammas:
                        if int(gamma[0]) == level_index_cs:
                            level_energy_cs = float(gamma[2])
                            depop_cs = float(gamma[6]) * (1+float(gamma[8]))
                            level_depop_cs += depop_cs

                            #f,x,dx,y,dy
                            d_depop_cs = Uncertainties().quad_error(depop_cs, float(gamma[6]), float(gamma[7]), (1+float(gamma[8])), float(gamma[9]))
                            d_level_depop_cs += float(d_depop_cs)**2
                    d_level_depop_cs = np.sqrt(d_level_depop_cs)

                    level_index_gs = 0
                    level_energy_gs = None
                    level_pop_gs = 0.0
                    d_level_pop_gs = 0.0
                    GROUND_STATE_FEEDING = False
                    for gamma in gammas:
                        if int(gamma[1]) == level_index_gs:
                            GROUND_STATE_FEEDING = True
                            level_energy_gs = float(gamma[3])
                            pop_gs = float(gamma[6]) * (1+float(gamma[8]))
                            level_pop_gs += pop_gs

                            #f,x,dx,y,dy
                            d_pop_gs = Uncertainties().quad_error(pop_gs, float(gamma[6]), float(gamma[7]), (1+float(gamma[8])), float(gamma[9]))
                            d_level_pop_gs += float(d_pop_gs)**2
                    d_level_pop_gs = np.sqrt(d_level_pop_gs)

                    if GROUND_STATE_FEEDING == False:
                        level_energy_gs = 0.0

                    diff = level_pop_gs - level_depop_cs
                    d_diff = np.sqrt(d_level_depop_cs**2 + d_level_pop_gs**2)

                    percent_diff = (diff/level_depop_cs)*100
                    d_percent_diff = (d_diff/level_depop_cs)*100 # Not included in list output

                    if d_diff > 0:
                        res = diff / d_diff

                        deadends.append([level_index_gs, level_energy_gs, level_pop_gs, d_level_pop_gs, level_index_cs, level_energy_cs, level_depop_cs, d_level_depop_cs, diff, d_diff, res, percent_diff])

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" intensity_balance(edata, \"S35\", intensity=<str>)")
            print("or:")
            print(" intensity_balance(edata, 16, 35, intensity=<str>)")
            return

        if DECAY_SCHEME_EXISTS == False:
            print("No residual compound-nucleus decay scheme in EGAF for input arguments provided.")
            return
        if len(deadends) == 0:
            print("No gammas in decay scheme")
            return
        else:
            if CAPTURE_STATE == True:
                return deadends[0]
            else:
                print("Warning! Neutron-capture state not measured in EGAF.")
                print("Highest level observed in residual nucleus: {0} keV.".format(highest_level_if_no_cs))
                return deadends[0]
