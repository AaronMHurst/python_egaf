from .base_egaf import *

class Separation(Uncertainties):    
    __doc__="""Class to handle neutron and proton separation energies."""

    def __init__(self):
        BaseEGAF.__init__(self)
        Meta.__init__(self)
        Uncertainties.__init__(self)

    def get_residual_Sn_AME(self,list,*args):
        """Neutron-separation energy from the AME2020 Atomic Mass Evaluation.  
        All values are taken from Ref:

        [2021Wa16]: M. Wang et al., Chin. Phys. C 45, 030003 (2021).

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID of the compound nucleus produced in 
                            the (n,g) reaction passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns: 
            A tuple object with the following elements:

            [0]: Neutron-separation energy from AME2020 (float);
            [1]: Neutron-separation energy associated uncertainty (float).

        Examples:
            For 28Al produced in 27Al(n,g):
            get_residual_Sn_AME(edata,"Al28")
            get_residual_Sn_AME(edata,13,28)
        """
        self.list = list
        self.args = args
        Sn_AME = None
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True        
        for jdict in self.list:
            try:
                if (len(self.args) == 1 and str(args[0]) == jdict["nucleusID"]) or (len(self.args) == 2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    for each_q in jdict["recordQ"]:
                        Sn_AME = (each_q["energyNeutronSeparationAME2020"],each_q["dEnergyNeutronSeparationAME2020"])
            except ValueError:
                WRONG_INPUTS = True
        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_residual_Sn_AME(edata,\"Al28\")")
            print("or:")
            print(" get_residual_Sn_AME(edata,13,28)")
            return                 
        if Sn_AME == None:
            print("No residual nucleus in EGAF file for defined input.")
            return None
        else:
            return Sn_AME

    def get_residual_Sp_AME(self,list,*args):
        """Proton-separation energy from the AME2020 Atomic Mass Evaluation.  
        All values are taken from Ref:

        [2021Wa16]: M. Wang et al., Chin. Phys. C 45, 030003 (2021).

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID of the compound nucleus produced in 
                            the (n,g) reaction passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns: 
            A tuple object with the following elements:

            [0]: Proton-separation energy from AME2020 (float);
            [1]: Proton-separation energy associated uncertainty (float).

        Examples:
            For 28Al produced in 27Al(n,g):
            get_residual_Sp_AME(edata,"Al28")
            get_residual_Sp_AME(edata,13,28)
        """
        self.list = list
        self.args = args
        Sp_AME = None
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True                
        for jdict in self.list:
            try:
                if (len(self.args) == 1 and str(args[0]) == jdict["nucleusID"]) or (len(self.args) == 2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    for each_q in jdict["recordQ"]:
                        Sp_AME = (each_q["energyProtonSeparationAME2020"],each_q["dEnergyProtonSeparationAME2020"])
            except ValueError:
                WRONG_INPUTS = True
                
        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_residual_Sp_AME(edata,\"Al28\")")
            print("or:")
            print(" get_residual_Sp_AME(edata,13,28)")
            return
        if Sp_AME == None:
            print("No residual nucleus in EGAF file for defined input.")
            return None
        else:
            return Sp_AME        

    def get_residual_Sn_EGAF(self,list,*args):
        """Neutron-separation energy from EGAF.  Provided that the residual 
        compound nucleus produced in the (n,g) reaction contains primary gamma 
        rays, the neutron-separation energy will correspond to the 
        neutron-capture state and will be represented by  the highest-value 
        energy level in the decay scheme.

        Arguments:
            list: A list of EGAF-data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: The residual ID of the compound nucleus produced in 
                            the (n,g) reaction passed as a string argument.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns: 
            A tuple object with the following elements:

            [0]: Neutron-separation energy from AME2020 (float);
            [1]: Neutron-separation energy associated uncertainty (float).

        Examples:
            For 28Al produced in 27Al(n,g):
            get_residual_Sn_EGAF(edata,"Al28")
            get_residual_Sn_EGAF(edata,13,28)
        """
        self.list = list
        self.args = args
        Sn_EGAF = None
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True        
        for jdict in self.list:
            try:
                if (len(self.args) == 1 and str(args[0]) == jdict["nucleusID"]) or (len(self.args) == 2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    for each_q in jdict["recordQ"]:
                        Sn_EGAF = (each_q["energyNeutronSeparationEGAF"],each_q["dEnergyNeutronSeparationEGAF"])
            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" get_residual_Sn_EGAF(edata,\"Al28\")")
            print("or:")
            print(" get_residual_Sn_EGAF(edata,13,28)")
            return                        
        if Sn_EGAF == None:
            print("No residual nucleus in EGAF file for defined input.")
            return None
        else:
            return Sn_EGAF
        
    def get_all_separation_energies(self,list,str):
        """All separation energies according to user-specified input whether 
        values from AME2020 or EGAF are required.  The EGAF values correspond to
        the reported neutron-capture state and the AME2020 values are taken 
        from Ref:

        [2021Wa16]: M. Wang et al., Chin. Phys. C 45, 030003 (2021).

        Arguments:
            list: A list of EGAF-data JSON objects.
            str: One of the following case-insensitive string arguments is 
                 passed according to required output:

                 "proton"  : Proton-separation energies from AME2020.
                 "neutron" : Proton-separation energies from AME2020.
                 "egaf"    : Neutron-separation energies from EGAF.
        
        Returns:
            A dictionary object.  The key is a string object associated with the
            compound nucleus produced in the (n,g) reaction.  The value is a 
            two-element tuple containing the input-specified separation energy 
            and associated uncertainty:

            [0]: Proton (AME2020) or Neutron (AME2020 or EGAF) separation 
                 energy (float);
            [1]: Associated separation energy uncertainty (float).
        
        Examples:
            For all Sp from AME2020: 
            get_all_separation_energies(edata, "proton")
        
            For all Sn from AME2020: 
            get_all_separation_energies(edata, "neutron")

            For all Sn EGAF: 
            get_all_separation_energies(edata, "egaf")
        """
        self.str = str
        self.list = list
        egaf_residuals = BaseEGAF.egaf_residual_list(self,self.list)
        separation_dict = {}
        if self.str.lower()=="neutron":
            for residual in egaf_residuals:
                SnAME = Separation.get_residual_Sn_AME(self,self.list,residual)
                separation_dict.update({residual: (SnAME[0], SnAME[1])})
        elif self.str.lower()=="proton":
            for residual in egaf_residuals:
                SpAME = Separation.get_residual_Sp_AME(self,self.list,residual)
                separation_dict.update({residual: (SpAME[0], SpAME[1])})
        elif self.str.lower()=="egaf":
            for residual in egaf_residuals:
                SnEGAF = Separation.get_residual_Sn_EGAF(self,self.list,residual)
                separation_dict.update({residual: (SnEGAF[0], SnEGAF[1])})
        else:
            print("Parameter passed has no return value.")
            print("Acceptable strings are: 'neutron', 'proton', or 'egaf'")
            separation_dict = None

        if separation_dict == {} or separation_dict == None:
            print("No separation energy data available.")
            return
        else:
            return separation_dict
        
