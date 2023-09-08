import numpy as np
import pandas as pd
import json
import glob
import re
import os

class BaseEGAF(object):
    __doc__="""Base class to handle EGAF data sets."""

    def __init__(self,contents=None):
        self.contents = [] or None

    def load_egaf(self):
        """Function to assign all 245 JSON-formatted EGAF thermal neutron 
        capture (n,g) data sets to a list object variable.
        
        Arguments:
            No arguments are passed to this function.
        
        Returns:
            A list object containing all 3226 JSON-formatted ENSDF-decay data 
            sets.

        Example:
            
            import pyEGAF as egaf
            e = egaf.EGAF()
            edata = e.load_egaf()
        """
        print("Loading EGAF data sets, please wait...")
        
        from . import get_data
        EGAF_JSON_PATH = get_data('EGAF_JSON')
        json_egaf_list = [j for j in glob.glob("%s/*.json"%EGAF_JSON_PATH)]
        json_egaf_data = []

        JSON_COUNT = 0
        for json_file in json_egaf_list:
            JSON_COUNT += 1
            with open(json_file, mode='r') as jf:
                json_egaf_dict = json.loads(jf.read())
                json_egaf_data.append(json_egaf_dict)
            jf.close()

        if JSON_COUNT == 245:
            print("All {0} JSON-formatted EGAF data sets loaded.".format(JSON_COUNT))
        elif (JSON_COUNT > 0) and (JSON_COUNT < 245):
            print("{0} JSON-formatted EGAF data sets loaded.".format(JSON_COUNT))
            print("{0} JSON-formatted EGAF data sets are missing.".format(245-int(JSON_COUNT)))
        else:
            if JSON_COUNT == 0:
                print("{0} JSON-formatted EGAF data sets loaded.".format(JSON_COUNT))
        return json_egaf_data
    
    def sort_by_json_key(list,str='nucleusTargetID'):
        """Internal function: Sorts list in alphabetical order of target 
        nucleus ID."""
        return list['%s'%str]

    def sort_by_residual_A(list,str="nucleusA"):
        """Internal function: Sorts list in order of increasing mass (A) of 
        residual nucleus."""
        return list['%s'%str]

    def egaf_target_list(self,list):
        """Target nuclides for all (n,g) data sets.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.

        Returns: 
            A list of alphabetically-sorted target IDs for each (n,g) data set.
            Each target ID is a string object.

        Example:
            egaf_target_list(edata) 
        """
        self.list = list
        egaf_targets = sorted(self.list, key=BaseEGAF.sort_by_json_key)
        list_of_targets = []
        for jdict in egaf_targets:
            list_of_targets.append(jdict['nucleusTargetID'])
        return list_of_targets

    def egaf_residual_list(self,list):
        """Residual nuclides (A+1) corresponding to each target (A) (n,g) data 
        set.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.

        Returns: 
            A list of alphabetically-sorted residual (A+1) nuclides produced 
            from each target (A) (n,g) data set.  Each residual ID is a string 
            object.

        Example:
            egaf_residual_list(edata) 
        """
        self.list = list
        egaf_targets = sorted(self.list, key=BaseEGAF.sort_by_json_key)
        list_of_residuals = []
        for jdict in egaf_targets:
            list_of_residuals.append(jdict['nucleusID'])
        return list_of_residuals

    def egaf_target_residual_dict(self,list):
        """Residual (A+1) - target (A) pairs for each EGAF (n,g) data set.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.

        Returns: 
            A dictionary object with a key containing the residual (A+1) 
            nucleus ID in a string.  The dictionary value is given in a tuple 
            with the following elements:

            [0]: Target (A) nucleus ID <str>;
            [1]: Target atomic (Z) number <int>;
            [2]: Mass (A) number of residual compound nucleus <int>.


        Example:
            egaf_target_residual_dict(edata)
        """
        self.list = list
        egaf_residuals = sorted(self.list, key=BaseEGAF.sort_by_residual_A)
        dict_of_residuals = {}
        for jdict in egaf_residuals:
            dict_of_residuals.update({jdict['nucleusID']:(jdict['nucleusTargetID'], jdict['nucleusZ'], jdict['nucleusA'])})
        return dict_of_residuals

    
class Meta(BaseEGAF):
    __doc__="""Class to handle the auxiliary meta data in the EGAF file."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def num_primaries(self,list,*args):
        """Finds number of primary gammas (i.e., gamma rays that originate at 
        the neutron-capture state) in the residual decay scheme of a compound 
        nucleus produced in a thermal (n,g) reaction.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: A string object corresponding to the residual ID of 
                            the (A+1) compound nucleus produced in an (n,g) 
                            reaction.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns:
            An integer object corresponding to the number of primary gamma rays 
            in the residual decay scheme of the compound nucleus.

        Example:
            For the 28Si(n,g)29Si reaction:

            num_primaries(edata,"Si29")
            num_primaries(edata,14,29)
        """
        self.list = list
        self.args = args
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True
        
        num_primaries_gammas = None
        target = None
        residual = None
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    residual = jdict["nucleusID"]
                    target = jdict["nucleusTargetID"]
                    num_primaries_gammas = jdict["numberPrimaryGammas"]

                    print("{0}(n,g){1}".format(target,residual))
                    print("Number of primaries = {0}".format(num_primaries_gammas))
            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" num_primaries(edata,\"Cl36\")")
            print("or:")
            print(" num_primaries(edata,17,36)")
            return
                
        if num_primaries_gammas != None:
            return num_primaries_gammas
        else:
            print("No (n,g) data for input residual nucleus.")
            return
        
    def num_secondaries(self,list,*args):
        """Finds number of secondary gammas (i.e., gamma rays associated with   
        transitions between low-lying levels) in the residual decay scheme of a 
        compound nucleus produced in a thermal (n,g) reaction.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: A string object corresponding to the residual ID of 
                            the (A+1) compound nucleus produced in an (n,g) 
                            reaction.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns:
            An integer object corresponding to the number of secondary gamma 
            rays in the residual decay scheme of the compound nucleus.

        Example:
            For the 28Si(n,g)29Si reaction:

            num_secondaries(edata,"Si29")
            num_secondaries(edata,14,29)
        """
        self.list = list
        self.args = args
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True        

        num_secondaries_gammas = None
        target = None
        residual = None
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    residual = jdict["nucleusID"]
                    target = jdict["nucleusTargetID"]
                    num_secondaries_gammas = jdict["numberSecondaryGammas"]

                    print("{0}(n,g){1}".format(target,residual))
                    print("Number of secondaries = {0}".format(num_secondaries_gammas))
            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" num_secondaries(edata,\"Cl36\")")
            print("or:")
            print(" num_secondaries(edata,17,36)")
            return
                
        if num_secondaries_gammas != None:
            return num_secondaries_gammas
        else:
            print("No (n,g) data for input residual nucleus.")
            return

    def num_gammas(self,list,*args):
        """Finds total number of gammas in the residual decay scheme of a 
        compound nucleus produced in a thermal (n,g) reaction.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: A string object corresponding to the residual ID of 
                            the (A+1) compound nucleus produced in an (n,g) 
                            reaction.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns:
            An integer object corresponding to the total number of gamma rays 
            in the residual decay scheme of the compound nucleus.

        Example:
            For the 28Si(n,g)29Si reaction:

            num_gammas(edata,"Si29")
            num_gammas(edata,14,29)
        """
        self.list = list
        self.args = args
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True        

        number_gammas = None
        target = None
        residual = None
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    residual = jdict["nucleusID"]
                    target = jdict["nucleusTargetID"]
                    number_gammas = jdict["totalNumberGammas"]

                    print("{0}(n,g){1}".format(target,residual))
                    print("Total number of gammas = {0}".format(number_gammas))

            except ValueError:
                WRONG_INPUTS = True
                
        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" num_gammas(edata,\"Cl36\")")
            print("or:")
            print(" num_gammas(edata,17,36)")
            return
                
        if number_gammas != None:
            return number_gammas
        else:
            print("No (n,g) data for input residual nucleus.")
            return

    def num_levels(self,list,*args):
        """Finds total number of levels in the residual decay scheme of a 
        compound nucleus produced in a thermal (n,g) reaction.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: A string object corresponding to the residual ID of 
                            the (A+1) compound nucleus produced in an (n,g) 
                            reaction.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns:
            An integer object corresponding to the total number of levels in  
            the residual decay scheme of the compound nucleus.

        Example:
            For the 28Si(n,g)29Si reaction:

            num_levels(edata,"Si29")
            num_levels(edata,14,29)
        """
        self.list = list
        self.args = args
        WRONG_INPUTS = False
        if len(args) == 0 or len(args) > 2:
            WRONG_INPUTS = True        

        number_levels = None
        target = None
        residual = None
        for jdict in self.list:
            try:
                if (len(args)==1 and str(args[0]) == jdict["nucleusID"]) or (len(args)==2 and int(args[0]) == jdict["nucleusZ"] and int(args[1]) == jdict["nucleusA"]):
                    residual = jdict["nucleusID"]
                    target = jdict["nucleusTargetID"]
                    number_levels = jdict["totalNumberLevels"]

                    print("{0}(n,g){1}".format(target,residual))
                    print("Total number of levels = {0}".format(number_levels))

            except ValueError:
                WRONG_INPUTS = True

        if WRONG_INPUTS == True:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" num_levels(edata,\"Cl36\")")
            print("or:")
            print(" num_levels(edata,17,36)")
            return
                
        if number_levels != None:
            return number_levels
        else:
            print("No (n,g) data for input residual nucleus.")
            return

    def get_stats(self,list,*args):
        """Decay-scheme statistics from EGAF associated with the (A+1) residual 
        compound nucleus produced in a thermal (n,g) reaction.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            args: Takes either 1 or 2 additional arguments:
            
                  (i) 1 args:
                  residual: A string object corresponding to the residual ID of 
                            the (A+1) compound nucleus produced in an (n,g) 
                            reaction.

                  (ii) 2 args:
                  Z: Atomic number passed as an integer argument.
                  A: Atomic mass of the residual compound nucleus passed as an 
                     integer argument.

        Returns:
            An list object with the following elements:

            [0]: Number of primary gammas in (A+1) compound nucleus (int);
            [1]: Number of secondary gammas in (A+1) compound nucleus (int);
            [2]: Total number of gammas in (A+1) compound nucleus (int);
            [3]: Number of levels in (A+1) compound nucleus (int);

        Example:
            For the 28Si(n,g)29Si reaction:

            get_stats(edata,"Si29")
            get_stats(edata,14,29)
        """
        self.list = list
        if len(args) == 1:
            p=self.num_primaries(self.list,args[0])
            s=self.num_secondaries(self.list,args[0])
            g=self.num_gammas(self.list,args[0])
            l=self.num_levels(self.list,args[0])
            meta_data = [p,s,g,l]
            if p == None and s == None and g == None and l == None:
                return
            else:
                return meta_data
        elif len(args) == 2:
            p=self.num_primaries(self.list,args[0],args[1])
            s=self.num_secondaries(self.list,args[0],args[1])
            g=self.num_gammas(self.list,args[0],args[1])
            l=self.num_levels(self.list,args[0],args[1])
            meta_data = [p,s,g,l]
            if p == None and s == None and g == None and l == None:
                return
            else:
                return meta_data
        else:
            print("Incorrect input sequence.")
            print("Pass arguments to function as:")
            print(" num_levels(edata,\"Cl36\")")
            print("or:")
            print(" num_levels(edata,17,36)")
            return 
            

class Uncertainties(Meta):
    __doc__="""Class to handle uncertainty calculations using EGAF data."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    
    def quad_error(self,f,x,dx,y,dy):
        """Internal function: Combines errors in quadrature for calculations 
        of the type:

        f = x*y
        f = x/y

        Arguments:
            f: Result.
            x: First variable.
            dx: First variable uncertainty.
            x: Second variable.
            dx: Second variable uncertainty.
            
        Returns:
            Associated uncertainty.
        """
        self.f = f
        self.x, self.dx = x, dx
        self.y, self.dy = y, dy
        try:
            return float(self.f)*np.sqrt((float(self.dx)/float(self.x))**2+(float(self.dy)/float(self.y))**2)
        except TypeError:
            if self.f == None:
                return 0
            elif self.x == None and self.y == None:
                return 0
            elif self.x == None and float(self.y) > 0.0:
                return float(self.f)*float(self.dy)/float(self.y)
            elif float(self.x) > 0.0 and self.y == None:
                return float(self.f)*float(self.dx)/float(self.x)
        except ZeroDivisionError:
            if self.x != None and float(self.x) > 0:
                return float(self.f)*float(self.dx)/float(self.x)
            elif self.y != None and float(self.y) > 0:
                return float(self.f)*float(dy)/float(self.y)
            else:
                return 0
