from .base_egaf import *
from .separation import Separation
from .cross_section import CrossSection
from .decay import Levels, Gammas
from .cap_gam import CapGam

from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

class RIPL(CapGam):
    __doc__="""Class to handle RIPL-formatted EGAF data sets."""
    
    def __init__(self):
        BaseEGAF.__init__(self)
        Meta.__init__(self)
        Uncertainties.__init__(self)
        Separation.__init__(self)
        CrossSection.__init__(self)
        Levels.__init__(self)
        Gammas.__init__(self)
        CapGam.__init__(self)
        
    def get_ripl(self,list,str,bool=False):
        """Display RIPL-formatted EGAF data in the console.  The RIPL-formatted 
        EGAF data set may also be written to file.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            str: The ID of the residual compound nucleus passed as a string
                 argument.
            bool: True: The corresponding RIPL file will be printed in the 
                        current working directory.
                  False: The RIPL will not be printed (default argument).

        Returns: 
            A list object corresponding to the appropriate RIPL-formatted EGAF 
            data set; an ASCII text dump of the RIPL-formatted data is also 
            displayed in the console.  The RIPL-formatted EGAF object may also 
            be written to file according to the value of the boolean argument 
            passed to the function.

        Example:
            To display the RIPL data for 28Si(n,g)29Si:
            get_ripl(edata, "Si29")

            To display the RIPL data for 28Si(n,g)29Si and write "EGAF_RIPL_Si29.dat"
            to file :
            get_ripl(edata, "Si29", True)
        """
        # Load the RIPL files to a list:
        from . import get_data
        EGAF_RIPL_PATH = get_data('EGAF_RIPL')
        ripl_egaf_list = [r for r in glob.glob("%s/*.dat"%EGAF_RIPL_PATH)]
        RIPL_COUNT = len(ripl_egaf_list)
        
        self.list = list
        self.str = str
        self.bool = bool

        Z_res = None
        A_res = None
        res_ID = None
        targ_ID = None
        for jdict in self.list:
            if self.str == jdict["nucleusID"]:
                Z_res = jdict["nucleusZ"]
                A_res = jdict["nucleusA"]
                res_ID = jdict["nucleusID"]
                targ_ID = jdict["nucleusTargetID"]

        if Z_res != None and A_res != None:

            ripl_list = []
            for ripl_file in ripl_egaf_list:
                with open(ripl_file, 'r') as rf:
                    RIPL_MATCH = False
                    for i,line in enumerate(rf):
                        if i == 0:
                            line_cols = line.split()
                            if int(A_res) == int(line_cols[1]) and int(Z_res) == int(line_cols[2]):
                                RIPL_MATCH = True
                        if RIPL_MATCH == True:
                            print(line.strip('\n'))
                            ripl_list.append(line)
                    rf.close()

            if self.bool == True:
                if len(ripl_list) > 0:
                    ripl_file = "EGAF_RIPL_{0}_NG_{1}.dat".format(targ_ID, res_ID)
                    ripl_out = open(ripl_file, mode='w')
                    for line in ripl_list:
                        ripl_out.write("{0}".format(line))
                    print("{0} written to current working directory.".format(ripl_file))
                    ripl_out.close()

            with suppress_stdout():
                return ripl_list

        else:
            print("No match found for RIPL-formatted EGAF data set.")
            return
        

class JSONFile(RIPL):
    __doc__="""Class to handle JSON-structured EGAF data sets."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def get_json(self,list,str,bool=False):
        """Display JSON-formatted EGAF data in the console.  The JSON-formatted 
        EGAF data set may also be written to file.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            str: The ID of the residual compound nucleus passed as a string
                 argument.
            bool: True: The corresponding JSON file will be printed in the 
                        current working directory.
                  False: The JSON will not be printed (default argument).

        Returns: 
            A dictionary object containing to the appropriate JSON-structured 
            EGAF data.  The JSON-formatted EGAF object may also be written to 
            file depending on the value of the boolean argument passed to the 
            function.

        Example:
            To get the JSON data for 28Si(n,g)29Si:
            get_json(edata, "Si29")

            To get the JSON data for 28Si(n,g)29Si and write 
            "EGAF_JSON_Si28_NG_Si29.json" to file :
            get_json(edata, "Si29", True)
        """
        # Load the JSON files to a list:
        from . import get_data
        EGAF_JSON_PATH = get_data('EGAF_JSON')
        json_egaf_list = [j for j in glob.glob("%s/*.json"%EGAF_JSON_PATH)]
        JSON_COUNT = len(json_egaf_list)

        self.list = list
        self.str = str
        self.bool = bool

        res_ID = None
        targ_ID = None
        for jdict in self.list:
            if self.str == jdict["nucleusID"]:
                res_ID = jdict["nucleusID"]
                targ_ID = jdict["nucleusTargetID"]

                if self.bool == True:
                    with open("EGAF_JSON_{0}_NG_{1}.json".format(targ_ID,res_ID), "w") as jf:
                        json.dump(jdict, jf, indent=4, ensure_ascii=False)
                        jf.close()
                        print("{0} written to current working directory.".format(jf.name))

                return jdict
            
        if res_ID == None and targ_ID == None:
            print("No match found for JSON-formatted EGAF data set.")
            return
        

class ENSDF(JSONFile):
    __doc__="""Class to handle JSON-structured EGAF data sets."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def get_ensdf(self,list,str,bool=False):
        """Display ENSDF-formatted EGAF data in the console.  The 
        ENSDF-formatted EGAF data set may also be written to file.

        Arguments:
            list: A list of EGAF (n,g) data JSON objects.
            str: The ID of the residual compound nucleus passed as a string
                 argument.
            bool: True: The corresponding ENSDF file will be printed in the 
                        current working directory.
                  False: The ENSDF will not be printed (default argument).

        Returns: 
            A list object containing to the appropriate ENSDF-formatted EGAF 
            data set.  The ENSDF-formatted EGAF object may also be written to 
            file depending on the value of the boolean argument passed to the 
            function.

        Example:
            To get the ENSDF-formatted EGAF data for 28Si(n,g)29Si:
            get_ensdf(edata, "Si29")

            To get the ENSDF-formatted EGAF data for 28Si(n,g)29Si and write 
            "EGAF_ENSDF_28SI_NG_29SI.ens" to file :
            get_ensdf(edata, "Si29", True)
        """
        # Load the ENSDF files to a list:
        from . import get_data
        EGAF_ENSDF_PATH = get_data('EGAF_ENSDF')
        ensdf_egaf_list = [j for j in glob.glob("%s/*.ens"%EGAF_ENSDF_PATH)]
        ENSDF_COUNT = len(ensdf_egaf_list)

        self.list = list
        self.str = str
        self.bool = bool

        chem_symb = re.search(r'\D+', self.str).group(0)
        res_A = re.findall(re.compile(r'\d+'), self.str)[0]
        targ_A = "{0}".format(int(res_A) - 1)

        ds_pattern = targ_A+chem_symb.upper() + '_ng_' + res_A+chem_symb.upper()

        filename = None
        ensdf_data = None
        ENSDF_MATCH = False
        for egaf_file in ensdf_egaf_list:
            if ds_pattern in egaf_file:
                ENSDF_MATCH = True
                filename = egaf_file
                
                with open(filename, mode='r') as rf:
                    ensdf_data = rf.readlines()
                    [print(line.strip('\n')) for line in ensdf_data]
                    rf.close()
                    
                    if self.bool == True:
                        with open("EGAF_ENSDF_{0}.ens".format(ds_pattern.upper()), mode="w") as wf:
                            for line in ensdf_data:
                                wf.write("{0}".format(line))
                            print("{0} written to current working directory.".format(wf.name))
                            wf.close()

        if ENSDF_MATCH == False:
            print("No match found for ENSDF-formatted EGAF data set.")
            return
        else:
            return ensdf_data
        
