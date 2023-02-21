from .base_egaf import *
from .separation import Separation
from .cross_section import CrossSection
from .decay import Levels, Gammas
from .analysis import Analysis

class CapGam(Analysis):
    __doc__="""Class to represent capture-gamma data in CapGam style."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def capgam(self, list, str, *args):
        """Function to display CapGam-style thermal neutron-capture gamma-ray
        data.  The (n,g) reaction, adopted total thermal neutron-capture cross 
        section, keynumber reference, and normalization transition are also 
        output to screen.
        
        Arguments:
            list: A list of EGAF-data JSON objects.
            str: A string object describing the residual nucleus produced in an 
                 (n,gamma) reaction.
            args: The optional string argument object "more".

        Returns: 
            A DataFrame object containing the CapGam-style representation of 
            gamma rays and associated properties.  The columns correspond to:

            [0]: DataFrame ID number;
            [1]: Transition type: primary or secondary gamma (str);
            [2]: Gamma-ray energy (float);
            [3]: Gamma-ray energy uncertainty (float);
            [4]: Relative intensity (float); 
            [5]: Relative intensity uncertainty (float).

            The relative intensities (RI) are normalized to the strongest 
            gamma-ray transition which is defined to have RI=100 (arb. units).

            If the optional string argument "more" is passed to the function,
            a more verbose output will be displayed with additional data 
            captured in the resulting DataFrame object:

            [0]: DataFrame ID number;
            [1]: Transition type: primary or secondary gamma (str);
            [2]: Associated initial level index of gamma transition (int);
            [3]: Associated final level index of gamma transition (int);
            [4]: Associated initial level energy of gamma transition (float);
            [5]: Associated final level energy of gamma transition (float);
            [6]: Gamma-ray energy (float);
            [7]: Gamma-ray energy uncertainty (float);
            [8]: Relative intensity (float); 
            [9]: Relative intensity uncertainty (float).

        Examples:
            For 28Si(n,g)29Si:

            capgam(egaf_data, "Si29")
            capgam(egaf_data, "Si29", "more")
        """
        self.list = list
        self.str = str
        decay_instance = Gammas()
        spe = decay_instance.get_gammas(self.list,self.str,intensity='isotopic')
        try:
            if len(spe) > 0:

                primaries = decay_instance.get_gamma_types(self.list,self.str,intensity='isotopic',gammas='primary')
                secondaries = decay_instance.get_gamma_types(self.list,self.str,intensity='isotopic',gammas='secondary')

                # Find adopted cross section and reference
                adopted_cs, d_adopted_cs = None, None
                unit_adopted_cs, ref_adopted_cs = None, None
                target_nucleus, residual_nucleus = None, None
                for jdict in self.list:
                    if self.str == jdict["nucleusID"]:
                        residual_nucleus = jdict["nucleusID"]
                        target_nucleus = jdict["nucleusTargetID"]
                        for each_n in jdict["neutronCaptureNormalization"]:
                            for each_nr in each_n["normalizationRecord"]:
                                adopted_cs = each_nr["adoptedTotalThermalCaptureCrossSection"]
                                d_adopted_cs = each_nr["dAdoptedTotalThermalCaptureCrossSection"]
                                unit_adopted_cs = each_nr["unitAdoptedCrossSection"]
                                ref_adopted_cs = each_nr["keyNumber"]

                print("Target nucleus: {0}".format(target_nucleus))
                print("Residual (compound nucleus): {0}".format(residual_nucleus))
                print("{0}(n,g){1}".format(target_nucleus,residual_nucleus))
                print("Total radiative thermal neutron-capture cross section = {0} {1} \xb1 {2}".format(adopted_cs, unit_adopted_cs, d_adopted_cs))
                print("Reference: {0}".format(ref_adopted_cs))
                #print("\n")

                # Extract maximum gamma intensity from a DataFrame object

                E = spe[:,4]
                I = spe[:,6]
                dI = spe[:,7]

                df = pd.DataFrame({'E':E, 'I':I, 'dI':dI})
                max_I = df['I'].loc[df['I'].idxmax()]
                E_at_max_I = df['E'].loc[df['I'].idxmax()]

                print("Maximum I = {0} b at E = {1} keV; RI = 100.".format(max_I, E_at_max_I))

                cg_data = []
                for s in spe:
                    level_i = s[0]
                    level_f = s[1]
                    level_energy_i = s[2]
                    level_energy_f = s[3]
                    gamma_energy = s[4]
                    d_gamma_energy =s[5]
                    gamma_intensity = s[6]
                    d_gamma_intensity = s[7]

                    gamma_type = None

                    cap_gam_I = (gamma_intensity/max_I) * 100.0
                    d_cap_gam_I = None
                    #try:
                    #    d_cap_gam_I = cap_gam_I * (d_gamma_intensity/gamma_intensity)
                    #except ZeroDivisionError:
                    #    d_cap_gam_I = 0.0
                    #except RuntimeWarning:
                    #    d_cap_gam_I = 0.0
                    if gamma_intensity > 0.0:
                        d_cap_gam_I = cap_gam_I * (d_gamma_intensity/gamma_intensity)
                    else:
                        d_cap_gam_I = 0.0

                    try:
                        for pri in primaries:
                            if int(level_i) == int(pri[0]) and int(level_f) == int(pri[1]):
                                gamma_type = pri[10]
                    except TypeError:
                        pass

                    try:
                        for sec in secondaries:
                            if int(level_i) == int(sec[0]) and int(level_f) == int(sec[1]):
                                gamma_type = sec[10]
                    except TypeError:
                        pass

                    cg_data.append([int(level_i), int(level_f), level_energy_i, level_energy_f, gamma_energy, d_gamma_energy, cap_gam_I, d_cap_gam_I, gamma_type])


                if len(cg_data) > 0:
                    cg_data = np.array(cg_data)
                    cg_i = cg_data[:,0].astype(int)
                    cg_f = cg_data[:,1].astype(int)
                    cg_Ei = cg_data[:,2].astype(float)
                    cg_Ef = cg_data[:,3].astype(float)
                    cg_E = cg_data[:,4].astype(float)
                    cg_dE = cg_data[:,5].astype(float)
                    cg_RI = cg_data[:,6].astype(float)
                    cg_dRI = cg_data[:,7].astype(float)
                    cg_gtype = cg_data[:,8]

                    df_cg = None
                    if len(args) == 0:
                        df_cg = pd.DataFrame({'Type':cg_gtype, 'E': cg_E, 'dE': cg_dE, 'RI': cg_RI, 'dRI': cg_dRI})
                    elif len(args) == 1:
                        if args[0].lower() == "more":
                            df_cg = pd.DataFrame({'Type':cg_gtype, 'i': cg_i, 'f': cg_f, 'E(i)': cg_Ei, 'E(f)': cg_Ef, 'E': cg_E, 'dE': cg_dE, 'RI': cg_RI, 'dRI': cg_dRI})
                        else:
                            print("Wrong argument!")
                            print("Do you want more output?")
                            print("The only optional argument accepted is \"more\".")
                            return

                    df_cg = df_cg.sort_values('E')

                    return df_cg
        
        except TypeError:
            if spe == None:
                print("No thermal neutron-capture gamma rays for defined input.")
                return
