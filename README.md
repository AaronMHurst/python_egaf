# pyEGAF

This project is a Python package that allows for the extraction and manipulation of thermal-neutron capture gamma-ray data from the Evaluated Gamma-ray Activation File (EGAF) library [[1]](#1).  The EGAF library is a database of &gamma;-ray energies and their corresponding partial &gamma;-ray cross sections from thermal-neutron capture measurements carried out with a guided neutron beam at the Budapest Research Reactor for 245 isotopes encompassing measurements of natural elemental samples for targets from *Z* = 1-83, 90, and 92, except for Tc (*Z* = 43) and Pm (*Z* = 61).  The database comprises a total of 8172 primary &gamma; rays and 29605 secondary &gamma; rays (a total of 37777 &gamma; rays) associated with 12564 levels.  The (*n*,&gamma;) targets and corresponding residual compound nuclides relevant to the EGAF project are summarized in the schematic of the nuclear chart shown in the figure below.

![EGAF Nuclides](EGAF_nuclides.png?raw=true "Schematic of nuclear chart relevant to EGAF nuclides")

The `pyEGAF` package provides users with a convenient means for access and analysis of the thermal neutron-capture data in EGAF including decay schemes and associated properties for all compound nuclides, as well as a capability to search for by &gamma;-ray energy for forensics applications.  After building, installation, and testing of the project, refer to the examples provided in the `Jupyter` notebook for an overview regarding some of the available methods.

# Building and installation

This project should be built and installed by executing the `installation.sh` script at the terminal command line of the project directory:

```Bash
$ git clone https://github.com/AaronMHurst/python_egaf.git
$ cd python_egaf
$ sh installation.sh
```

# Testing

A suite of Python modules containing 208 unit tests have been written for this project and are located in the `tests` folder.  To run the test suite and ensure they work with the Python environment, run `tox` in the project directory where the `tox.ini` file is also located:

```Bash
$ tox -r
```

This project has the following Python-package dependencies: `numpy`, `pandas`, and `pytest`.  The test session is automatically started after building against the required Python environment.

# Running the software

After the installation the `pyEGAF` scripts can be ran from any location by importing the package and making an instance of the `EGAF` class:

```Bash
$ python
>>> import pyEGAF as egaf
>>> e = egaf.EGAF()
```

Most methods also require passing the EGAF `JSON` source data set as a list object:

```
>>> edata = e.load_egaf()
```

A `Jupyter` notebook is provided illustrating use of the various methods for access and manipulation of the EGAF data.  Additionally, a few analysis methods commonly adopted in the analysis of thermal-neutron capture data are also included in the `pyEGAF` software package.  Launch the notebook provided in the `notebook` folder and execute the cells to run through the example-use cases.  This notebook also has a `matplotlib` Python-package dependency.

# Docstrings

All `pyEGAF` classes and functions have supporting docstrings.  Please refer to the individual dosctrings for more information on any particular function including how to use it.  The dosctrings for each method generally have the following structure:

* A short explanation of the function.
* A list and description of arguments that need to be passed to the function.
* The return value of the function.
* An example(s) invoking use of the function.

# EGAF source data sets

Although the `pyEGAF` methods already provide greatly enhanced user access to the EGAF data, the original data sets are also bundled with this software package for convenience and to allow users to curate data in a bespoke manner should they prefer.  The data sets are provided in the following three formats:

* Evaluated Nuclear Structure Data File (ENSDF);
* Reference Input Parameter Library (RIPL);
* JavaScript Object Notation (JSON).

Each of these formats are described below.

## ENSDF format

The original EGAF data sets were prepared in accordance with the mixed-record 80-character column format of the Evaluated Nuclear Structure Data File (ENSDF) [[2]](#2).  These ENSDF-formatted files are maintained online by the International Atomic Energy Agency [[3]](#3).  The relevant fields of the `Normalization`, `Level`, and `Gamma` records that are commonly adopted in the EGAF data sets are explained in the manual [[2]](#2).  In addition, `Comment` records are also frequently encountered in EGAF data sets.  The ENSDF-formatted EGAF data sets can be accessed from the project folder by changing into the following directory and listing its contents:

```Bash
$ cd pyEGAF/EGAF_ENSDF
$ ls
```

Alternatively, individual files can also be accessed using `pyEGAF` methods by passing the EGAF data set list object and the *residual compound nucleus* produced in an (*n*,&gamma;), for example, <sup>28</sup>Si(*n*,&gamma;)<sup>29</sup>Si:

```Bash
>>> ensdf = e.get_ensdf(edata, "Si29")
```

## RIPL format

Because many nuclear reaction codes source decay-scheme information in a particular Reference Input Parameter Library (RIPL) [[4]](#4) format, representative RIPL-translated data sets have also been generated for each corresponding EGAF data set and are bundled with the software.  The RIPL-formatted EGAF data sets are located in the `python_egaf/pyEGAF/EGAF_RIPL` directory.  These files can also be accessed from the interpreter, for example, <sup>28</sup>Si(*n*,&gamma;)<sup>29</sup>Si:

```Bash
>>> ripl = e.get_ripl(edata, "Si29")
```

The proton- and neutron-separation energies in the RIPL headers are taken from the 2020 Atomic Mass Evaluation [[5]](#5).

## JSON format

All original EGAF data sets have been translated into a representative JavaScript Object Notation (JSON) format using an intuitive syntax to describe the quantities sourced from the primary and continuation records of the ENSDF-formatted data sets.  The JSON-formatted data sets are also bundled with the software package and are located in `python_egaf/pyEGAF/EGAF_JSON`.  Again, individual data sets can be accessed through the interpreter, for example, <sup>28</sup>Si(*n*,&gamma;)<sup>29</sup>Si:

```Bash
>>> jfile = e.get_json(edata, "Si29")
```

The JSON data structures support the following data types:

* *string*
* *number*
* *boolean*
* *null*
* *object* (JSON object)
* *array*

The JSON-formatted EGAF schema is explained in the tables below:

| JSON key | Explanation |
| --- | --- |
|    `"nucleusID"` | A string type describing the compound nucleus `<symbol><mass>`.|
|    `"datasetType"` | A string type to identify the data set.|
|    `"nucleusZ"` | A number type denoting the atomic number of the compound nucleus.|
|    `"nucleusA"` | A number type denoting the mass number of the compound nucleus.|
|    `"nucleusN"` | A number type denoting the neutron number of the compound nucleus.|
|    `"nucleusTargetZ"` | A number type denoting the atomic number of the target nucleus.|
|    `"nucleusTargetA"` | A number type denoting the mass number of the target nucleus.|
|    `"nucleusTargetN"` | A number type denoting the neutron number of the target nucleus.|
|    `"nucleusTargetElement"` | A string type (one- or two-character) denoting the chemical element ID.|
|    `"nucleusTargetID"` | A string type describing the target nucleus `<symbol><mass>`.|
|    `"numberPrimaryGammas"` | A number type denoting the number of primary &gamma; rays.|
|    `"numberSecondaryGammas"` | A number type denoting the number of secondary &gamma; rays.|
|    `"totalNumberLevels"` | A number type denoting the total number of levels in the decay scheme.|
|    `"totalNumberGammas"` | A number type denoting the total number of &gamma; rays in the decay scheme.|
|    `"unitEnergy"` | A string type to indicate the units of the energy quantities.|
| `"recordQ"` | An array type containing information from the *Q*-value record for the compound nucleus.|
| `"neutronCaptureNormalization"` | An array type containing normalization information for the compound nucleus.|
| `"levelScheme"` | An array type containing decay-scheme information for the compound nucleus.|

The JSON arrays are described below:

### `"recordQ"` array

| JSON key | Explanation |
| --- | --- |
|           `"energyNeutronSeparationAME2020"` | A number type denoting the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus.|
|           `"dEnergyNeutronSeparationAME2020"` | A number type denoting the uncertainty for the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus.|
|            `"energyProtonSeparationAME2020"` | A number type denoting the AME2020 [[5]](#5) proton-separation energy of the compound nucleus.|
|            `"dEnergyProtonSeparationAME2020"` | A number type denoting the uncertainty for the AME2020 [[5]](#5) proton-separation energy of the compound nucleus.|
|            `"energyNeutronSeparationENSDF"` | A number type denoting the ENSDF neutron-separation energy of the compound nucleus.|
|            `"energyProtonSeparationENSDF"` | A number type denoting the ENSDF proton-separation energy of the compound nucleus.|
|            `"energyNeutronSeparationEGAF"` | A number type denoting the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus.|
|            `"dEnergyNeutronSeparationEGAF"` | A number type denoting the uncertainty for the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus.|

### `"neutronCaptureNormalization"` array

This array contains the `"normalizationRecord"` JSON object, an array with the following contents:

| JSON key | Explanation |
| --- | --- |
| `"multiplierIsotopicCorrection"` | A number type corresponding to elemental-isotopic conversion factor. |
|   `"dMultiplierIsotopicCorrection"` | A number type corresponding to the uncertainty for the elemental-isotopic conversion factor. |
|   `"naturalIsotopicAbundance"` | A number type corresponding to the isotopic abundance of the EGAF target.|
|   `"dNaturalIsotopicAbundance"` | A number type corresponding to the uncertainty for the isotopic abundance of the EGAF target.|
|   `"adoptedTotalThermalCaptureCrossSection"` | A number type corresponding to the adopted total thermal neutron capture cross section. |
|   `"dAdoptedTotalThermalCaptureCrossSection"` | A number type corresponding to the uncertainty for the adopted total thermal neutron capture cross section. |
|   `"unitAdoptedCrossSection"` | A string type to indicate the units of the cross-section quantities. |
|   `"keyNumber"` | A string type corresponding to the keynumber reference of the adopted cross section. |

### `"levelScheme"` array

| JSON key | Explanation |
| --- | --- |
|            `"levelIndex"` | A number type (integer) corresponding to unique index associated with an energy level. |
|            `"levelEnergy"` | A number type (float) corresponding to the level excitation energy. |
|            `"dLevelEnergy"` | A number type (float) corresponding to the uncertainty of the level excitation energy. |
|            `"levelIsIsomer"` | A boolean type to flag levels with isomeric properties. |
|            `"isomerDecay"` | An array type corresponding to the isomer-decay properties of the level. |
|           `"numberOfSpins"` | A number type (integer) corresponding to the number of spin-parity permutations of the level. |
|	    `"spins"` | An array type corresponding to the spin-parity information associated with the level. |
|           `"numberOfGammas"` | A number type (integer) corresponding to the number of deexcitation &gamma; rays belonging to the levels. |
|	    `"gammaDecay"` | An array type corresponding to the &gamma;-decay properties of the level. |

### `"isomerDecay"` array

| JSON key | Explanation |
| --- | --- |
| `"halfLifeBest"` | A number type representing the halflife in *best* units from original data set.|
| `"dHalfLifeBest"` | A number type representing the associated uncertainty on the halflife in *best* units. |
| `"unitHalfLifeBest"` | A string type to indicate the *best* halflife units.|
| `"halfLifeConverted"` | A number type representing the halflife *converted* to units of seconds.|
| `"dHalfLifeConverted"` | A number type representing the associated uncertainty on the halflife *converted* to seconds.|
| `"unitHalfLifeConverted"` | A string type to indicate the *converted* halflife units.|


### `"spins"` array

| JSON key | Explanation |
| --- | --- |
|  `"spinIndex"`| A number type (integer) associated with the indexed sequence of the spin-parity permutations.|
|  `"spinReal"`| A number type (float) corresponding to the real spin value of the level. |
|  `"spinIsTentative"` | A boolean type to flag tentative spin assignments. |
|  `"spinIsLimit"` | A boolean type to flag levels with spin values expressed as limits. |
|  `"spinLimits"` | A string type representing the associated spin limits of the level; a `null` value is given if the level does not have any spin limits. |
|  `"parity"` | A number type (integer) that represents the parity of the level: -1 (negative &pi;), 1 (positive &pi;), 0 (no &pi; assignment). |
|  `"paritySign"` | A string type referring to the parity (*"negative"*, *"positive"*, or *null*) of the level. |
|  `"parityIsTentative"` | A boolean type to flag tentative parity assignments. |

### `"gammaDecay"` array

In this version of the EGAF database internal conversion coefficients have not, in general, been included in the source data sets and have consequently been set to zero for both total and individual-shell contributions.  A future release of the database will include BrIcc-calculated [[6]](#6) values for both total conversion as well as contributions from individual shells.

| JSON key | Explanation |
| --- | --- |
|  `"gammaEnergy"` | A number type corresponding to the &gamma;-ray energy.|
|  `"dGammaEnergy"` | A number type corresponding to the associated uncertainty of the &gamma;-ray energy.|
|  `"levelIndexInitial"` | A number type (integer) corresponding to the index of the initial level associated with the &gamma;-ray transition.|
|  `"levelIndexFinal"` | A number type (integer) corresponding to the index of the final level associated with the &gamma;-ray transition.|
|  `"levelEnergyInitial"` | A number type (float) corresponding to the excitation energy of the initial level associated with the &gamma;-ray transition.|
|  `"levelEnergyFinal"` | A number type (float) corresponding to the excitation energy of the final level associated with the &gamma;-ray transition.|
|  `"gammaTransitionType"` | A string type denoting *"primary"* or *"secondary"* &gamma;-ray transition types.|
|  `"gammaFeedsGroundState"` | A boolean type to flag &gamma;-ray transitions that feed the ground state.|
|  `"gammaAbsoluteIntensities"` | An array type containing &gamma;-ray intensity information associated with the transition.|
|  `"multipolarity"` | A string type (or *null*) describing the multipolarity (*"M1", "E1", "M2", "E2", "M1+E2",* etc.) of the &gamma;-ray transition.|
|  `"multipolarityIsTentative"`| A boolean type to flag tentative multipolarity assignments.|
|  `"multipolarityIsAssumed"`| A boolean type to indicate evaluator-assumed multipolarity assignments.|
|  `"mixingRatio"`| A number type (or *null*) corresponding to the &gamma;-ray mixing ratio where known.|
|  `"dMixingRatio"`| A number type (or *null*) corresponding to the associated uncertainty of the &gamma;-ray mixing ratio.|
|  `"mixingRatioSign"`| A string type (or *null*) corresponding to the sign (*"positive"* or *"negative"*) of the &gamma;-ray mixing ratio.
|  `"calculatedTotalInternalConversionCoefficient"`| A number type corresponding to the calculated total internal-conversion coefficient associated with the &gamma;-ray transition.|
|  `"dCalculatedTotalInternalConversionCoefficient"`| A number type corresponding to the associated uncertainty of the calculated total internal-conversion coefficient.|
|  `"calculatedAtomicShellConversionCoefficients"` | An array type containing the atomic subshell internal-conversion coefficient data.|

### `"gammaAbsoluteIntensities"` array

The source EGAF data sets contain partial elemental &gamma;-ray cross sections &sigma;<sub>&gamma;<sub>*e*</sub></sub> (JSON key: `"partialElementalCrossSection"`) in the `RI` field of the `Gamma` record [[2]](#2).  The isotopically-corrected partial &gamma;-ray cross sections &sigma;<sub>&gamma;<sub>*i*</sub></sub> (JSON key: `"partialIsotopicCrossSection"`) are derived according to

$$ \sigma_{\gamma_{i}} = \sigma_{\gamma_{e}} M, \quad (1) $$

where *M* is the photon intensity multiplier from the `NR` field of the `Normalization` record [[2]](#2) (JSON key: `"multiplierIsotopicCorrection"` from the `"normalizationRecord"` JSON object).  The &gamma;-ray populations per neutron capture *P*<sub>&gamma;</sub> are then deduced from

$$ P_{\gamma} = \frac{\sigma_{\gamma_{i}}}{\sigma_{0}}. \quad (2) $$

where &sigma;<sub>0</sub> is the adopted total thermal neutron capture cross section (JSON key: `"adoptedTotalThermalCaptureCrossSection"` from the `"normalizationRecord"` JSON object).

| JSON key | Explanation |
| --- | --- |
| `"partialElementalCrossSection"` | A number type corresponding to the *elemental* partial &gamma;-ray cross section.|
|  `"dPartialElementalCrossSection"` | A number type corresponding to the associated uncertainty of the *elemental* partial &gamma;-ray cross section.|
|  `"partialIsotopicCrossSection"`| A number type corresponding to the partial &gamma;-ray cross section corrected for *isotopic* abundance.|
|  `"dPartialIsotopicCrossSection"`| A number type corresponding to the associated uncertainty of the *isotopically-corrected* partial &gamma;-ray cross section.|
|  `"populationPerNeutronCapture"`| A number type corresponding to the population per neutron capture for the &gamma;-ray transition.|
|  `"dPopulationPerNeutronCapture"`| A number type corresponding to the associated uncertainty of the population per neutron capture for the &gamma;-ray transition.|

### `"calculatedAtomicShellConversionCoefficients"` array

| JSON key | Explanation |
| --- | --- |
|  `"calculatedInternalConversionCoefficientAtomicShellK"`| A number type corresponding to the calculated *K*-shell internal-conversion coefficient.|
|  `"dCalculatedInternalConversionCoefficientAtomicShellK"`| A number type corresponding to the associated uncertainty of the calculated *K*-shell internal-conversion coefficient.|
|  `"calculatedInternalConversionCoefficientAtomicShellL"`| A number type corresponding to the calculated *L*-shell internal-conversion coefficient.|
|  `"dCalculatedInternalConversionCoefficientAtomicShellL"`| A number type corresponding to the associated uncertainty of the calculated *L*-shell internal-conversion coefficient.|
|  `"calculatedInternalConversionCoefficientAtomicShellM"`| A number type corresponding to the calculated *M*-shell internal-conversion coefficient.|
|  `"dCalculatedInternalConversionCoefficientAtomicShellM"`| A number type corresponding to the associated uncertainty of the calculated *M*-shell internal-conversion coefficient.|
|  `"calculatedInternalConversionCoefficientAtomicShellN"`| A number type corresponding to the calculated *N*-shell internal-conversion coefficient.|
|  `"dCalculatedInternalConversionCoefficientAtomicShellN"`| A number type corresponding to the associated uncertainty of the calculated *N*-shell internal-conversion coefficient.|
|  `"calculatedInternalConversionCoefficientAtomicShellO"`| A number type corresponding to the calculated *O*-shell internal-conversion coefficient.|
|  `"dCalculatedInternalConversionCoefficientAtomicShellO"`| A number type corresponding to the associated uncertainty of the calculated *O*-shell internal-conversion coefficient.|
|  `"calculatedInternalConversionCoefficientAtomicShellP"`| A number type corresponding to the calculated *P*-shell internal-conversion coefficient.|
|  `"dCalculatedInternalConversionCoefficientAtomicShellP"`| A number type corresponding to the associated uncertainty of the calculated *P*-shell internal-conversion coefficient.|
|  `"calculatedInternalConversionCoefficientAtomicShellQ"`| A number type corresponding to the calculated *Q*-shell internal-conversion coefficient.|
|  `"dCalculatedInternalConversionCoefficientAtomicShellQ"`| A number type corresponding to the associated uncertainty of the calculated *Q*-shell internal-conversion coefficient.|

## References
<a id="1">[1]</a>
R.B. Firestone *et al*.,
*"Database of Prompt Gamma Rays from Slow Thermal Neutron Capture for Elemental Analysis"*,
IAEA STI/PUB/1263, 251 (2007).

<a id="2">[2]</a>
J.K. Tuli,
*"Evaluated Nuclear Structure Data File"*,
BNL-NCS-51655-01/02-Rev (2001).

<a id="3">[3]</a>
Evaluated Gamma-ray Activation File (EGAF);
https://www-nds.iaea.org/pgaa/egaf.html

<a id="4">[4]</a>
R. Capote *et al*.,
*"RIPL - Reference Input Parameter Library for Calculation of Nuclear Reactions and Nuclear Data Evaluations"*,
Nucl. Data Sheets **110**, 3107 (2009).

<a id="5">[5]</a>
M. Wang, W.J. Huang, F.G. Kondev, G. Audi, S. Naimi,
*"The AME2020 atomic mass evaluation"*,
Chin. Phys. C **45**, 030003 (2021).

<a id="6">[6]</a>
T. Kibedi, T.W. Burrows, M.B. Trzhaskovskaya, P.M. Davidson, C.W. Nestor Jr.,
*"Evaluation of theoretical conversion coefficients using BrIcc"*,
Nucl. Intrum. Methods Phys. Res. Sect. A **589**, 202 (2008).
