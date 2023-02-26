# pyEGAF

Python package that allows for the extraction and manipulation of thermal-neutron capture gamma-ray data from the Evaluated Gamma-ray Activation File (EGAF) library [[1]](#1).

![EGAF Nuclides](EGAF_nuclides.png?raw=true "Schematic of nuclear chart relevant to EGAF nuclides")

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

A `Jupyter` notebook is provided illustrating use of the various methods for access and manipulation of the EGAF data.  Additionally, a few analysis methods commonly adopted in the analysis of thermal-neutron capture data are also included in the `pyEGAF` software package.  Launch the notebook provided in the `notebook` folder and execute the cells to run through the example use cases.  This notebook also has a `matplotlib` Python-packagae dependency.

# Docstrings

All `pyEGAF` classes and functions have supporting docstrings.  Please refer to the individual dosctrings for more information on any particular function including how to use it.  The dosctrings for each method generally follow the following template:

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

The JSON data structure is explained in the tables below:

| JSON key | Meaning |
| --- | --- |
|    `"nucleusID"` | A string describing the compound nucleus <symbol><mass> |
|    `"datasetType"` | A string to identify the data set |
|    `"nucleusZ"` | A number denoting the atomic number of the compound nucleus |
|    `"nucleusA"` | A number denoting the mass number of the compound nucleus |
|    `"nucleusN"` | A number denoting the neutron number of the compound nucleus |
|    `"nucleusTargetZ"` | A number denoting the atomic number of the target nucleus |
|    `"nucleusTargetA"` | A number denoting the mass number of the target nucleus |
|    `"nucleusTargetN"` | A number denoting the neutron number of the target nucleus |
|    `"nucleusTargetElement"` | A one- or two-character string denoting the chemical element ID |
|    `"nucleusTargetID"` | A string describing the target nucleus <symbol><mass> |
|    `"numberPrimaryGammas"` | A number denoting the number of primary &gamma; rays |
|    `"numberSecondaryGammas"` | A number denoting the number of secondary &gamma; rays |
|    `"totalNumberLevels"` | A number denoting the total number of levels in the decay scheme |
|    `"totalNumberGammas"` | A number denoting the total number of &gamma; rays in the decay scheme |
|    `"unitEnergy"` | A string to indicate the units of the energy quantities |
| `"recordQ"` | An array containing information from the *Q*-value record for the compound nucleus
| `"neutronCaptureNormalization"` | An array containing normalization information for the compound nucleus |
| `"levelScheme"` | An array contating decay-scheme information for the compound nucleus |

The JSON arrays are desbcribed below:

### `"recordQ"` array

| JSON key | Meaning |
| --- | --- |
|           `"energyNeutronSeparationAME2020"` | A number denoting the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus |
|           `"dEnergyNeutronSeparationAME2020"` | A number denoting the uncertainty for the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus |
|            `"energyProtonSeparationAME2020"` | A number denoting the AME2020 [[5]](#5) proton-separation energy of the compound nucleus |
|            `"dEnergyProtonSeparationAME2020"` | A number denoting the uncertainty for the AME2020 [[5]](#5) proton-separation energy of the compound nucleus |
|            `"energyNeutronSeparationENSDF"` | A number denoting the ENSDF neutron-separation energy of the compound nucleus |
|            `"energyProtonSeparationENSDF"` | A number denoting the ENSDF proton-separation energy of the compound nucleus |
|            `"energyNeutronSeparationEGAF"` | A number denoting the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus |
|            `"dEnergyNeutronSeparationEGAF"` | A number denoting the uncertainty for the AME2020 [[5]](#5) neutron-separation energy of the compound nucleus |

### `"neutronCaptureNormalization"` array

This array contains the `"normalizationRecord"` JSON object, an array with the following contents:

| JSON key | Meaning |
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

| JSON key | Meaning |
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

| JSON key | Meaning |
| --- | --- |
| `"halfLifeBest"` | A number type representing the halflife in *best* units from original data set.|
| `"dHalfLifeBest"` | A number type representing the associated uncertainty on the halflife in *best* units. |
| `"unitHalfLifeBest"` | A string type to indicate the *best* halflife units.|
| `"halfLifeConverted"` | A number type representing the halflife *converted* to units of seconds.|
| `"dHalfLifeConverted"` | A number type representing the associated uncertainty on the halflife *converted* to seconds.|
| `"unitHalfLifeConverted"` | A string type to indicate the *converted* halflife units.|


### `"spins"` array

| JSON key | Meaning |
| --- | --- |
|  `"spinIndex"`| A number type (integer) associated with the indexed sequence of the spin-parity perumations.|
|  `"spinReal"`| A number type (float) corresponding to the real spin value of the level. |
|  `"spinIsTentative"` | A boolean type to flag tentative spin assignments. |
|  `"spinIsLimit"` | A boolean type to flag levels with spin values expressed as limits. |
|  `"spinLimits"` | A string type representing the associated spin limits of the level; a `null` value is given if the level does not have any spin limits. |
|  `"parity"` | A number type (integer) that represents the parity of the level: -1 (negative &pi;), 1 (positive &pi;), 0 (no &pi; assignment). |
|  `"paritySign"` | A string type referring to the parity (*"negative"*, ""positive"*, or *null*) of the level. |
|  `"parityIsTentative"` | A boolean type to flag tentative parity assignments. |

### `"gammaDecay"` array

| JSON key | Meaning |
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
|  `"calculatedTotalInternalConversionCoefficient"`| A number type corresponding to the BrIcc-calculated [[6]](#6) total internal-conversion coefficient associated with the &gamma;-ray transition.|
|  `"dCalculatedTotalInternalConversionCoefficient"`| A number type corresponding to the associated uncertainty of the BrIcc-calculated [[6]](#6) total internal-conversion coefficient.|
|  `"calculatedAtomicShellConversionCoefficients"` | An array type containing the atomic subshell internal-conversion coefficient data.|

### `"gammaAbsoluteIntensities"` array

| JSON key | Meaning |
| --- | --- |
| "partialElementalCrossSection" | A number type corresponding to the *elemental* partial &gamma;-ray cross section.|
|  "dPartialElementalCrossSection" | A number type corresponding to the associated uncertainty of the *elemental* partial &gamma;-ray cross section.|
|  "partialIsotopicCrossSection"| A number type corresponding to the partial &gamma;-ray cross section corrected for *isotopic* abundance.|
|  "dPartialIsotopicCrossSection"| A number type corresponding to the associated uncertainty of the *isotopically-corrected* partial &gamma;-ray cross section.|
|  "populationPerNeutronCapture"| A number type corresponding to the population per neutron capture for the &gamma;-ray transition.|
|  "dPopulationPerNeutronCapture"| A number type corresponding to the associated uncertainty of the population per neutron capture for the &gamma;-ray transition.|

### `"calculatedAtomicShellConversionCoefficients"` array

| JSON key | Meaning |
| --- | --- |
|  `"calculatedInternalConversionCoefficientAtomicShellK"`|
|  `"dCalculatedInternalConversionCoefficientAtomicShellK"`|
|  `"calculatedInternalConversionCoefficientAtomicShellL"`|
|  `"dCalculatedInternalConversionCoefficientAtomicShellL"`|
|  `"calculatedInternalConversionCoefficientAtomicShellM"`|
|  `"dCalculatedInternalConversionCoefficientAtomicShellM"`|
|  `"calculatedInternalConversionCoefficientAtomicShellN"`|
|  `"dCalculatedInternalConversionCoefficientAtomicShellN"`|
|  `"calculatedInternalConversionCoefficientAtomicShellO"`|
|  `"dCalculatedInternalConversionCoefficientAtomicShellO"`|
|  `"calculatedInternalConversionCoefficientAtomicShellP"`|
|  `"dCalculatedInternalConversionCoefficientAtomicShellP"`|
|  `"calculatedInternalConversionCoefficientAtomicShellQ"`|
|  `"dCalculatedInternalConversionCoefficientAtomicShellQ"`|

## References
<a id="1">[1]</a>
R.B. Firestone *et al*.,
*"EGAF: Measurement and Analysis of Gamma-ray Cross Sections"*,
Nucl. Data Sheets **119**, 79 (2014);
https://www.doi.org/10.1016/j.nds.2014.08.024

<a id="2">[2]</a>
J.K. Tuli,
*"Evaluated Nuclear Structure Data File"*,
BNL-NCS-51655-01/02-Rev (2001).

<a id="3">[3]</a>
Evaluated Gamma-ray Activation File (EGAF);
https://www-nds.iaea.org/pgaa/egaf.html

<a id="4">[4]</a>
R. Capote *et al*.,
*"Reference Input Parameter Library"*,
Nucl. Data Sheets **110**, 3107 (2009).

<a id="5">[5]</a>
M. Wang, W.J. Huang, F.G. Kondev, G. Audi, S. Naimi,
*The AME2020 atomic mass evaluation*,
Chin. Phys. C **45**, 030003 (2021).
