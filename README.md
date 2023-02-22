# pyEGAF

Python package that allows for the extraction and manipulation of thermal-neutron capture gamma-ray data from the Evaluated Gamma-ray Activation File (EGAF) library [[1]](#1).

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

A `Jupyter` notebook is provided illustrating use of the various methods for access and manipulation of the EGAF data.  Additionally, a few analysis methods commonly adopted in the analysis of thermal-neutron capture data are also included in the `pyEGAF` software package.

# EGAF source data sets

Although the `pyEGAF` methods already provide greatly enhanced user access to the EGAF data, the original data sets are also bundled with this software package.  The data sets are provided in the following three formats:

* ENSDF
* RIPL
* JSON

Each of these formats are described below.





## References
<a id="1">[1]</a>
R.B. Firestone, *et al*.,
*"EGAF: Measurement and Analysis of Gamma-ray Cross Sections"*,
Nucl. Data Sheets **119**, 79 (2014);
https://www.doi.org/10.1016/j.nds.2014.08.024

