import setuptools

requirements = ["numpy", "pandas"]

setuptools.setup(
    name="pyEGAF",
    version="0.1.0",
    url="https://github.com/AaronMHurst/python_egaf",
    author="Aaron M. Hurst",
    author_email="amhurst@berkeley.edu",
    description="Allows for extraction and manipulation of thermal-neutron capture gamma-ray data from the EGAF library.",
    long_description=open('README.md').read(),
    #packages=setuptools.find_packages(),
    packages=setuptools.find_namespace_packages(),
    #install_requires=[],
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    include_package_data=True,
    package_data={'': ['EGAF_JSON/*.json', 'EGAF_RIPL/*.dat', 'EGAF_ENSDF/*.ens']},
)
