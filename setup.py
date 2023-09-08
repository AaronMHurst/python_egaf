import setuptools

requirements = ["numpy", "pandas"]

setuptools.setup(
    name="pyEGAF",
    version="0.1.0",
    url="https://github.com/AaronMHurst/python_egaf",
    author="Aaron M. Hurst",
    author_email="amhurst@berkeley.edu",
    description="Allows for interaction, manipulation, and analysis of thermal-neutron capture gamma-ray data from the EGAF library.",
    long_description=open('README.md').read(),
    license_files=('LICENSE'),
    #packages=setuptools.find_packages(),
    packages=setuptools.find_namespace_packages(),
    #install_requires=[],
    install_requires=requirements,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    include_package_data=True,
    package_data={'': ['EGAF_JSON/*.json', 'EGAF_RIPL/*.dat', 'EGAF_ENSDF/*.ens']},
)
