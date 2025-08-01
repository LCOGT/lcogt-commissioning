import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lcocommissioning",
    version="2.1.13",
    author="Daniel Harbeck",
    author_email="dharbeck@lco.global",
    description="Tool to characterize CCD detectors and other commissioning tasks for the LCO observatory.",
    long_description=long_description,
    url="https://github.com/LCOGT/lcogt-commissioning",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points = {
        'console_scripts': ['submitNRESObservation = lcocommissioning.nres.submitNRESObservation:main',
                            'plot_nres_spec = lcocommissioning.nres.plot_nres_1d:main',
                            'noisegainmef = lcocommissioning.noisegainrawmef:main',
                            'submitXtalkObservation = lcocommissioning.submitXtalkObservation:main',
                            'analysegainhistory = noisegaincrawler.analysegainhistory:main',
                            'crawlnoisegain = noisegaincrawler.crawl_noisegain:main',
                            'submitNamedModeTest = lcocommissioning.submitNameModeTest:main',
                            'sinistrocrosstalk = lcocommissioning.crosstalk:main',
                            'focuscurve = lcocommissioning.focus.focuscurve:main',
                            'effocus = lcocommissioning.focus.ef_focuscalibration:main',
                            'focusstack = lcocommissioning.focus.focustrends:main',
                            'submit_floyds_calibration = lcocommissioning.floyds.submitFloydsCalibration:main',
                            'submit_floyds_observation = lcocommissioning.floyds.submitFloydsObservation:main',
                            'submit_muscat_observation = lcocommissioning.muscat.submitMuscatObservation:main',
                            'submit_cdk_observation = lcocommissioning.cdk.submit_cdk_observation:main',
                            'get_sbigframe = lcocommissioning.sbig.get_sbig_frame:main',
                            'get_qhyframe = lcocommissioning.cmostest.getqhyccdframe:main',
                            'get_archonframe = lcocommissioning.archon.archonexpose:main',
                            'elongation_assessment = lcocommissioning.cdk.elongation_assessment:main',
                            'muscat_focus =  lcocommissioning.muscat.muscat_focus:main',
                            ],

    },
    zip_safe=False,

)
