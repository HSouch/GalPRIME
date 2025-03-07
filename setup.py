from setuptools import setup, find_packages


__version__ = "3.0.1"

setup(
    name="galprime",
    version=__version__,
    description="GALaxy Profile Recovery from Images of Model Emission",
    author="Harrison Souchereau",
    author_email="harrison.souchereau@yale.edu",
    url="https://github.com/hsouch/galprime",
    packages=find_packages(),
    install_requires=[
        "astropy>=5",
        "photutils>=1.10",
        "matplotlib>=3.8",
        "numpy>=1.26",
        "scikit-image>=0.19",
        "tqdm",
        "configobj",
        "pebble",
    ],

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    scripts=[
        "bin/run_galprime",
        "bin/bg_selector",
        "bin/gp_combine"
    ],
)

