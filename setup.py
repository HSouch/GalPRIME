import galprime
from setuptools import setup, find_packages

setup(
    name='galprime',
    version=galprime.__version__,
    description='GALaxy Profile Recovery from Images of Model Emission',
    author='Harrison Souchereau',
    author_email='harrison.souchereau@yale.edu',
    url='https://github.com/hsouch/galprime',
    packages=find_packages(),
    install_requires=[
        # Add your dependencies here
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
)
