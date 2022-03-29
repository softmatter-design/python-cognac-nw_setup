"""Minimal setup file for tasks project."""

from setuptools import setup, find_packages

setup(
    name='nw_setup',
    version='0.1.2',
    license='proprietary',
    description='Module Experiment',

    author='hsasaki',
    author_email='hsasaki@softmatters.net',
    url='None.com',

    packages=find_packages(where='src'),
    package_dir={'': 'src'},

    entry_points={
        "console_scripts": [
          'nw_setup = nw_setup.nw_setup:main'
        ]
    }
)