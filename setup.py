"""Minimal setup file for tasks project."""

from setuptools import setup, find_packages

setup(
    name='nw_setup',
    version='0.3.0',
    license='proprietary',
    description='Module Experiment',

    author='hsasaki',
    author_email='hsasaki@softmatters.net',
    url='https://github.com/softmatter-design/python-cognac-nw_setup/',

    packages=find_packages(where='src'),
    package_dir={'': 'src'},

    entry_points={
        "console_scripts": [
          'nw_setup = nw_setup.nw_setup:main',
          'nw_setup_multi = nw_setup.nw_setup_multi:main',
          'evaluate_nw = chain_evaluation.evaluate_all:evaluate'
        ]
    }
)
