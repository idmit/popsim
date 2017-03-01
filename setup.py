# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license_ = f.read()

setup(
    name='popsim',
    version='0.1.0',
    description='A tool to simulate progeny from a plink ped file',
    long_description=readme,
    author='Ivan Dmitrievsky',
    author_email='ivan.dmitrievsky+opensource@gmail.com',
    url='https://gitlab.com/idmit/popsim',
    install_requires=[
        'click',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'popsim=popsim.core:root'
        ]
    },
    license=license_,
    packages=find_packages(exclude=('tests', 'docs'))
)
