import os
from setuptools import setup, find_packages

setup(
    name='homer',
    description='CSE185 Project',
    author='Hwa Yeon Lee',
    author_email='hyl023@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "myhomer=myhomer.myhomer:main"
        ],
    },
)