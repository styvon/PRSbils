#!/usr/bin/env python3

import sys
from setuptools import setup, find_packages

VERSION = '1.0.0'  # keep same with PRSbils/__init__.py

if sys.version_info < (3, 6):
    sys.exit('Python >=3.6 is required for PRSbils')

with open('README.md', encoding="utf8") as f:
    readme = f.read().split('--------------------')[-1]

with open('requirements.txt') as f:
    rq = []
    for line in f:
        line = line.strip()
        rq.append(line.split('==')[0])


if __name__ == '__main__':
    setup(
        name='prsbils',
        version=VERSION,
        url='http://github.com/styvon/PRSbils',
        author='Yongwen Zhuang',
        author_email='zyongwen@umich.edu',
        license='LICENSE',
        description='Polygenic risk score with bilevel continuous shrinkage',
        long_description=readme,
        long_description_content_type='text/markdown',
        url='http://parl.ai/',
        python_requires='>=3.6',
        packages=find_packages(exclude=('data', 'docs', 'test')),
        install_requires=reqs,
        include_package_data=True,
        package_data={'': ['*.txt', '*.md', '*.opt']},
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
        test_suite='nose.collector',
        tests_require=['nose'],
    )