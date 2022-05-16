#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


setup(
    author="Yuxuan Zhuang",
    author_email='yuxuan.zhuang@dbb.su.se',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="test for parallel analysis for ensemble simulations",
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    name='ENPMDATests',
    packages=find_packages(),
    package_dir={'ENPMDATests': 'ENPMDATests'},
    install_requires=[
        'pytest>=3.3.0',
        'hypothesis',
    ],
    url='https://github.com/yuxuanzhuang/ENPMDA',
    version='0.1.0',
    zip_safe=False,
)
