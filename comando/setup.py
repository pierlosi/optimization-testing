"""Setup for the COMANDO package."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHOR: Marco Langiu
from setuptools import setup, find_packages


# Read the version from the back of 'comando/__init__.py'
with open('comando/__init__.py', 'rb') as f:
    f.seek(-42, 2)  # 42 is a good number!
    line = f.readlines()[-1].decode()  # sth. like "__version__ = 'a.b.c'\n"
    version = line.split("'")[1]

with open('README.md', 'r') as f:
    README = f.read()


# As proposed by Han Xiao in https://hanxiao.io/2019/11/07/
# A-Better-Practice-for-Managing-extras-require-Dependencies-in-Python/
def get_extra_requires(path, add_all=True):
    """Parse extra-requirements.txt for a {feature: requirements} map."""
    import re
    from collections import defaultdict

    with open(path) as fp:
        extra_deps = defaultdict(set)
        for k in fp:
            if k.strip() and not k.startswith('#'):
                tags = set()
                if ':' in k:
                    k, v = k.split(':')
                    tags.update(vv.strip() for vv in v.split(','))
                tags.add(re.split('[<=>]', k)[0])
                for t in tags:
                    extra_deps[t].add(k)

        # add tag `all` at the end
        if add_all:
            extra_deps['all'] = {vv for v in extra_deps.values() for vv in v}

    return extra_deps


setup(name='comando',
      version=version,
      author='Marco Langiu',
      author_email='m.langiu@fz-juelich.de',
      description='A next-generation modeling framework for component-oriented'
      ' modeling and optimization for nonlinear design and operation of '
      'integrated energy systems',
      long_description=README,
      long_description_content_type='text/markdown',
      license='MIT',
      url='https://jugit.fz-juelich.de/iek-10/core-projects/comando/comando',
      packages=find_packages(include=['comando', 'comando.*']),
      python_requires='>=3.7',
      install_requires=['pandas', 'sympy>=1.6'],
      extras_require=get_extra_requires('extra-requirements.txt'),
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Topic :: Scientific/Engineering'
      ])
