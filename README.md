[![Build Status](https://travis-ci.com/ladybug-tools/ladybug-geometry-polyskel.svg?branch=master)](https://travis-ci.com/ladybug-tools/ladybug-geometry-polyskel)
[![Coverage Status](https://coveralls.io/repos/github/ladybug-tools/ladybug-geometry-polyskel/badge.svg?branch=master)](https://coveralls.io/github/ladybug-tools/ladybug-geometry-polyskel)

[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Python 2.7](https://img.shields.io/badge/python-2.7-green.svg)](https://www.python.org/downloads/release/python-270/) [![IronPython](https://img.shields.io/badge/ironpython-2.7-red.svg)](https://github.com/IronLanguages/ironpython2/releases/tag/ipy-2.7.8/)

# ladybug-geometry-polyskel

A library with straight skeleton methods using ladybug-geometry.

## Credits

Ladybug-geometry-polyskel is a derivative work of the the [polyskel package](https://github.com/Botffy/polyskel)
by [Ármin Scipiades](https://github.com/Botffy) (@Bottfy), which is, itself, a Python 3
implementation of the straight skeleton algorithm as described by Felkel and Obdržálek
in their 1998 conference paper [Straight skeleton implementation](https://github.com/Botffy/polyskel/blob/master/doc/StraightSkeletonImplementation.pdf).

Key differences between Bottfy's original implementation and this package are:

* It has been modified for compatibility with both Python 2.7 and Python 3.7.
* The code as been re-stylized to conform to the PEP8 style guide.
* Modules have been added to extract core/perimeter polygons from the straight skeleton.

## Installation

```console
pip install -U ladybug-geometry-polyskel
```

## QuickStart

```python
import ladybug_geometry_polyskel

```

## [API Documentation](http://ladybug-tools.github.io/ladybug-geometry-polyskel/docs)

## Local Development
1. Clone this repo locally
```console
git clone git@github.com:ladybug-tools/ladybug-geometry-polyskel

# or

git clone https://github.com/ladybug-tools/ladybug-geometry-polyskel
```
2. Install dependencies:
```console
cd ladybug-geometry-polyskel
pip install -r dev-requirements.txt
pip install -r requirements.txt
```

3. Run Tests:
```console
python -m pytest tests/
```

4. Generate Documentation:
```console
sphinx-apidoc -f -e -d 4 -o ./docs ./ladybug_geometry_polyskel
sphinx-build -b html ./docs ./docs/_build/docs
```

## Copyright 

Ladybug Geometry Polyskel, Copyright (c) 2021, Ármin Scipiades, Ladybug Tools LLC
and other contributors. All rights reserved.
