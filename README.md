[![Build Status](https://travis-ci.com/ladybug-tools/ladybug-geometry-polyskel.svg?branch=master)](https://travis-ci.com/ladybug-tools/ladybug-geometry-polyskel)
[![Coverage Status](https://coveralls.io/repos/github/ladybug-tools/ladybug-geometry-polyskel/badge.svg?branch=master)](https://coveralls.io/github/ladybug-tools/ladybug-geometry-polyskel)

[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Python 2.7](https://img.shields.io/badge/python-2.7-green.svg)](https://www.python.org/downloads/release/python-270/) [![IronPython](https://img.shields.io/badge/ironpython-2.7-red.svg)](https://github.com/IronLanguages/ironpython2/releases/tag/ipy-2.7.8/)

# ladybug-geometry-polyskel

A library with straight skeleton methods using ladybug-geometry.

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
