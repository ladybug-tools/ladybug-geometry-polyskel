import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('cli-requirements.txt') as f:
    cli_requirements = f.read().splitlines()

setuptools.setup(
    name="ladybug-geometry-polyskel",
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author="Ladybug Tools",
    author_email="info@ladybug.tools",
    description="A library with poly skeleton methods using ladybug-geometry",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ladybug-tools/ladybug-geometry-polyskel",
    packages=setuptools.find_packages(exclude=["tests"]),
    install_requires=requirements,
    extras_require={'cli': cli_requirements},
    entry_points={
        "console_scripts": ["ladybug-geometry-polyskel = ladybug_geometry_polyskel.cli:main"]
    },
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: IronPython",
        "Operating System :: OS Independent"
    ],
    license="AGPL-3.0"
)
