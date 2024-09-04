import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = {}
with open("obstools_ipgp/version.py") as fp:
    exec(fp.read(), version)

setuptools.setup(
    name="obstools_ipgp",
    version=version['__version__'],
    author="Wayne Crawford",
    author_email="crawford@ipgp.fr",
    description="Ocean bottom seismometer evaluation/processing routines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WayneCrawford/obstools_ipgp",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
          'obspy>=1.1',
          'pyyaml>5.0',
          'jsonschema>=2.6',
          'jsonref>=0.2',
          'progress>=1.5',
          'tiskitpy>=0.4'
      ],
    entry_points={
         'console_scripts': [
             'plotPSDs=obstools_ipgp.plotPSDs:main',
             'obstest=obstools_ipgp.obstest:main'
         ]
    },
    python_requires='>=3.8',
    classifiers=(
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ),
    keywords='oceanography, marine, OBS'
)
