import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="dicomtools",
  version="0.0.4",
  author="Ryan Bosca",
  description="Some helpful tools for working with DICOM files",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/rjbosca/dicomtools",
  packages=setuptools.find_packages(),
  install_requires=[
    "pydicom"
  ],
  classifiers=[
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
  ]
)