# Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2020-01-31

Third unstable version.

### Changed
 - Established new version scheme with major, minor and patch number
 - Major revision of project structure and development process according to
   the bertha project (CI support, contributing and coding style guidelines,
   preparation for deployment, ...)
 - Major refactoring of solver classes (including separation of treatment
   of electromagnetic field and quantum mechanical systems)

### Removed
 - Dropped support for Intel Xeon Phi KNC
 - Dropped Predictor-Corrector and Runge-Kutta methods
 - Dropped support for MATLAB data export

## [0.2] - 2018-09-27

Second unstable version.

### Added
 - General description of materials, in particular quantum mechanical
   properties
 - Alternative numerical methods such as those based on the Rodrigues approach
 - Writer for HDF5 data files

## [0.1] - 2017-07-18

First unstable version, served as proof of concept.

### Added
 - CMake build system, modular project structure
 - Base classes that describe the device and simulation scenario
 - Elementary numerical methods
 - Parallelization support for OpenMP and CUDA
 - Writer for MATLAB data files
