# Changelog

## [4.1.0] - 2020-11-07

### Added
  - Added option `-dd-recursive` to enable a recursive search for FASTQ files in the data-directory.
  - Added `--log-level` option to control the verbosity of log messages.
  - Added verbose versions of command-line options, e.g. `--data-directory` for `-dd`.
  - Added the ability to control output formats: Currently XLSX and/or JSON.

### Removed
  - Removed the `-do` (output directory) command-line option.

### Changed
  - The location of all output files and folders are now specified using a positional argument representing a prefix for those files. The prefix defaults to `output`.
  - Changed from one-file script to a more complex script structure. A `setup.py` file was added to allow easy installation of the script and its dependencies.

### Fixed
  - Fixed `-pl`/`-plate-layout` not being used in the calculation of well positions.


## 4.0.0 - 2014-04-16

Initial release.


[Unreleased]: https://github.com/laeblab/hamplicons/compare/v4.0.1...HEAD
[4.0.1]: https://github.com/laeblab/hamplicons/compare/v4.0.0...v4.0.1
