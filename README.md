# HDMMODE

Materials for Multiobjective Differential Evolution for Higher-Dimensional Multimodal Multiobjective Optimization.

## Files

| File | Description |
|---|---|
| [HDMMODE.zip](HDMMODE.zip) | HDMMODE archive |
| [HDMMF.zip](HDMMF.zip) | HDMMF archive |
| [JAS-2024-0096_print.pdf](JAS-2024-0096_print.pdf) | Paper PDF |

## Reading the materials

This page indexes the files currently stored in the repository. The archive inventory and selected entry points are described below; runtime dependencies remain unverified. GitHub file search does not search inside ZIP archives; consult the archive contents and accompanying documents before running code.

The original packages and documents remain the source materials. If browsable source files are added later, identify the archive version they came from.

## Archive layout inspected on 2026-09-15

The archive directories and selected entry-point text were inspected without running experiments:

| Archive | Internal files | Entry points / contents |
|---|---:|---|
| `HDMMODE.zip` | 11 | MATLAB algorithm sources; `HDMMODE.m` defines a class derived from `ALGORITHM` |
| `HDMMF.zip` | 30 | 15 MATLAB problem classes and 15 reference MAT files; `HDMMF/NMMF1.m` derives from `PROBLEM` |

The sources use PlatEMO-style `ALGORITHM` / `PROBLEM` interfaces. The compatible platform version, installation and numerical results have not been validated by this index. The problem files load the accompanying reference MAT files, which are retained with the package.

Paths above are inside the original archives. Matching reference data in another benchmark package does not establish identical problem source code or justify removing a package version.
