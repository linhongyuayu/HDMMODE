# HDMMODE

Materials for Multiobjective Differential Evolution for Higher-Dimensional Multimodal Multiobjective Optimization.

## Files

| File | Description |
|---|---|
| [HDMMODE.zip](HDMMODE.zip) | HDMMODE archive |
| [HDMMF.zip](HDMMF.zip) | HDMMF archive |
| [JAS-2024-0096_print.pdf](JAS-2024-0096_print.pdf) | Paper PDF |

## Browsable source and complete member index

- [source/](source/): source browsing copies under `source/<archive-stem>/<original-member-path>`.
- [Complete member index](ARCHIVE_INDEX.md) / [JSON index](ARCHIVE_INDEX.json): all 41 archive members, their sizes and SHA256 hashes, including links to 26 source copies.

The original ZIPs are frozen artifacts; `source/` copies preserve the exact member bytes. Running the code still requires the matching archive data, working directory and dependencies; standalone execution has not been verified.

## Archive layout inspected on 2026-09-15

The archive directories and selected entry-point text were inspected without running experiments:

| Archive | Internal files | Entry points / contents |
|---|---:|---|
| `HDMMODE.zip` | 11 | MATLAB algorithm sources; `HDMMODE.m` defines a class derived from `ALGORITHM` |
| `HDMMF.zip` | 30 | 15 MATLAB problem classes and 15 reference MAT files; `HDMMF/NMMF1.m` derives from `PROBLEM` |

The sources use PlatEMO-style `ALGORITHM` / `PROBLEM` interfaces. The compatible platform version, installation and numerical results have not been validated by this index. The problem files load the accompanying reference MAT files, which are retained with the package.

Matching reference data in another benchmark package does not establish identical problem source code or justify removing a package version.
