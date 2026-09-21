# HDMMODE

<!-- solver-policy-20260922 -->
> **2026-09-22 求解器决定：** 今后不再使用 Gurobi，也不再要求许可证或续期。采用当前项目已验证的替代器；尚未迁移的旧入口保持停用。历史结果及求解器标注保留。本段即本仓当前求解器约束；不改写历史实验记录。
<!-- /solver-policy-20260922 -->

**[逐文件路径与分类](FILEMAP.md) · [机器可读清单](FILEMAP.csv)**

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

## Directory roles and external dependencies (2026-09-21)

| Location | Role |
|---|---|
| [source/HDMMODE/](source/HDMMODE/) | Algorithm class and selection/operator/helper browsing copies |
| [source/HDMMF/HDMMF/](source/HDMMF/HDMMF/) | 15 benchmark problem classes |
| `HDMMF.zip/HDMMF/*_Reference_PSPF_data.mat` | 15 reference datasets, kept in the original archive alongside their corresponding class versions |
| ZIP files / paper PDF | Frozen source/data releases and publication material |

[HDMMODE.m](source/HDMMODE/HDMMODE.m) inherits from external `ALGORITHM`; the benchmark classes inherit from external `PROBLEM`. Neither platform base class is supplied here. For a future run, use a compatible PlatEMO environment and retain each MAT file beside its matching extracted problem class on the MATLAB path. [NMMF1.m](source/HDMMF/HDMMF/NMMF1.m) loads its data by filename at lines 53, 59 and 67; `source/` alone does not contain these data. Helper calls such as `pdist2` in [Crowding.m](source/HDMMODE/Crowding.m) also require a compatible MATLAB environment. No exact platform/toolbox version has been established.

All 41 ZIP members passed CRC/path/index-hash checks; the 26 committed source copies match their original member bytes. No MATLAB execution or result validation was performed. This release and a later benchmark package remain separate versions even when reference MAT bytes match; no live project-memory file is implied by this archive.
