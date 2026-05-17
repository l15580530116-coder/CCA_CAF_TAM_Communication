# GitHub Release Checklist

**Date**: 2026-05-18

## Pre-Release Checks

| # | Check | Status |
|---|-------|--------|
| 1 | README.md exists | YES |
| 2 | LICENSE exists (MIT) | YES |
| 3 | .gitignore exists | YES |
| 4 | scripts/ copied (47 files: R + Python) | YES |
| 5 | gene_sets/ copied (5 files) | YES |
| 6 | environment/ with requirements + versions | YES |
| 7 | docs/ with workflow, data, reproducibility | YES |
| 8 | No raw data copied (data/ excluded) | YES |
| 9 | No submission_CSBJ*/ folders copied | YES |
| 10 | No author private documents | YES |
| 11 | No large RDS/TAR/GZ files | YES |
| 12 | No .Rhistory, .RData, __pycache__ | YES |

## Post-Release

| # | Action |
|---|--------|
| 13 | Run `git init` and commit |
| 14 | Create GitHub repository |
| 15 | Push to GitHub |
| 16 | Update manuscript Data/Code Availability with repository URL |
| 17 | Generate Zenodo DOI for repository snapshot |

## Repository Size Check

Estimated size: < 500 KB (scripts + docs only, no data).

---

*Ready for git init and push.*
