── R CMD check results ─────────────────────────────────── estimators 0.8.3 ────
Duration: 2m 30.1s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔


── Windows Devel Check ───────────────────────────────────────────────────────── 

0 errors ✔ | 0 warnings ✔ | 1 note 

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2024-04-23 as issues were not corrected in time.

These issues have been fixed. Specifically, the error was due to an automated
test in which I did not include a `set.seed()` for reproducibility, allowing a
very small but positive probability of failure in the test. It is now solved.

── Rhub Check ────────────────────────────────────────────────────────────────── 

The following Rhub checks have been executed successfully: 

1 [VM] linux          R-* (any version)             ubuntu-latest on GitHub
2 [VM] macos          R-* (any version)             macos-13 on GitHub
3 [VM] macos-arm64    R-* (any version)             macos-latest on GitHub
4 [VM] windows        R-* (any version)             windows-latest on GitHub
5 [CT] atlas          R-devel (2024-05-04 r86521)   Fedora Linux 38
