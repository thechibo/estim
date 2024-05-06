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
