── R CMD check results ─────────────────────────────────── estimators 0.8.3 ────
Duration: 5m 20.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔


── Archived Package ──────────────────────────────────────────────────────────── 

CRAN repository db overrides:
     X-CRAN-Comment: Archived on 2024-04-23 as issues were not corrected
       in time.

- Have these issues been fixed?
- Yes, the error was due to an automated test in which I did not include a
`set.seed()` for reproducibility, allowing a very small but positive probability
of failure in the test. It is now solved.
