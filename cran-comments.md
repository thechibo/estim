── R CMD check results ─────────────────────────────────── estimators 0.8.4 ────
Duration: 2m 30.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔


── Windows Devel Check ───────────────────────────────────────────────────────── 

0 errors ✔ | 0 warnings ✔ | 1 note 

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2024-04-23 as issues were not corrected in time.

These issues have been fixed. Specifically, the error was due to an automated
test in which I did not include a `set.seed()` for reproducibility, allowing a
very small but positive probability of failure in the test. It is now solved.

── Rhub Check ────────────────────────────────────────────────────────────────── 

The following checks have been executed successfully with `rhub::rhub_check()`: 

 1 [VM] linux          R-* (any version)                     ubuntu-latest on GitHub
 2 [VM] macos          R-* (any version)                     macos-13 on GitHub
 3 [VM] macos-arm64    R-* (any version)                     macos-latest on GitHub
 4 [VM] windows        R-* (any version)                     windows-latest on GitHub
 5 [CT] atlas          R-devel (2024-05-06 r86526)           Fedora Linux 38 (Container Image)
 6 [CT] clang-asan     R-devel (2024-05-06 r86526)           Ubuntu 22.04.4 LTS
 7 [CT] clang16        R-devel (2024-05-04 r86521)           Ubuntu 22.04.4 LTS
 8 [CT] clang17        R-devel (2024-05-04 r86521)           Ubuntu 22.04.4 LTS
 9 [CT] clang18        R-devel (2024-05-04 r86521)           Ubuntu 22.04.4 LTS
10 [CT] donttest       R-devel (2024-05-04 r86521)           Ubuntu 22.04.4 LTS
11 [CT] gcc13          R-devel (2024-05-06 r86526)           Fedora Linux 38 (Container Image)
12 [CT] intel          R-devel (2024-05-06 r86526)           Fedora Linux 38 (Container Image)
13 [CT] mkl            R-devel (2024-05-06 r86526)           Fedora Linux 38 (Container Image)
14 [CT] nold           R-devel (2024-05-06 r86526)           Ubuntu 22.04.4 LTS
16 [CT] ubuntu-clang   R-devel (2024-05-06 r86526)           Ubuntu 22.04.4 LTS
17 [CT] ubuntu-gcc12   R-devel (2024-05-06 r86526)           Ubuntu 22.04.4 LTS
18 [CT] ubuntu-next    R-4.4.0 (patched) (2024-05-06 r86526) Ubuntu 22.04.4 LTS
19 [CT] ubuntu-release R-4.4.0 (2024-04-24)                  Ubuntu 22.04.4 LTS
20 [CT] valgrind       R-devel (2024-05-06 r86526)           Fedora Linux 38 (Container Image)
