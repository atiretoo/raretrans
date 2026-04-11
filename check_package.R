# =============================================================
# raretrans — CRAN pre-submission check script
# Run this from the package root directory
# =============================================================

# ---- 1. Housekeeping ----------------------------------------
cat("===== Step 1: Clean up old artefacts =====\n")
unlink("raretrans/", recursive = TRUE)
unlink("doc",        recursive = TRUE)
unlink("Meta",       recursive = TRUE)

# ---- 2. Document with roxygen2 ------------------------------
cat("\n===== Step 2: Roxygenise =====\n")
devtools::document()

# ---- 3. Basic devtools::check() (== R CMD check --as-cran) --
cat("\n===== Step 3: devtools::check (--as-cran) =====\n")
chk <- devtools::check(cran = TRUE, manual = TRUE, vignettes = TRUE)
print(chk)

# ---- 4. Build the tarball ------------------------------------
cat("\n===== Step 4: Build source tarball =====\n")
pkg_tar <- devtools::build()
cat("Tarball:", pkg_tar, "\n")

# ---- 5. R CMD check on the tarball (strictest) ---------------
cat("\n===== Step 5: R CMD check on tarball =====\n")
rcmdcheck::rcmdcheck(pkg_tar, args = c("--as-cran", "--no-manual"),
                     error_on = "warning")

# ---- 6. Spell check -----------------------------------------
cat("\n===== Step 6: Spell check =====\n")
spell <- devtools::spell_check()
print(spell)

# ---- 7. Check URLs / DOIs -----------------------------------
cat("\n===== Step 7: URL checks =====\n")
urlchecker::url_check()

# ---- 8. Win-builder (optional, uncomment to submit) ----------
# cat("\n===== Step 8: Win-builder =====\n")
# devtools::check_win_devel()
# devtools::check_win_release()

# ---- 9. Summary ----------------------------------------------
cat("\n===== Done =====\n")
cat("Errors:  ", length(chk$errors),   "\n")
cat("Warnings:", length(chk$warnings), "\n")
cat("Notes:   ", length(chk$notes),    "\n")
cat("\nIf 0 errors, 0 warnings, 0 notes — ready for CRAN!\n")
cat("Submit with:  devtools::release()\n")
