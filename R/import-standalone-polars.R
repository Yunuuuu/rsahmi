# Standalone file: do not edit by hand
# Source: <https://github.com/Yunuuuu/standalone/blob/main/R/standalone-polars.R>
# ----------------------------------------------------------------------
#
# ---
# repo: Yunuuuu/standalone
# file: standalone-polars.R
# last-updated: 2025-02-26
# license: https://unlicense.org
# dependencies: [standalone-pkg.R]
# imports: [rlang]
# ---

# This file contains several helper functions for working with the polars
# package.

# Note: Since `polars` is not available on `CRAN`, if you plan to publish your
# package on `CRAN`, please include the following in your `DESCRIPTION` file:
# `Additional_repositories: https://community.r-multiverse.org` and run
# `usethis::use_package("polars", "Suggests")`. For each user-defined function
# that uses polars, be sure to add `use_polars()` as the first line in the
# function

# ## Changelog
# 2025-02-26:
# - Add `install_polars`
# - Add `use_polars`
# - Add `series_lapply`
#
# nocov start

pl <- NULL

# Helper function to install `polars` package
#
# Helper function to install `polars` package from
# [`r-multiverse`](https://community.r-multiverse.org).
install_polars <- function() {
    orepos <- getOption("repos")
    options(repos = c("https://community.r-multiverse.org", orepos))
    on.exit(options(repos = orepos), add = TRUE)
    install_pkgs("polars")
}

use_polars <- function(reason) {
    if (is.null(pl)) {
        if (!is_installed("polars")) {
            if (missing(reason)) {
                reason <- sprintf("to use `%s` package", pkg_nm())
            }
            orepos <- getOption("repos")
            options(repos = c("https://community.r-multiverse.org", orepos))
            on.exit(options(repos = orepos), add = TRUE)
            rlang::check_installed("polars", reason = reason)
        }
        utils::assignInNamespace("pl",
            getExportedValue("polars", "pl"),
            ns = pkg_namespace()
        )
    }
}

#' @param .x A [Series][polars::Series] object.
#' @param .fn Must call `$collect_in_background()` method.
#' @param ... Additional arguments passed to `.fn`.
#' @param .progress Additional arguments passed to [cli::cli_progress_bar].
#' @keywords internal
#' @noRd
series_lapply <- function(.x, .fn, ..., .progress, .threads = 1L) {
    # .threads <- min(as.integer(.threads), pl$thread_pool_size())
    poos_size <- max(.threads, 1L)
    n <- .x$len()
    if (is_installed("cli")) {
        bar <- rlang::inject(
            getExportedValue("cli", "cli_progress_bar")(
                !!!.progress, total = n, clear = FALSE
            )
        )
    }
    out <- vector("list", n)
    # pools: store the index of result
    # NA means this pool can be used
    pools <- rep_len(NA_integer_, poos_size)
    i <- pool <- 1L
    while (i <= n || !all(is.na(pools))) {
        handle_index <- .subset(pools, pool)
        # if i > n, we skip add task
        if (i <= n && is.na(handle_index)) {
            # this pool can be used, we add task into this pool
            # `.fn()` must return with `$collect_in_background()` method
            # For Series, cannot subset with `[[`
            out[[i]] <- .fn(.x$slice(i - 1L, 1L), ...)
            pools[pool] <- i
            i <- i + 1L
        }

        # if there is a task in this pool
        # we check if we can release this pool
        if (!is.na(handle_index)) {
            polars_handle <- .subset2(out, handle_index)
            if (polars_handle$is_finished()) {
                # collect result from this pool and release this pool
                out[[handle_index]] <- polars_handle$join()
                pools[pool] <- NA_integer_
                if (is_installed("cli")) {
                    getExportedValue("cli", "cli_progress_update")(
                        inc = 1L, id = bar
                    )
                }

                # this pool has been released, so we directly
                # step into next cycle and re-use this pool
                next
            }
        }

        # search next pool
        if (pool == poos_size) {
            pool <- 1L
        } else {
            pool <- pool + 1L
        }
    }
    out
}

# nocov end
