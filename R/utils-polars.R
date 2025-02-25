pl <- NULL

use_polars <- function() {
    if (is.null(pl)) {
        orepos <- getOption("repos")
        options(repos = c("https://community.r-multiverse.org", orepos))
        on.exit(options(repos = orepos))
        rlang::check_installed("polars", "to use `rsahmi` package")
        pl <<- polars::pl
    }
}

#' @param .x A series object.
#' @param .fn Must call `$collect_in_background()` method.
#' @keywords internal
#' @noRd
polars_lapply <- function(.x, .fn, ..., .progress, .threads = 2L) {
    # .threads <- min(as.integer(.threads), pl$thread_pool_size())
    poos_size <- max(.threads, 1L)
    n <- .x$len()
    bar <- do.call(cli::cli_progress_bar, c(
        .progress, list(total = n, clear = FALSE)
    ))
    results <- vector("list", n)
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
            results[[i]] <- .fn(.x$slice(i - 1L, 1L), ...)
            pools[pool] <- i
            i <- i + 1L
        }

        # if there is a task in this pool
        # we check if we can release this pool
        if (!is.na(handle_index)) {
            polars_handle <- .subset2(results, handle_index)
            if (polars_handle$is_finished()) {
                # collect result from this pool and release this pool
                results[[handle_index]] <- polars_handle$join()
                pools[pool] <- NA_integer_
                cli::cli_progress_update(inc = 1L, id = bar)

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
    results
}
