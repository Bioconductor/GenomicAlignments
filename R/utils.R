### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


### 3 equivalent implementations for this:
###   (a) x %in% x[duplicated(x)]
###   (b) duplicated(x) | duplicated(x, fromLast=TRUE)
###   (c) xx <- match(x, x); ans <- xx != seq_along(xx); ans[xx] <- ans; ans
### Comparing the 3 implementations on an integer vector of length 12 millions:
###   (a) is the most memory efficient;
###   (b) is a little bit faster than (a) (by only 8%) but uses between 12-14%
###       more memory;
###   (c) is as fast as (a) but uses about 30% more memory.
has_duplicates <- function(x)
{
    x %in% x[duplicated(x)]
}

