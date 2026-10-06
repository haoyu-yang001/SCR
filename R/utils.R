# Internal helpers shared by SCR(), scr_next(), smart() and scr_imbalance().

# Coerce continuous covariates to a numeric matrix (NULL allowed).
.as_covariate_matrix <- function(covariate, name = "covariate", allow_null = TRUE) {
  if (is.null(covariate)) {
    if (!allow_null) stop(sprintf("`%s` must not be NULL.", name), call. = FALSE)
    return(NULL)
  }
  if (is.data.frame(covariate)) {
    ok <- vapply(covariate, function(col) is.numeric(col) || is.logical(col), logical(1))
    if (!all(ok)) {
      stop(sprintf(
        "`%s` must contain numeric columns only; pass categorical covariates through `categorical`. Non-numeric column(s): %s",
        name, paste(names(covariate)[!ok], collapse = ", ")
      ), call. = FALSE)
    }
    covariate <- as.matrix(covariate)
  }
  if (is.null(dim(covariate))) covariate <- matrix(covariate, ncol = 1)
  if (!is.numeric(covariate) && !is.logical(covariate)) {
    stop(sprintf("`%s` must be a numeric matrix or data frame.", name), call. = FALSE)
  }
  storage.mode(covariate) <- "double"
  if (anyNA(covariate)) stop(sprintf("`%s` contains missing values.", name), call. = FALSE)
  covariate
}

# Coerce categorical covariates to a data frame (NULL allowed).
.as_categorical_df <- function(categorical, name = "categorical") {
  if (is.null(categorical)) return(NULL)
  categorical <- as.data.frame(categorical, stringsAsFactors = FALSE)
  if (anyNA(categorical)) stop(sprintf("`%s` contains missing values.", name), call. = FALSE)
  categorical
}

# Integer stratum code for each row: strata are the observed combinations of
# all categorical covariates (within-stratum imbalance, Section 2.2).
.strata_codes <- function(categorical) {
  if (is.null(categorical)) return(NULL)
  key <- do.call(paste, c(lapply(categorical, as.character), sep = "\r"))
  as.integer(factor(key))
}

# Validate an (optional) vector of prior assignments coded 0/1.
.check_assignment <- function(assignment, n, levels = c(0L, 1L)) {
  if (is.null(assignment) || (length(assignment) == 1L && is.na(assignment))) {
    return(integer(0))
  }
  if (is.data.frame(assignment)) assignment <- assignment[[ncol(assignment)]]
  if (anyNA(assignment)) stop("`assignment` contains missing values.", call. = FALSE)
  if (length(assignment) > n) {
    stop("`assignment` is longer than the number of units.", call. = FALSE)
  }
  if (!all(assignment %in% levels)) {
    if (identical(levels, c(0L, 1L)) && all(assignment %in% c(1, 2))) {
      stop("`assignment` must be coded 0/1 (1 = treatment, 0 = control). ",
           "SCR 1.0.0 used 1/2 coding as input; use `assignment - 1` to convert.",
           call. = FALSE)
    }
    stop(sprintf("`assignment` must take values in {%s}.", paste(levels, collapse = ", ")),
         call. = FALSE)
  }
  as.integer(assignment)
}

.check_prob <- function(q, name) {
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < 0.5 || q > 1) {
    stop(sprintf("`%s` must be a single number in [0.5, 1].", name), call. = FALSE)
  }
}

.check_d <- function(d) {
  if (!is.numeric(d) || length(d) != 1L || is.na(d) || d < 0) {
    stop("`d` must be a single non-negative number.", call. = FALSE)
  }
}

.check_w <- function(w) {
  if (!is.numeric(w) || length(w) != 1L || is.na(w) || w < 0 || w > 1) {
    stop("`w` must be a single number in [0, 1].", call. = FALSE)
  }
}

# Running mean and scatter matrix (Welford), so cov(x_1, ..., x_i) is
# available at every step without recomputing from scratch.
.welford_add <- function(state, x) {
  state$k <- state$k + 1L
  delta <- x - state$mean
  state$mean <- state$mean + delta / state$k
  state$M2 <- state$M2 + tcrossprod(delta, x - state$mean)
  state
}

.welford_cov <- function(state) {
  if (state$k < 2L) return(matrix(0, length(state$mean), length(state$mean)))
  state$M2 / (state$k - 1L)
}

# Draw a 0/1 assignment with P(T = 1) = prob. Exactly one uniform is drawn per
# unit (even when prob is 0 or 1) so that SCR() and repeated scr_next() calls
# consume the random number stream identically.
.draw <- function(prob) as.integer(stats::runif(1) < prob)

# Is an imbalance difference numerically zero?
.is_tie <- function(diff, scale) abs(diff) <= 1e-10 * (1 + scale)
