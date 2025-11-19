#######################################################
# Helper scripts to fit synthetic controls to simulations
#######################################################

#' Make a V matrix from a vector (or null)
make_V_matrix <- function(t0, V) {
  if(is.null(V)) {
        V <- diag(rep(1, t0))
    } else if(is.vector(V)) {
        if(length(V) != t0) {
          stop(paste("`V` must be a vector with", t0, "elements or a", t0, 
                     "x", t0, "matrix"))
        }
        V <- diag(V)
    } else if(ncol(V) == 1 & nrow(V) == t0) {
        V <- diag(c(V))
    } else if(ncol(V) == t0 & nrow(V) == 1) {
        V <- diag(c(V))
    } else if(nrow(V) == t0) {
    } else {
        stop(paste("`V` must be a vector with", t0, "elements or a", t0, 
                     "x", t0, "matrix"))
    }

  return(V)
}

#' Fit synthetic controls on outcomes after formatting data
#' @param synth_data Panel data in format of Synth::dataprep
#' @param V Matrix to scale the obejctive by
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Synth weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#' }
fit_synth_formatted <- function(synth_data, V = NULL, warm_start = NULL) {

  t0 <- dim(synth_data$Z0)[1]
  V  <- make_V_matrix(t0, V)

  sol <- synth_qp(
    X1 = synth_data$X1,
    X0 = t(synth_data$X0),
    V  = V,
    warm_start = warm_start
  )

  weights <- sol$x

  l2_imbalance <- sqrt(sum((synth_data$Z0 %*% weights - synth_data$Z1)^2))
  uni_w <- matrix(1/ncol(synth_data$Z0), nrow=ncol(synth_data$Z0), ncol=1)
  unif_l2_imbalance <- sqrt(sum((synth_data$Z0 %*% uni_w - synth_data$Z1)^2))
  scaled_l2_imbalance <- l2_imbalance / unif_l2_imbalance

  return(list(
    weights = weights,
    warm_status = sol$warm_status,
    l2_imbalance = l2_imbalance,
    scaled_l2_imbalance = scaled_l2_imbalance
  ))
}



#' Solve the synth QP directly
#' @param X1 Target vector
#' @param X0 Matrix of control outcomes
#' @param V Scaling matrix
#' @noRd
synth_qp <- function(X1, X0, V, warm_start = NULL) {

  Pmat <- X0 %*% V %*% t(X0)
  qvec <- - t(X1) %*% V %*% t(X0)

  n0 <- nrow(X0)

  A <- rbind(rep(1, n0), diag(n0))
  l <- c(1, numeric(n0))
  u <- c(1, rep(1, n0))

  # Quiet OSQP (no iteration spam)
  settings <- osqp::osqpSettings(
    verbose = FALSE,   # NO printing
    eps_rel = 1e-6,    # looser tolerance so fewer checks
    eps_abs = 1e-6
  )

  # Build solver model
  model <- osqp::osqp(
    P    = Pmat,
    q    = as.numeric(qvec),
    A    = A,
    l    = l,
    u    = u,
    pars = settings
  )

  ###############################################################
  #        UNIVERSAL WARM-START + ACCEPTANCE DETECTOR
  ###############################################################
  warm_status <- "NOT_USED"

  if (!is.null(warm_start)) {

    warm_start <- as.numeric(warm_start)

    if (length(warm_start) != n0) {
      warm_status <- paste0("REJECTED_LENGTH(", length(warm_start), " vs ", n0, ")")

    } else {

      # ---- Try newest API first ----
      applied <- FALSE
      if ("WarmStart" %in% names(model)) {
        model$WarmStart(x = warm_start)
        applied <- TRUE
      } else if ("warm_start" %in% names(model)) {
        model$warm_start(x = warm_start)
        applied <- TRUE
      } else if ("warm_start_x" %in% names(model)) {
        model$warm_start_x(x = warm_start)
        applied <- TRUE
      } else {
        warm_status <- "REJECTED_NO_API"
      }

      if (applied) {
        # check whether warm start stuck
        internal_x <- tryCatch(model$x, error = function(e) NULL)

        if (!is.null(internal_x) && length(internal_x) == length(warm_start)) {
          if (max(abs(internal_x - warm_start)) < 1e-6) {
            warm_status <- "ACCEPTED"
          } else {
            warm_status <- "REJECTED_SOLVER_MODIFIED_X"
          }
        } else {
          warm_status <- "REJECTED_INTERNAL_EMPTY"
        }
      }
    }
  }

  ###############################################################

  sol <- model$Solve()

  return(list(
    x = sol$x,
    warm_status = warm_status
  ))
}




