#modified random effects function from metafor
# developed by Dr. Chris Fleming Jan 2023

ranef2 <- function(object, level, digits, transf, targs, verbose=FALSE, ...)
{
 mstyle <- metafor:::.get.mstyle("crayon" %in% .packages())

 metafor:::.chkclass(class(object), must="rma.mv")

 x <- object

 na.act <- getOption("na.action")

 if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
  stop(mstyle$stop("Unknown 'na.action' specified under options()."))

 if (missing(level))
  level <- x$level

 if (missing(digits)) {
  digits <- metafor:::.get.digits(xdigits=x$digits, dmiss=TRUE)
 } else {
  digits <- metafor:::.get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
 }

 if (missing(transf))
  transf <- FALSE

 if (missing(targs))
  targs <- NULL

 level <- metafor:::.level(level)

 if (x$test == "z")
  crit <- qnorm(level/2, lower.tail=FALSE)

 ### TODO: check computations for user-defined weights

 if (!is.null(x$W))
  stop(mstyle$stop("Extraction of random effects not available for models with non-standard weights."))

 ddd <- list(...)

 metafor:::.chkdots(ddd, c("expand"))

 expand <- ifelse(is.null(expand), FALSE, isTRUE(ddd$expand)) # TODO: make this an option?

 #########################################################################

 out <- NULL

 if (verbose)
  message(mstyle$message("\nComputing inverse marginal var-cov and hat matrix ... "), appendLF = FALSE)

 ### compute inverse marginal var-cov and hat matrix

 W     <- chol2inv(chol(x$M))
 stXWX <- chol2inv(chol(as.matrix(t(x$X) %*% W %*% x$X)))
 Hmat  <- x$X %*% stXWX %*% crossprod(x$X,W)

 if (verbose)
  message(mstyle$message("Done!"))

 ### compute residuals

 ei <- c(x$yi - x$X %*% x$beta) ### use this instead of resid(), since this guarantees that the length is correct

 ### create identity matrix

 if (x$sparse) {
  I <- metafor:::Diagonal(x$k)
 } else {
  I <- diag(x$k)
 }

 if (x$withS) {

  # u^ = DZ'W(y - Xb) = DZ'We, where W = M^-1
  # note: vpred = var(u^ - u)

  out <- vector(mode="list", length=x$sigma2s)
  names(out) <- x$s.names

  for (j in seq_len(x$sigma2s)) {

   if (verbose)
    message(mstyle$message(paste0("Computing BLUPs for '", paste0("~ 1 | ", x$s.names[j]), "' term ... ")), appendLF = FALSE)

   if (x$Rfix[j]) {
    if (x$sparse) {
     D <- x$sigma2[j] * metafor:::Matrix(x$R[[j]], sparse=TRUE)
    } else {
     D <- x$sigma2[j] * x$R[[j]]
    }
   } else {
    if (x$sparse) {
     D <- x$sigma2[j] * metafor:::Diagonal(x$s.nlevels[j])
    } else {
     D <- x$sigma2[j] * diag(x$s.nlevels[j])
    }
   }

   DZtW <- D %*% t(x$Z.S[[j]]) %*% W
   pred <- as.vector(DZtW %*% cbind(ei))
   pred[abs(pred) < 100 * .Machine$double.eps] <- 0
   #vpred <- D - (DZtW %*% x$Z.S[[j]] %*% D - DZtW %*% x$X %*% stXWX %*% t(x$X) %*% W %*% x$Z.S[[j]] %*% D)
   vpred <- D - (DZtW %*% (I - Hmat) %*% x$Z.S[[j]] %*% D) # this one is the same as ranef.rma.uni() for standard RE/ME models
   #vpred <- DZtW %*% (I - Hmat) %*% x$Z.S[[j]] %*% D # = var(u^)
   #vpred <- D - (DZtW %*% x$Z.S[[j]] %*% D) # same as lme4::ranef()
   #vpred <- DZtW %*% x$Z.S[[j]] %*% D

   # if (is.element(x$test, c("knha","adhoc","t"))) {
   #  ddf <- metafor:::.ddf.calc(x$dfs, k=x$k, p=x$p, mf.s=x$mf.s[[j]], beta=FALSE)
   #  # crit <- qt(level/2, df=ddf, lower.tail=FALSE)
   # }

   # se <- sqrt(diag(vpred))
   # pi.lb <- c(pred - crit * se)
   # pi.ub <- c(pred + crit * se)

   # pred <- data.frame(intrcpt=pred, se=se, pi.lb=pi.lb, pi.ub=pi.ub)

   if (na.act == "na.omit") {

    names(pred) <- x$s.levels[[j]]
    dimnames(vpred) <- list(names(pred),names(pred))

    out[[j]] <- list(est=pred,COV=vpred)
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {

    ### determine which levels were removed

    s.levels.r <- !is.element(x$s.levels.f[[j]], x$s.levels[[j]])

    NAs <- rep(NA, x$s.nlevels.f[j])
    tmp <- NAs
    tmp[!s.levels.r] <- pred
    pred <- tmp

    names(pred) <- x$s.levels.f[[j]]

    out[[j]] <- pred

   }

   if (expand) {

    rows <- c(x$Z.S[[j]] %*% seq_along(x$s.levels[[j]]))
    pred <- pred[rows]
    rnames <- x$s.levels[[j]][rows]

    names(pred) <- metafor:::.make.unique(x$s.levels[[j]][rows])
    dimnames(vpred) <- list(names(pred),names(pred))

    out[[j]] <- list(est=pred,COV=vpred)
   }

   if (verbose)
    message(mstyle$message("Done!"))

  }

 }

 if (x$withG) {

  if (is.element(x$struct[1], c("GEN","GDIAG"))) {
   if (verbose)
    message(mstyle$message("Computation of BLUPs not currently available for struct=\"GEN\"."))
  } else {

   if (verbose)
    message(mstyle$message(paste0("Computing BLUPs for '", paste(x$g.names, collapse=" | "), "' term ... ")), appendLF = FALSE)

   G <- (x$Z.G1 %*% x$G %*% t(x$Z.G1)) * tcrossprod(x$Z.G2)
   GW <- G %*% W
   pred  <- as.vector(GW %*% cbind(ei))
   pred[abs(pred) < 100 * .Machine$double.eps] <- 0
   #vpred <- G - (GW %*% G - GW %*% x$X %*% stXWX %*% t(x$X) %*% W %*% G)
   vpred <- G - (GW %*% (I - Hmat) %*% G)

   # if (is.element(x$test, c("knha","adhoc","t"))) {
   #  ddf <- metafor:::.ddf.calc(x$dfs, k=x$k, p=x$p, mf.g=x$mf.g[[2]], beta=FALSE)
   #  # crit <- qt(level/2, df=ddf, lower.tail=FALSE)
   # }

   # se <- sqrt(diag(vpred))
   # pi.lb <- c(pred - crit * se)
   # pi.ub <- c(pred + crit * se)

   # pred <- data.frame(intrcpt=pred, se=se, pi.lb=pi.lb, pi.ub=pi.ub)

   nvars <- ncol(x$mf.g)

   if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
    r.names <- paste(formatC(x$ids[x$not.na], format="f", digits=0, width=max(nchar(x$ids[x$not.na]))), x$mf.g[[nvars]], sep=" | ")
   } else {
    #r.names <- paste(x$mf.g[[1]], x$mf.g[[2]], sep=" | ")
    r.names <- paste(sprintf(paste0("%", max(nchar(paste(x$mf.g[[1]]))), "s", collapse=""), x$mf.g[[1]]), x$mf.g[[nvars]], sep=" | ")
   }

   is.dup <- duplicated(r.names)

   pred <- pred[!is.dup]
   vpred <- vpred[!is.dup,!is.dup]

   names(pred) <- r.names[!is.dup]
   dimnames(vpred) <- list(names(pred),names(pred))

   if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
    #r.order <- order(x$mf.g[[nvars]][!is.dup], seq_len(x$k)[!is.dup])
    r.order <- seq_len(x$k)
   } else {
    r.order <- order(x$mf.g[[2]][!is.dup], x$mf.g[[1]][!is.dup])
   }

   pred <- pred[r.order]
   vpred <- vpred[r.order,r.order]

   out <- c(out, list(est=pred,COV=vpred))
   #names(out)[length(out)] <- paste(x$g.names, collapse=" | ")
   names(out)[length(out)] <- paste0(x$formulas[[1]], collapse="")

   if (verbose)
    message(mstyle$message("Done!"))

  }

 }

 if (x$withH) {

  if (is.element(x$struct[2], c("GEN","GDIAG"))) {
   if (verbose)
    message(mstyle$message("Computation of BLUPs not currently available for struct=\"GEN\"."))
  } else {

   if (verbose)
    message(mstyle$message(paste0("Computing BLUPs for '", paste(x$h.names, collapse=" | "), "' term ... ")), appendLF = FALSE)

   H <- (x$Z.H1 %*% x$H %*% t(x$Z.H1)) * tcrossprod(x$Z.H2)
   HW <- H %*% W
   pred  <- as.vector(HW %*% cbind(ei))
   pred[abs(pred) < 100 * .Machine$double.eps] <- 0
   #vpred <- H - (HW %*% H - HW %*% x$X %*% stXWX %*% t(x$X) %*% W %*% H)
   vpred <- H - (HW %*% (I - Hmat) %*% H)

   # if (is.element(x$test, c("knha","adhoc","t"))) {
   #  ddf <- metafor:::.ddf.calc(x$dfs, k=x$k, p=x$p, mf.h=x$mf.h[[2]], beta=FALSE)
   #  # crit <- qt(level/2, df=ddf, lower.tail=FALSE)
   # }

   # se <- sqrt(diag(vpred))
   # pi.lb <- c(pred - crit * se)
   # pi.ub <- c(pred + crit * se)

   # pred <- data.frame(intrcpt=pred, se=se, pi.lb=pi.lb, pi.ub=pi.ub)

   nvars <- ncol(x$mf.h)

   if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
    r.names <- paste(formatC(x$ids[x$not.na], format="f", digits=0, width=max(nchar(x$ids[x$not.na]))), x$mf.h[[nvars]], sep=" | ")
   } else {
    #r.names <- paste(x$mf.h[[1]], x$mf.h[[2]], sep=" | ")
    r.names <- paste(sprintf(paste0("%", max(nchar(paste(x$mf.h[[1]]))), "s", collapse=""), x$mf.h[[1]]), x$mf.h[[nvars]], sep=" | ")
   }

   is.dup <- duplicated(r.names)

   pred <- pred[!is.dup]
   vpred <- vpred[!is.dup,!is.dup]

   names(pred) <- r.names[!is.dup]
   dimnames(vpred) <- list(names(pred),names(pred))

   if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
    #r.order <- order(x$mf.h[[nvars]][!is.dup], seq_len(x$k)[!is.dup])
    r.order <- seq_len(x$k)
   } else {
    r.order <- order(x$mf.h[[2]][!is.dup], x$mf.h[[1]][!is.dup])
   }

   pred <- pred[r.order]
   vpred <- vpred[r.order,r.order]

   out <- c(out, list(est=pred,COV=vpred))
   #names(out)[length(out)] <- paste(x$h.names, collapse=" | ")
   names(out)[length(out)] <- paste0(x$formulas[[2]], collapse="")

   if (verbose)
    message(mstyle$message("Done!"))

  }

 }

 if (verbose)
  cat("\n")

 #########################################################################

 ### if requested, apply transformation function

 if (is.function(transf)) {
  if (is.null(targs)) {
   out <- lapply(out, transf)
  } else {
   out <- lapply(out, transf, targs)
  }
  out <- lapply(out, function(x) x[,-2,drop=FALSE])
  transf <- TRUE
 }

 ### make sure order of intervals is always increasing

 #tmp <- .psort(pi.lb, pi.ub)
 #pi.lb <- tmp[,1]
 #pi.ub <- tmp[,2]

 #########################################################################

 if (is.null(out)) {
  return()
 } else {
  return(out)
 }

}
