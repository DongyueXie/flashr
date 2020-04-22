# @title Update precision parameter
#
# @description Updates the estimated precision to increase the value of
#   the objective function.
#
# @inheritParams flash
#
# @param f A flash object.
#
# @return An updated flash object.
#
flash_update_precision = function(data,
                                  f,
                                  var_type) {
  R2 = flash_get_R2(data, f)
  f$tau = compute_precision(R2, data, var_type)
  return(f)
}


compute_precision = function(R2, data, var_type) {
  if (data$anyNA) {
    R2[data$missing] = NA
  }

  if (var_type == "by_column") {
    tau = mle_precision_by_column(R2)
  }
  else if (var_type == "by_row") {
    tau = t(mle_precision_by_column(t(R2)))
  }
  else if (var_type == "constant") {
    tau = mle_precision_constant(R2)
  }
  else if (var_type == "zero") {
    tau = 1 / data$S^2
  }else if(var_type == 'by_column+known'){
    tau = 1 / mle_var_by_column_known(R2,data$S)
  }

  if (is.matrix(tau) && data$anyNA) {
    tau[data$missing] = 0
  }

  return(tau)
}

# @title estimate sigma^2 x_i\sim N(0,\sigma^2+s_i^2)
#
# @param y2 squared x_i
# @param varx0 known variance s_i^2
# @return a scalar
#
varx_est = function(y2,varx0,bound=c(-1e5,1e5)){
  normaleqn=function(varx,y2,varx0){
    return(sum(y2/(varx+varx0)^2)-sum(1/(varx+varx0)))
  }
  tt = try(uniroot(normaleqn,bound,y2=y2,varx0=varx0),silent = TRUE)
  if(class(tt)=='try-error'){
    return(0)
  }else{
    return(pmax(tt$root,0))
  }
}

# R2: squared residual matrix
# S: known standard error matrix
# return variance matrix
mle_var_by_column_known = function(R2,S){
  for(i in 1:ncol(R2)){
    S[,i] = (S[,i]^2+varx_est(R2[,i],S[,i]^2))
  }
  S
}
# @title MLE for precision (separate parameter for each column)
#
# @param R2 An n by p matrix of squared residuals (with NAs for missing).
#
# @return An n by p matrix of precisions (separate value for each column).
#
mle_precision_by_column = function (R2) {
  sigma2 = colMeans(R2, na.rm = TRUE)  # a p vector

  # If a value of tau becomes numerically negative, set it to a
  # small positive number.
  tau = pmax(1/sigma2, .Machine$double.eps)
  return(outer(rep(1, nrow(R2)), tau))
}


# @title MLE for precision (single value)
#
# @param R2 An n by p matrix of squared residuals (with NAs for missing).
#
# @return An n by p matrix of precisions (a single value).
#
mle_precision_constant = function(R2) {
  sigma2 = mean(R2, na.rm = TRUE)  # a scalar

  tau = pmax(1/sigma2, .Machine$double.eps)
  return(tau)
}
