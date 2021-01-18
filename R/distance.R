#' compute beta-diversity from raw OTU data
#'
#' @param OTU_raw: The raw OTU data
#' @param meth: The string indicates the method of beta-div computation
#' @return The distance matrix of beta-diversity
#' @import vegan
#' @export
beta.dev <- function(OTU_raw, meth){
  Yi = as.matrix(OTU_raw)
  dat = Yi / rowSums(Yi)
  
  ### 2. compute beta-div for the OTU relative abundance ###
  d = vegdist(
    dat,
    method = meth,
    binary = FALSE,
    diag = FALSE,
    upper = FALSE,
    na.rm = FALSE
  )
  d.mat = dist2mat(d, 128)
  return(d.mat)
}

#' compute beta-diversity from raw OTU data
#'
#' @param df_cts: The dataframe with only continuous variables
#' @param meth: The string indicates the method of pairwise distance computation
#' @return A list with pairwise distance for each continuous variable
#' @import vegan
#' @export
pairwise.dist <- function(df_cts, meth){
  cl = ncol(df_cts)
  l = list()
  for (i in 1:cl) {
    cont.covariate = df_cts[, i]
    dat.cont.cov = as.matrix(na.omit(cont.covariate))
    ### compute euclidean distance g() for the cov ###
    gfun.cov = vegdist(
      dat.cont.cov,
      method = meth,
      binary = FALSE,
      diag = FALSE,
      upper = FALSE,
      na.rm = FALSE
    )
    gfun.cov.mat = dist2mat(gfun.cov, 128)
    l[[i]] <- gfun.cov.mat
  }
  return(l)
}