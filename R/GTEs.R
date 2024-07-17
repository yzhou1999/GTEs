#' Compute the group technical effects.
#'
#' @param X Input data matrix.
#' @param meta Input metadata (data.frame).
#' @param g_factor Group variable (s).
#' @param b_factor Batch variable (s).
#' @param do.scale Whether to perform scaling.
#'
#' @export
#'
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @useDynLib GTEs
Run.GroupTechEffects = function(X, meta, g_factor, b_factor, do.scale = F) {
  meta[, g_factor] = as.character(meta[, g_factor])
  if (do.scale) X = scale_data(X, do.center = F)

  message("Compute the group technical effects!")
  s1 = Sys.time()
  G = group_onehot(meta, g_factor)
  BG = group_onehot(meta, c(g_factor, b_factor))
  res = gene_GroupTechEffects(G = G, BG = BG, X = as.matrix(X) )

  colnames(res$GroupTechEffects) = sort(unique(meta[, g_factor]))
  rownames(res$GroupTechEffects) = rownames(X)
  scale_factor = length(which(colSums(res$GroupTechEffects) > 0))
  res$GroupTechEffects = res$GroupTechEffects / scale_factor

  res$OverallTechEffects = as.vector(res$OverallTechEffects)
  names(res$OverallTechEffects) = rownames(X)
  res$OverallTechEffects = res$OverallTechEffects / scale_factor
  s2 = Sys.time()
  print(s2-s1)
  message("Done!")
  return(res)
}


#' Select HTGs under a group variable.
#'
#' @param GTEs GTEs object.
#' @param bins Bins.
#' @param gte.ratio Ratio of selected HTGs to overall GTEs.
#'
#' @importFrom dplyr `%>%`
#'
#' @export
Select.HTGs <- function(GTEs, bins = 0.1, gte.ratio = 0.95) {
  htgs_list <- lapply(1:ncol(GTEs$GroupTechEffects), function(i) select_htgs(GTEs$GroupTechEffects[, i], bins, gte.ratio))
  names(htgs_list) <- colnames(GTEs$GroupTechEffects)
  each_dfs = lapply(1:length(htgs_list),
                    function(i) data.frame(gte = unname(GTEs$GroupTechEffects[, i][htgs_list[[i]] ]),
                                           htg = htgs_list[[i]],
                                           data = names(htgs_list)[i]) )
  all_df = do.call(rbind, each_dfs)
  mm = table(all_df[, "htg"]) %>% as.data.frame()
  colnames(mm) = c("htg", "Freq")
  mm$htg = as.character(mm$htg)
  mm$gte = as.vector(tapply(all_df$gte, all_df$htg, sum)[as.character(mm$htg)])
  mm = mm[order(-mm$Freq, -mm$gte), ]
  all_gtes = GTEs$OverallTechEffects[mm$htg]
  other_gtes = sort(GTEs$OverallTechEffects[!names(GTEs$OverallTechEffects) %in% mm$htg], decreasing = T)
  all_gtes = c(all_gtes, other_gtes)
  all_htgs = select_htgs(all_gtes, bins, gte.ratio, is.sort = F)
  message(paste0("Find ", length(all_htgs), " HTGs, GTEs proportion: ", sum(all_gtes[all_htgs])/sum(all_gtes)))
  all_htgs
}


#' Select HTGs using GTEs vector.
#'
#' @param gtes Named GTEs vector.
#' @param bins Bins.
#' @param gte.ratio Ratio of selected HTGs to overall GTEs.
#' @param is.sort Whether to sort genes by GTEs from largest to smallest.
#'
select_htgs <- function(gtes, bins = 0.1, gte.ratio = 0.95, is.sort = T) {
  n <- 1/bins
  n_genes <- length(gtes)
  split_ns <- split(1:n_genes, ceiling(1:n_genes / ceiling(n_genes/n)))
  quantile_ns <- sapply(split_ns, max)
  GTEs <- gtes
  if (is.sort) GTEs <- sort(gtes, decreasing = T)
  ratio_GTE <- cumsum(GTEs)[quantile_ns] / sum(GTEs)
  htg_ratio <- which(ratio_GTE >= gte.ratio)[1] * bins
  htgs <- names(GTEs[GTEs > stats::quantile(GTEs, 1-htg_ratio)])
  htgs
}


#' Compute one-hot matrix for given data frame and variable (s)
#'
#' @param x Input data frame.
#' @param ivar Variable (s) for one-hot computation.
#'
#' @importFrom stats model.matrix
group_onehot <- function(x, ivar) {
  if (length(unique(x[, ivar])) == 1) {
    matrix(1, nrow = length(x[, ivar]), ncol = 1)
  } else {
    x <- data.frame(ivar = x[, ivar])
    x <- Reduce(paste0, x)
    model.matrix(~ 0 + x)
  }
}


#' Scale data matrix
#'
#' @param data.x Input data matrix.
#' @param do.center Whether center the row values. (default TRUE)
#' @param do.scale Whether scale the row values. (default TRUE)
#' @param row.means The provided row means to center. (default NULL)
#' @param row.sds The provided row standard deviations to scale. (default NULL)
#'
#' @import Matrix
scale_data <- function(data.x,
                       do.center = T,
                       do.scale = T,
                       row.means = NULL,
                       row.sds = NULL) {
  if (do.center) {
    if (is.null(row.means)) {
      data_mean <- Matrix::rowMeans(data.x)
    } else {
      data_mean <- row.means
    }
    data.x <- data.x - sapply(1:ncol(data.x), function(i) data_mean)
  }
  if (do.scale) {
    if (is.null(row.sds)) {
      data_stddev <- matrixStats::rowSds(as.matrix(data.x))
    } else {
      data_stddev <- row.sds
    }
    index <- which(data_stddev > 0)
    data.x[index, ] <- data.x[index, ] / sapply(1:ncol(data.x), function(i) data_stddev[index])
  }

  data.x
}

