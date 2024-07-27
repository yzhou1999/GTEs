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
#' @useDynLib GTE
Run.GroupTechEffects = function(X, meta, g_factor, b_factor, do.scale = F) {
  meta[, g_factor] = as.character(meta[, g_factor])
  if (do.scale) X = scale_data(X, do.center = F)

  message("Compute the group technical effects!")
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
  message("Done!")
  return(res)
}


#' Select HTGs under a group variable.
#'
#' @param GTE GTE result.
#' @param bins Bins.
#' @param gte.ratio Ratio of selected HTGs to overall GTE.
#'
#' @importFrom dplyr `%>%`
#'
#' @export
Select.HTGs <- function(GTE, bins = 0.1, gte.ratio = 0.95) {
  htgs_list <- lapply(1:ncol(GTE$GroupTechEffects), function(i) select_htgs(GTE$GroupTechEffects[, i], bins, gte.ratio))
  names(htgs_list) <- colnames(GTE$GroupTechEffects)
  each_dfs = lapply(1:length(htgs_list),
                    function(i) data.frame(gte = unname(GTE$GroupTechEffects[, i][htgs_list[[i]] ]),
                                           htg = htgs_list[[i]],
                                           data = names(htgs_list)[i]) )
  all_df = do.call(rbind, each_dfs)
  mm = table(all_df[, "htg"]) %>% as.data.frame()
  colnames(mm) = c("htg", "Freq")
  mm$htg = as.character(mm$htg)
  mm$gte = as.vector(tapply(all_df$gte, all_df$htg, sum)[as.character(mm$htg)])
  mm = mm[order(-mm$Freq, -mm$gte), ]
  all_gte = GTE$OverallTechEffects[mm$htg]
  other_gte = sort(GTE$OverallTechEffects[!names(GTE$OverallTechEffects) %in% mm$htg], decreasing = T)
  all_gte = c(all_gte, other_gte)
  all_htgs = select_htgs(all_gte, bins, gte.ratio, is.sort = F)
  message(paste0("Find ", length(all_htgs), " HTGs, GTE proportion: ", sum(all_gte[all_htgs])/sum(all_gte)))
  all_htgs
}


#' Select HTGs using GTE vector.
#'
#' @param gte Named GTE vector.
#' @param bins Bins.
#' @param gte.ratio Ratio of selected HTGs to overall GTE.
#' @param is.sort Whether to sort genes by GTE from largest to smallest.
#'
select_htgs <- function(gte, bins = 0.1, gte.ratio = 0.95, is.sort = T) {
  n <- 1/bins
  n_genes <- length(gte)
  split_ns <- split(1:n_genes, ceiling(1:n_genes / ceiling(n_genes/n)))
  quantile_ns <- sapply(split_ns, max)
  GTE <- gte
  if (is.sort) GTE <- sort(gte, decreasing = T)
  ratio_GTE <- cumsum(GTE)[quantile_ns] / sum(GTE)
  htg_ratio <- which(ratio_GTE >= gte.ratio)[1] * bins
  htgs <- names(GTE[GTE > stats::quantile(GTE, 1-htg_ratio)])
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

