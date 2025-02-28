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
#' @examples
#' # X is a normalized expression matrix with rows as features and columns as cells.
#'
#' # meta is a data.frame with columns containing metadata such as cell type, batch, etc.
#'
#' data_file <- system.file("extdata", "example_data.rds", package = "GTEs")
#' example_data <- readRDS(data_file)
#' meta_file <- system.file("extdata", "example_meta.rds", package = "GTEs")
#' example_meta <- readRDS(meta_file)
#' GTE_ct <- Run.GroupTechEffects(example_data, example_meta,
#'                                g_factor = "CellType",
#'                                b_factor = "Batch")
#' @return A list containing the overall GTE ($OverallTechEffects) and the GTE ($GroupTechEffects) of each subgroup under the group variable.
#' @useDynLib GTEs
Run.GroupTechEffects = function(X, meta, g_factor, b_factor, do.scale = FALSE) {
  meta[, g_factor] = as.character(meta[, g_factor])
  if (do.scale) X = scale_data(X, do.center = FALSE)

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


#' Select highly batch-sensitive genes (HBGs) under a group variable.
#'
#' @param GTE GTE result.
#' @param bins Bins.
#' @param gte.ratio Ratio of selected HBGs to the total GTE.
#'
#' @export
#'
#' @importFrom dplyr `%>%`
#' @examples
#' # GTE is the result of Run.GroupTechEffects function.
#' data_file <- system.file("extdata", "GTE_ct.rds", package = "GTEs")
#' GTE_ct <- readRDS(data_file)
#' HBGs <- Select.HBGs(GTE_ct)
#' @return Identified HBGs.
#' @useDynLib GTEs
Select.HBGs <- function(GTE, bins = 0.1, gte.ratio = 0.95) {
  hbgs_list <- lapply(1:ncol(GTE$GroupTechEffects), function(i) select_hbgs(GTE$GroupTechEffects[, i], bins, gte.ratio))
  names(hbgs_list) <- colnames(GTE$GroupTechEffects)
  each_dfs = lapply(1:length(hbgs_list),
                    function(i) data.frame(gte = unname(GTE$GroupTechEffects[, i][hbgs_list[[i]] ]),
                                           hbg = hbgs_list[[i]],
                                           data = names(hbgs_list)[i]) )
  all_df = do.call(rbind, each_dfs)
  mm = table(all_df[, "hbg"]) %>% as.data.frame()
  colnames(mm) = c("hbg", "Freq")
  mm$hbg = as.character(mm$hbg)
  mm$gte = as.vector(tapply(all_df$gte, all_df$hbg, sum)[as.character(mm$hbg)])
  mm = mm[order(-mm$Freq, -mm$gte), ]
  all_gte = GTE$OverallTechEffects[mm$hbg]
  other_gte = sort(GTE$OverallTechEffects[!names(GTE$OverallTechEffects) %in% mm$hbg], decreasing = TRUE)
  all_gte = c(all_gte, other_gte)
  all_hbgs = select_hbgs(all_gte, bins, gte.ratio, is.sort = FALSE)
  message(paste0("Find ", length(all_hbgs), " HBGs, GTE proportion: ", sum(all_gte[all_hbgs])/sum(all_gte)))
  all_hbgs
}


#' Select HBGs using GTE vector.
#'
#' @param gte Named GTE vector.
#' @param bins Bins.
#' @param gte.ratio Ratio of selected HBGs to overall GTE.
#' @param is.sort Whether to sort genes by GTE from largest to smallest.
#'
select_hbgs <- function(gte, bins = 0.1, gte.ratio = 0.95, is.sort = TRUE) {
  n <- 1/bins
  n_genes <- length(gte)
  split_ns <- split(1:n_genes, ceiling(1:n_genes / ceiling(n_genes/n)))
  quantile_ns <- sapply(split_ns, max)
  GTE <- gte
  if (is.sort) GTE <- sort(gte, decreasing = TRUE)
  ratio_GTE <- cumsum(GTE)[quantile_ns] / sum(GTE)
  hbg_ratio <- which(ratio_GTE >= gte.ratio)[1] * bins
  hbgs <- names(GTE[GTE > stats::quantile(GTE, 1-hbg_ratio)])
  hbgs
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
                       do.center = TRUE,
                       do.scale = TRUE,
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
