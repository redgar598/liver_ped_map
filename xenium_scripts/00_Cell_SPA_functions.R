

readBIDCell_adapted <- function(data_dir,
                                meta_idx = NULL,
                                tiff_path = NULL) {
  
  cell_outputs <- read.csv(data_dir)
  
  meta_idx <- c(1:8, (ncol(cell_outputs)-2):ncol(cell_outputs))
  
  data <- cell_outputs[, -meta_idx]
  meta <- cell_outputs[, meta_idx]
  
  duplicated_cell_id <- duplicated(cell_outputs$cell_id)
  
  if (sum(duplicated_cell_id) > 0) {
    warning(paste("There are", sum(duplicated_cell_id),
                  "cells with duplicated cell id"))
    data <- aggregate.Matrix(data, cell_outputs$cell_id)
    rownames(data) <- paste("Cell", rownames(data), sep = "_")
    meta <- meta[!duplicated_cell_id, ]
    rownames(meta) <- paste("Cell", meta$cell_id, sep = "_")
    meta <- meta[rownames(data), ]
  } else {
    rownames(data) <- paste("Cell", cell_outputs$cell_id, sep = "_")
  }
  
  data <- as(as.matrix(data), "CsparseMatrix")
  
  spe <- SpatialExperiment(assay = list(counts = Matrix::t(data)),
                           colData = data.frame(meta))
  
  spe <- .initialise_CellSPA_list(spe)
  spe <- .add_dataset_metrics(spe)
  
  tif_output <- readTiffOutput(tiff_path)
  tif_output <- tif_output[tif_output$cell_id %in% spe$cell_id, ]
  
  
  spe@metadata$CellSegOutput <- tif_output
  
  if (!all(spe$cell_id %in% spe@metadata$CellSegOutput$cell_id)) {
    warning("There are cells that are not in the tiff file,
                    please check your tiff file")
    spe <- spe[, spe$cell_id %in% spe@metadata$CellSegOutput$cell_id]
  }
  
  spe <- .cal_baseline(spe)
  spe@metadata$CellSPA$method <- "BIDCell"
  
  spe
}



#' Subset Spatial experiment object by column
#'
#'
#' @param spe A Spatial Experiment object
#' @param col_idx The column index to keep.
#'
#' @return A SpatialExperiment object
#' @export
subset <- function(spe, col_idx) {
  spe <- spe[, col_idx]
  if (!is.null(spe@metadata$CellSegOutput)) {
    keep <- spe@metadata$CellSegOutput$cell_id %in% spe$cell_id
    spe@metadata$CellSegOutput <- spe@metadata$CellSegOutput[keep, ]
  }
  
  return(spe)
}

.add_metrics <- function(spe, metrics_name, type) {
  spe@metadata$CellSPA$metrics[[type]] <- sort(unique(c(spe@metadata$CellSPA$metrics[[type]],
                                                        metrics_name)))
  return(spe)
}


.add_dataset_metrics <- function(spe) {
  spe@metadata$CellSPA$dataset_level_metrics <- list(num_cells = ncol(spe),
                                                     num_genes = nrow(spe))
  spe <- .add_metrics(spe, c("num_cells", "num_genes"), "dataset_level")
  return(spe)
}


.initialise_CellSPA_list <- function(spe) {
  spe@metadata$CellSPA <- list()
  spe@metadata$CellSPA$metrics <- vector("list", 4)
  names(spe@metadata$CellSPA$metrics) <- c("dataset_level",
                                           "gene_level",
                                           "cell_level",
                                           "bin_level")
  return(spe)
}


# From Matrix.utils::aggregate.Matrix() function, which has been archived on CRAN

aggregate.Matrix <- function(x, groupings = NULL, form = NULL, fun = 'sum', ...) {
  if (!methods::is(x,'Matrix')) {
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  }
  if (fun=='count') {
    x <- x != 0
  }
  groupings2 <- groupings
  if (!methods::is(groupings2,'data.frame')) {
    groupings2 <- as(groupings2,'data.frame')
  }
  
  groupings2 <- data.frame(lapply(groupings2, as.factor))
  groupings2 <- data.frame(interaction(groupings2,sep = '_'))
  colnames(groupings2) <- 'A'
  if(is.null(form)) {
    form <- stats::as.formula('~0+.')
  }
  
  form <- stats::as.formula(form)
  mapping <- dMcast(groupings2,form)
  colnames(mapping) <- substring(colnames(mapping),2)
  result <- Matrix::t(mapping) %*% x
  
  
  attr(result,'crosswalk') <- grr::extract(groupings,match(rownames(result),groupings2$A))
  return(result)
}

# From Matrix.utils::dMcast() function, which has been archived on CRAN

dMcast <- function(data, formula, fun.aggregate='sum',
                   value.var=NULL, as.factors=FALSE,
                   factor.nas=TRUE, drop.unused.levels=TRUE) {
  values <- 1
  if( !is.null(value.var)) {
    values<-data[,value.var]
  }
  
  alltms <- stats::terms(formula, data = data)
  response <- rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm <- attr(alltms,"term.labels")
  interactionsIndex <- grep(':',tm)
  interactions <- tm[interactionsIndex]
  simple <- setdiff(tm,interactions)
  i2 <- strsplit(interactions,':')
  newterms <- unlist(lapply(i2, function (x) {
    paste("paste(", paste(x,collapse=','),",","sep='_'",")")
  }))
  newterms <- c(simple,newterms)
  newformula <- stats::as.formula(paste('~0+',paste(newterms,collapse='+')))
  allvars <- all.vars(alltms)
  data <- data[,c(allvars),drop=FALSE]
  if (as.factors) {
    data <- data.frame(lapply(data,as.factor))
  }
  
  characters <- unlist(lapply(data,is.character))
  data[,characters] <- lapply(data[,characters,drop = FALSE], as.factor)
  factors <- unlist(lapply(data, is.factor))
  #Prevents errors with 1 or fewer distinct levels
  data[, factors] <- lapply(data[, factors, drop = FALSE], function (x) {
    if(factor.nas) {
      if(any(is.na(x))) {
        levels(x) <- c(levels(x),'NA')
        x[is.na(x)] <- 'NA'
      }
    }
    
    if (drop.unused.levels) {
      if(nlevels(x) != length(stats::na.omit(unique(x)))) {
        x <- factor(as.character(x))
      }
    }
    
    
    y <- stats::contrasts(x, contrasts = FALSE, sparse = TRUE)
    attr(x,'contrasts') <- y
    return(x)
  })
  #Allows NAs to pass
  attr(data,'na.action') <- stats::na.pass
  result <- Matrix::sparse.model.matrix(newformula, data, drop.unused.levels = FALSE, row.names=FALSE)
  brokenNames <- grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames] <- lapply(colnames(result)[brokenNames], function (x) {
    x <- gsub('paste(',replacement = '', x = x, fixed = TRUE)
    x <- gsub(pattern = ', ',replacement='_', x = x, fixed = TRUE)
    x <- gsub(pattern = '_sep = \"_\")', replacement = '', x=x, fixed = TRUE)
    return(x)
  })
  
  result <- result*values
  # if (isTRUE(response>0)) {
  #     responses = all.vars(terms(as.formula(paste(response, '~0'))))
  #     result <- aggregate.Matrix(result, data[,responses,drop=FALSE],
  #                                fun = fun.aggregate)
  # }
  return(result)
}


#' Calculate the baseline metrics of segmented cell from cell segmentation output
#'
#'
#' @param spe A Spatial Experiment object
#' @param metrics A character vector indicating the metrics to be calculated
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @importFrom SingleCellExperiment counts
#' @return A SpatialExperiment object
#' @export

calBaselineCellMetrics <- function(spe,
                                   metrics = c("total_transciprts", "total_genes",
                                               "total_cells", "meanExprsPct_cells"),
                                   rerun = TRUE,
                                   verbose = TRUE) {
  
  metrics <- match.arg(metrics, c("total_transciprts", "total_genes",
                                  "total_cells", "meanExprsPct_cells"), several.ok = TRUE)
  if (!rerun) {
    metrics <- setdiff(metrics, unlist(spe@metadata$CellSPA$metrics))
  }
  
  if (verbose) {
    print(paste("Metrics to run: ", paste(metrics, collapse = ", ")))
  }
  
  spe <- .cal_baseline(spe, metrics = metrics)
  return(spe)
}

#' Calculate the baseline metrics of cell shapes from cell segmentation output
#'
#'
#' @param spe A Spatial Experiment object
#' @param metrics A character vector indicating the metrics to be calculated
#' @param use_BPPARAM a BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @import BiocParallel

#' @return A SpatialExperiment object
#' @export

calBaselineCellShapeMetrics <- function(spe,
                                        metrics = c("cell_area", "elongation",
                                                    "compactness", "eccentricity",
                                                    "sphericity", "solidity",
                                                    "convexity", "circularity"),
                                        use_BPPARAM = BiocParallel::SerialParam(),
                                        rerun = TRUE,
                                        verbose = TRUE) {
  
  metrics <- match.arg(metrics, c("cell_area", "elongation",
                                  "compactness", "eccentricity",
                                  "sphericity", "solidity",
                                  "convexity", "circularity"), several.ok = TRUE)
  
  if (!rerun) {
    metrics <- setdiff(metrics, unlist(spe@metadata$CellSPA$metrics))
    if (verbose) {
      print(paste("Metrics to run: ", paste(metrics, collapse = ", ")))
    }
    
  }
  
  if (verbose) {
    BiocParallel::bpprogressbar(use_BPPARAM) <- TRUE
  }
  
  spe <- .cal_baseline_tiff(spe, metrics = metrics, use_BPPARAM = use_BPPARAM,
                            verbose = verbose)
  
  return(spe)
}




#' Calculate the baseline metrics of cell from cell segmentation output
#'
#'
#' @param spe A Spatial Experiment object
#' @param metrics A character vector indicating the metrics to be calculated
#' @param use_BPPARAM a BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @import BiocParallel

#' @return A SpatialExperiment object
#' @export

calBaselineAllMetrics <- function(spe,
                                  metrics = c("total_transciprts", "total_genes",
                                              "total_cells", "meanExprsPct_cells",
                                              "cell_area", "elongation",
                                              "compactness", "eccentricity",
                                              "sphericity", "solidity",
                                              "convexity", "circularity",
                                              "density"),
                                  use_BPPARAM = BiocParallel::SerialParam(),
                                  rerun = TRUE,
                                  verbose = TRUE) {
  
  non_cell_shape_metrics <- intersect(metrics,
                                      .CellSPARenvir[["non_cell_shape_baseline_metrics"]])
  
  
  if (length(non_cell_shape_metrics) != 0) {
    spe <- calBaselineCellMetrics(spe,
                                  metrics = non_cell_shape_metrics,
                                  rerun = rerun,
                                  verbose = verbose)
  }
  
  cell_shape_metrics <- intersect(metrics,
                                  .CellSPARenvir[["cell_shape_baseline_metrics"]])
  
  if (length(non_cell_shape_metrics) != 0 & !is.null(spe@metadata$CellSegOutput)) {
    spe <- calBaselineCellShapeMetrics(spe,
                                       metrics = cell_shape_metrics,
                                       use_BPPARAM = use_BPPARAM,
                                       rerun = rerun,
                                       verbose = verbose)
  }
  
  if ("density" %in% metrics & !is.null(spe@metadata$CellSegOutput)) {
    spe <- density(spe)
  }
  
  return(spe)
}



#' @importFrom SummarizedExperiment colData rowData
#' @importFrom Matrix rowMeans rowSums colSums

.cal_baseline <- function(spe, metrics = c("total_transciprts", "total_genes",
                                           "total_cells", "meanExprsPct_cells")) {
  if (!"CellSPA" %in% names(spe@metadata)) {
    spe@metadata$CellSPA <- list()
  }
  
  
  if ("total_transciprts" %in% metrics) {
    
    total_transciprts <- Matrix::colSums(counts(spe))
    colData(spe)$total_transciprts <- total_transciprts
    spe@metadata$CellSPA$metrics <- c(spe@metadata$CellSPA$metrics,
                                      "total_transciprts")
    spe <- .add_metrics(spe, "total_transciprts", "cell_level")
  }
  
  if ("total_genes" %in% metrics) {
    total_genes <- Matrix::colSums(counts(spe) != 0)
    colData(spe)$total_genes <- total_genes
    spe <- .add_metrics(spe, "total_genes", "cell_level")
  }
  
  if ("total_cells" %in% metrics) {
    total_cells <- Matrix::rowSums(counts(spe) != 0)
    SummarizedExperiment::rowData(spe)$total_cells <- total_cells
    spe <- .add_metrics(spe, "total_cells", "gene_level")
  }
  
  if ("meanExprsPct_cells" %in% metrics) {
    meanExprsPct_cells <- Matrix::rowMeans(counts(spe) != 0)
    SummarizedExperiment::rowData(spe)$meanExprsPct_cells <- meanExprsPct_cells
    spe <- .add_metrics(spe, "meanExprsPct_cells", "gene_level")
  }
  
  
  
  return(spe)
}



#' Generating polygon
#'
#'
#' @param spe A Spatial Experiment object
#' @param use_BPPARAM a BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @import BiocParallel

#' @return A SpatialExperiment object
#'
#' @export

generatePolygon <- function(spe,
                            use_BPPARAM = BiocParallel::SerialParam(),
                            rerun = TRUE,
                            verbose = TRUE) {
  
  
  tiff_res <- spe@metadata$CellSegOutput
  tiff_res_list <- split(tiff_res[, c(1:2)], tiff_res$cell_id)
  
  if (verbose) {
    BiocParallel::bpprogressbar(use_BPPARAM) <- TRUE
  }
  
  
  if (is.null(spe@metadata$CellSPA$poly_objects) | rerun) {
    spe@metadata$CellSPA$poly_objects <- BiocParallel::bplapply(tiff_res_list, function(x) {
      get_ashape_chull_poly(x)
    }, BPPARAM = use_BPPARAM)
  }
  
  
  return(spe)
}




#' @import BiocParallel

.cal_baseline_tiff <- function(spe,
                               metrics = c("cell_area", "elongation",
                                           "compactness", "eccentricity",
                                           "sphericity", "solidity",
                                           "convexity", "circularity"),
                               use_BPPARAM = SerialParam(),
                               verbose = verbose) {
  tiff_res <- spe@metadata$CellSegOutput
  tiff_res_list <- split(tiff_res[, c(1:2)], tiff_res$cell_id)
  # filter <- unlist(lapply(tiff_res_list, nrow)) <= 3
  # tiff_res_list <- tiff_res_list[!filter]
  
  if ("cell_area" %in% metrics) {
    if (!is.null(tiff_res)) {
      cell_area <- table(tiff_res$cell_id)
      cell_area <- as.numeric(cell_area[as.character(spe$cell_id)])
      colData(spe)$cell_area <- cell_area
      spe <- .add_metrics(spe, "cell_area", "cell_level")
    }
  }
  
  
  if (is.null(spe@metadata$CellSPA$poly_objects)) {
    if (verbose) {
      print("Generating ashape and chull...")
    }
    spe@metadata$CellSPA$poly_objects <- BiocParallel::bplapply(tiff_res_list, function(x) {
      get_ashape_chull_poly(x)
    }, BPPARAM = use_BPPARAM)
    
  }
  
  poly_objects <- spe@metadata$CellSPA$poly_objects
  
  for (i in setdiff(metrics, "cell_area")) {
    spe <- .run_cell_shape_metrics(spe,
                                   poly_objects,
                                   method = i,
                                   use_BPPARAM = use_BPPARAM,
                                   verbose = verbose)
  }
  
  return(spe)
}


density <- function(spe) {
  colData(spe)$density <- spe$total_transciprts/spe$cell_area
  spe <- .add_metrics(spe, "density", "cell_level")
  return(spe)
}

.CellSPARenvir <- new.env(parent = emptyenv())
.CellSPARenvir[["cell_level_baseline_metrics"]] <- c("total_transciprts",
                                                     "total_genes",
                                                     "cell_area", "elongation",
                                                     "compactness", "eccentricity",
                                                     "sphericity", "solidity",
                                                     "convexity", "circularity")


.CellSPARenvir[["cell_shape_baseline_metrics"]] <- c("cell_area", "elongation",
                                                     "compactness", "eccentricity",
                                                     "sphericity", "solidity",
                                                     "convexity", "circularity")


.CellSPARenvir[["non_cell_shape_baseline_metrics"]] <- c("total_transciprts",
                                                         "total_genes",
                                                         "total_cells",
                                                         "meanExprsPct_cells")

