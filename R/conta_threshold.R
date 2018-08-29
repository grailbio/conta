# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Filter outlier samples from conta simulations
#'
#' Filters out samples with extreme avg log llr without any simulated
#' contamination, and the top X percent of samples to set thresholds more
#' aggresively.
#'
#' @param sim_data Simulation data as read in by `read_sim_data()`
#' @param extreme_level Avg LLR level above which samples will be dropped,
#'                      default of 1.5 selected as a very high score
#'                      only observed in contaimination positive controls.
#' @param filter_quantile Quantile level above which samples will be filtered
#'                        to set more aggressive thresholds. Defaults to 0.98,
#'                        throwing out the top 2 percent of samples.
#' @param mad_thresh Filter out samples that are this many MADs away from the
#'                   median, applied after the extreme level and quantile filters.
#'                   Defaults to 3 MADs
#' @return Data frame of simulation data with samples flagged as
#'         outliers removed.
#' @importFrom dplyr filter
#' @export
filter_data <- function(sim_data,
                        extreme_level = 1.5,
                        filter_quantile = 0.98,
                        mad_thresh = 3) {
  if (nrow(dplyr::filter(sim_data, conta_level == 0)) == 0) {
    warning("WARNING: Can't filter data without 0 contamination level")
  }
  # Drop any samples with super high avg llr with 0 simulated contamination,
  # largely to remove positive controls
  outliers <- dplyr::filter(sim_data, conta_level == 0,
                            avg_log_lr > extreme_level) %>%
              dplyr::select(sample) %>%
              unique()
  outlier_samples <- outliers$sample
  sim_data <- dplyr::filter(sim_data, !(sample %in% outlier_samples))

  # Drop top 2% of samples in terms avg log llr with 0 simulated contamination
  no_conta <- dplyr::filter(sim_data, conta_level == 0)
  cutoff <- quantile(no_conta$avg_log_lr, probs = c(0.98))
  top2perc <- dplyr::filter(no_conta, avg_log_lr >= cutoff) %>%
              dplyr::select(sample) %>%
              unique()
  top2perc_samples <- top2perc$sample
  sim_data <- dplyr::filter(sim_data, !(sample %in% top2perc_samples))

  # Get mean and SD of remaining uncontaminated samples, filter any samples
  # with avg_log_lr >= mean + sd_thresh * sd
  no_conta <- dplyr::filter(sim_data, conta_level == 0)
  sim_median <- median(no_conta$avg_log_lr)
  sim_mad <- mad(no_conta$avg_log_lr)

  mad_filtered <- dplyr::filter(no_conta, avg_log_lr > sim_median +
                                  mad_thresh * sim_mad) %>%
                  dplyr::select(sample) %>%
                  unique()
  mad_samples <- mad_filtered$sample
  sim_data <- dplyr::filter(sim_data, !(sample %in% mad_samples))

  message(paste("Filtering out",
                length(c(outlier_samples, top2perc_samples, mad_samples)),
                "samples."))
  message(paste("Outliers:", paste(outlier_samples, collapse = " , ")))
  message(paste("Top 2%:", paste(top2perc_samples, collapse = " , ")))
  message(paste("MAD outliers:", paste(mad_samples, collapse = " , ")))

  return(sim_data)
}

#' Make boxplot of conta level vs avg log llr
#'
#' Makes a boxplot of conta level vs avg log llr using ggplot
#'
#' @param sim_data Simulation data as read in by `read_sim_data()`
#' @return ggplot object of the plot
#' @export
plot_avg_llr_by_conta <- function(sim_data) {
  g <- ggplot(sim_data) +
       geom_boxplot(aes(x = conta_level, y = avg_log_lr)) +
       labs(x = "Simulated contamination level", y = "Avg llr")
  return(g)
}


#' Get sensitivity and specificity at specified levels
#'
#' Uses a pROC::roc object to calculate specificity, sensitivity,
#' and the associated threshold at any number of specified specificity
#' or sensitivity levels. Will find the closest specificity or sensitivity
#' that is at least as large as the specified value.
#'
#' @param roc Roc object from `pROC::roc`
#' @param specs Vector of specificity values
#' @param senss Vector of sensitivity values
#' @return Dataframe with 3 columns: specificty, sensitivity, and threshold
#' @importFrom pROC roc
#' @export
thresholds_by_sens_spec <- function(roc,
                                    specs = c(0.9, 0.95, 0.99, 1.0),
                                    senss = c(0.9, 0.95, 0.99, 1.0)) {
  results <- data.frame(specificity = numeric(),
                        sensitivity = numeric(),
                        threshold = numeric())
  # For each given specificty, find the smallest value that is greater
  # than or equal to the given value
  for (spec in specs) {
    index <- min(which(roc$specificities >= spec))
    results <- rbind(results,
                     data.frame(
                       specificity = roc$specificities[index],
                       sensitivity = roc$sensitivities[index],
                       threshold = roc$thresholds[index]))
  }

  # For each given sensitivity, find the smallest value that is greater
  # then or equal to the given value
  for (sens in senss) {
    index <- max(which(roc$sensitivities >= sens))
    results <- rbind(results,
                     data.frame(
                       specificity = roc$specificities[index],
                       sensitivity = roc$sensitivities[index],
                       threshold = roc$thresholds[index]))
  }
  return(results)
}

#' Get senstivity and specificity for each conta level
#'
#' For each simulated conta level, get the sensitivty and threshold at
#' any number of provided specificity values.
#'
#' @param sim_data Simulation data as read in by `read_sim_data()`
#' @param target_specs Vector of target specificity values
#' @return Data frame with four columns: conta_level, specificity,
#'         sensitivity, and threshold
#' @importFrom pROC roc
#' @importFrom dplyr filter
#' @export
sens_spec_by_level <- function(sim_data, target_specs = c(0.99)) {
  all_data <- data.frame(conta_level = numeric(),
                         specificity = numeric(),
                         sensitivity = numeric(),
                         threshold = numeric())

  if (nrow(dplyr::filter(sim_data, conta_level == 0)) == 0) {
    stop("ERROR: Can't get sens/spec without 0 contamination level data")
  }

  if (nrow(dplyr::filter(sim_data, conta_level == 0)) == nrow(sim_data)) {
    stop("ERROR: Can't get sens/spec with only 0 contamination level data")
  }

  # Add true label to indicate simulated contamination
  sim_data$true_label <- sim_data$conta_level != 0

  # Want to get stats for each simulated conta level other than 0
  for (this_level in unique(sim_data$conta_level)) {
    if (this_level == 0) {
      next
    }

    # Get the cases and controls for this conta level and construct the roc obj
    this_conta_data <- filter(sim_data, conta_level %in% c(0, this_level))
    this_roc <- pROC::roc(this_conta_data$true_label,
                          this_conta_data$avg_log_lr)

    this_data <- thresholds_by_sens_spec(this_roc, specs = target_specs,
                                         senss = c())
    this_data$conta_level <- this_level

    all_data <- rbind(all_data, this_data)
  }

  # Stop conta level from becoming a factor
  all_data$conta_level <- as.numeric(as.character(all_data$conta_level))
  return(all_data)
}

#' Make ROC curves by conta level
#'
#' Makes a ROC plot with one curve per conta level. Optionally add
#' vertical dashed lines at given specificty values.
#'
#' @param sim_data Simulation data as read in by `read_sim_data()`
#' @param annotate_specs Vector of specificity values to draw vertical lines at
#' @return ggplot object containing the ROC plot
#' @importFrom dplyr filter
#' @importFrom pROC roc
#' @export
make_roc_plot_by_conta_level <- function(sim_data, annotate_specs = NULL) {
  rocs <- list()

  # Add true label to indicate simulated contamination and split in to
  # contaminated and non-contaminated datasets
  sim_data$true_label <- sim_data$conta_level != 0
  non_contam <- dplyr::filter(sim_data, conta_level == 0)
  contam <- dplyr::filter(sim_data, conta_level != 0)


  # For each simulated contamination level, create the ROC object
  conta_levels <- unique(contam$conta_level)
  for (this_conta_level in conta_levels) {
    this_data <- rbind(dplyr::filter(contam, conta_level == this_conta_level),
                       non_contam)
    this_roc <- pROC::roc(this_data$true_label, this_data$avg_log_lr)
    rocs[[this_conta_level]] <- this_roc
  }

  # ggroc uses the name attribute to color the curves
  names(rocs) <- conta_levels
  gg <- ggroc(rocs) +
        theme(axis.text = element_text(size = 12, color = "Black")) +
        guides(color = guide_legend(title = "Contamination level"))
  if (!is.null(annotate_specs)) {
        gg <- gg + geom_vline(xintercept = annotate_specs,
                              color = "black",
                              linetype = "dashed")
  }

  return(gg)
}

#' Fit a linear model to predict sensitivity from contamination level
#'
#' Fits a polynomial linear regression on log(simulated contamination level)
#' to predict sensitivity. Writes a plot of the fit as well as returns the
#' model and an estimated LoD given a target sensitivity.
#'
#' @param sens_by_level Data frame of sensitivities at a fixed specificity
#'                      split out by contamination level.
#' @param target_sens Target sensitivity level for LoD estimation
#' @param fname Path to write the final plot
#' @return A list with two elements: model and lod, where `model` is the
#'         fit object returned by `lm()` and `lod` is the estimated LoD
#'         for the given target sensitivity.
#' @importFrom dplyr mutate
#' @export
fit_sens_by_conta_level <- function(sens_by_level,
                                    target_sens = 0.95) {
  sens_by_level <- dplyr::mutate(sens_by_level,
                                 log_cl.sq = log(conta_level) ^ 2,
                                 log_cl.cb = log(conta_level) ^ 3)
  sens_lm <- lm(sensitivity ~ log(conta_level) + log_cl.sq + log_cl.cb,
                data = sens_by_level)

  prediction_levels <- unique(sens_by_level$conta_level)
  predicted_sens <- predict(sens_lm, newdata = sens_by_level,
                            interval = "predict")

  sens_idx <- which(predicted_sens > target_sens)[1]
  lod <- prediction_levels[sens_idx]
  sens_lod <- predicted_sens[sens_idx]

  png(fname, height = 800, width = 1200, res = 144)
  plot(sensitivity ~ conta_level, sens_by_level, type = "b")
  lines(prediction_levels, predicted_sens[, 1], col = "red")
  lines(prediction_levels, predicted_sens[, 2], col = "red", lty = 2, lwd = 2)
  lines(prediction_levels, predicted_sens[, 3], col = "red", lty = 2, lwd = 2)
  points(lod,
         sens_lod,
         col = "red", pch = "X")
  dev.off()

  return(list(lod = lod, model = sens_lm))
}
