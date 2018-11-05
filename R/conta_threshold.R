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
#'                   Defaults to 3 MADs.
#' @return Data frame of simulation data with samples flagged as
#'         outliers removed.
#' @importFrom dplyr filter %>%
#' @export
filter_data <- function(sim_data,
                        extreme_level = 1.5,
                        filter_quantile = 0.98,
                        mad_thresh = 3) {
  `%>%` <- dplyr::`%>%`
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

  # Drop top X% of samples in terms avg log llr with 0 simulated contamination
  no_conta <- dplyr::filter(sim_data, conta_level == 0)
  cutoff <- quantile(no_conta$avg_log_lr, probs = c(filter_quantile))
  top <- dplyr::filter(no_conta, avg_log_lr >= cutoff) %>%
         dplyr::select(sample) %>%
         unique()
  top_samples <- top$sample
  sim_data <- dplyr::filter(sim_data, !(sample %in% top_samples))

  # Get median and MAD of remaining uncontaminated samples, filter any samples
  # with avg_log_lr >= median + mad_thresh * MAD
  no_conta <- dplyr::filter(sim_data, conta_level == 0)
  sim_median <- median(no_conta$avg_log_lr)
  sim_mad <- mad(no_conta$avg_log_lr)

  mad_filtered <- dplyr::filter(no_conta,
                                avg_log_lr > (sim_median + mad_thresh * sim_mad)) %>%
                  dplyr::select(sample) %>%
                  unique()
  mad_samples <- mad_filtered$sample
  sim_data <- dplyr::filter(sim_data, !(sample %in% mad_samples))

  message(paste("Filtering out",
                length(c(outlier_samples, top_samples, mad_samples)),
                "samples."))
  message(paste(sprintf("Extreme outliers (avg_log_lr > %0.2f):", extreme_level),
                paste(outlier_samples, collapse = " , ")))
  message(sprintf("Top %0.2f%%: %s\n",
                  (1 - filter_quantile) * 100,
                  paste(top_samples, collapse = " , ")))
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
  sim_data$conta_level <- factor(sim_data$conta_level)
  g <- ggplot(sim_data) +
       geom_boxplot(aes(x = conta_level, y = avg_log_lr, color = conta_level)) +
       labs(x = "Simulated contamination level", y = "Avg llr")
  return(g)
}

#' Produce a pROC::roc object from simulation data
#'
#' Returns a pROC::roc object generated from the given simulation data,
#' expects simulation data in the format returned by `read_sim_results()`.
#'
#' @param sim_data Dataframe of simulation data as produced by `read_sim_results()`
#' @return A pROC::roc object describing contamination classification from the
#'         provided simulation results
#' @importFrom pROC roc
#' @export
sim_data_to_roc <- function(sim_data) {
  # Add true label column
  sim_data$true_label <- sim_data$conta_level != 0

  # Construct roc object, explicitly specifying that we expect avg_log_lr
  # to be lower for uncontaminated samples
  return(pROC::roc(sim_data$true_label,
                   sim_data$avg_log_lr,
                   direction = '<'))
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
    this_conta_data <- dplyr::filter(sim_data, conta_level %in% c(0, this_level))
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
#' @importFrom pROC roc ggroc
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
  for (i in seq_along(conta_levels)) {
    this_conta_level <- conta_levels[i]
    this_data <- rbind(dplyr::filter(contam, conta_level == this_conta_level),
                       non_contam)
    this_roc <- pROC::roc(this_data$true_label, this_data$avg_log_lr)
    rocs[[i]] <- this_roc
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
#' @param fname Path to write the final plot
#' @param target_sens Target sensitivity level for LoD estimation
#' @return A list with two elements: model and lod, where `model` is the
#'         fit object returned by `lm()` and `lod` is the estimated LoD
#'         for the given target sensitivity.
#' @importFrom dplyr mutate
#' @export
fit_sens_by_conta_level <- function(sens_by_level,
                                    fname,
                                    target_sens = 0.95) {
  sens_lm <- lm(sensitivity ~ log(conta_level) + I(log(conta_level) ^ 2) + I(log(conta_level) ^ 3),
                data = sens_by_level)

  prediction_levels <- seq(from = min(sens_by_level$conta_level),
                           to = max(sens_by_level$conta_level),
                           length.out = 100)
  predicted_sens <- predict(sens_lm,
                            newdata = list(conta_level = prediction_levels),
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

  return(list(lod = lod, sens = sens_lod, model = sens_lm))
}


#' Run the conta thresholding workflow
#'
#' Runs the full conta thresholding workflow using a provided set of
#' simulated contamination data. Workflow consists of filtering outliers
#' from the data using three filters (see `conta::filter_data()` for details)
#' followed by calculation of sensitivities, specificities, and LoD at a
#' variety of cutoffs. Generates tables and diagnostic plots to the specified
#' output directory.
#'
#' @param sim_data Simulation data as read in by `read_sim_data()`.
#' @param outdir Path to directory to write outputs to.
#' @param extreme_level Avg LLR level above which samples will be dropped,
#'                      default of 1.5 selected as a very high score
#'                      only observed in contaimination positive controls.
#' @param filter_quantile Quantile level above which samples will be filtered
#'                        to set more aggressive thresholds. Defaults to 0.98,
#'                        throwing out the top 2 percent of samples.
#' @param mad_thresh Filter out samples that are this many MADs away from the
#'                   median, applied after the extreme level and quantile filters.
#'                   Defaults to 3 MADs.
#' @param senss Vector of sensitivity values.
#' @param specs Vector of specificity values.
#' @return None, writes out TSVs and plots to `outdir`.
#' @importFrom readr write_tsv
#' @export
conta_threshold <- function(sim_data,
                            outdir,
                            extreme_thresh = 1.5,
                            quantile_thresh = 0.98,
                            mad_thresh = 3,
                            senss = c(0.95, 0.98, 0.99, 1.0),
                            specs = c(0.95, 0.98, 0.99, 1.0)) {

  sim_data_filt <- conta::filter_data(sim_data,
                                      extreme_level = extreme_thresh,
                                      filter_quantile = quantile_thresh,
                                      mad_thresh = mad_thresh)

  filt_data_roc <- conta::sim_data_to_roc(sim_data_filt)

  sens_spec_all_data <- conta::thresholds_by_sens_spec(filt_data_roc,
                                                       specs = c(0.95, 0.98, 0.99, 1.0),
                                                       senss = c(0.95, 0.98, 0.99, 1.0))
  sens_spec_by_level <- conta::sens_spec_by_level(sim_data_filt,
                                                  target_specs = c(0.95, 0.98, 0.99, 1.0))
  # Write out table of sens and spec and thresholds
  readr::write_tsv(sens_spec_all_data,
                   file.path(outdir, "sens_spec_all_data.tsv"))
  # Write out table of sens and spec and thresholds
  readr::write_tsv(sens_spec_by_level,
                   file.path(outdir, "sens_spec_by_level.tsv"))
  lods <- tapply(seq_len(nrow(sens_spec_by_level)),
                 sens_spec_by_level$specificity,
                 function(indices) {
                   this_data <- sens_spec_by_level[indices, ]
                   this_spec <- unique(this_data$specificity)

                   # Remove levels where sens hit 100% to prevent a bad fit
                   plot_out <- file.path(outdir,
                                         sprintf("sens_fit_spec_%0.2f.png",
                                                 this_spec))
                   fit_ret <- conta::fit_sens_by_conta_level(this_data, plot_out)
                   return(data.frame(lod = fit_ret$lod,
                                     spec = this_spec,
                                     sens = fit_ret$sens))
                 })
  lods <- do.call(rbind, lods)


  # Write out table with lods
  readr::write_tsv(lods,
                   file.path(outdir, "lod_table.tsv"))

  # Plotting
  # ROC plot
  roc_plot <- conta::make_roc_plot_by_conta_level(sim_data_filt)
  png(file.path(outdir, "roc_by_conta_level.png"),
      height = 800,
      width = 1000,
      res = 144)
  print(roc_plot)
  dev.off()

  # Box plot of avg llr by conta level
  llr_boxplot <- conta::plot_avg_llr_by_conta(sim_data_filt)
  png(file.path(outdir, "boxplot_conta_level.png"),
      height = 800,
      width = 1000,
      res = 144)
  print(llr_boxplot)
  dev.off()
}
