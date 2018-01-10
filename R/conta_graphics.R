
#' Plot log likelihoods and mark the top score
#'
#' Each contamination fraction (cf) value that was tested will be plotted
#' and the top scoring cf highlighted with a red cross. This plot enables
#' us to see the distribution of log likelihoods.
#'
#' @param cf_range range of contamination fraction
#' @param grid_lr initial grid that was searched
#' @param opt_val optimized result
#' @param save_dir location to save the results
#' @param sample sample name
#'
#' @return none
#'
#' @export
plot_max_likelihood <- function(cf_range, grid_lr, opt_val, save_dir, sample) {

  out_file <- file.path(save_dir, paste(sample, "likelihood.png", sep = "."))
  png(filename = out_file, width = 1280, height = 720)
  cf_range <- c(cf_range, opt_val$maximum)
  grid_lr <- c(grid_lr, opt_val$objective)[order(cf_range)]
  cf_range <- cf_range[order(cf_range)]
  plot(cf_range, grid_lr, log = "x", las = 2, type = "b",
       ylab = "average log-likelihood ratio",
       xlab = "contamination level",
       main = paste(sample, "cf_est =", round(opt_val$maximum, 5)),
       cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.25)
  abline(v = opt_val$maximum, lwd = 2, col = "red", lty = 2)
  abline(h = opt_val$objective, lwd = 2, col = "red", lty = 2)
  msg.trap <- capture.output(suppressMessages(dev.off()))
}

#' Calculate a set of likelihood metrics across SNPs
#'
#' @param save_dir folder location to write files and plots
#' @param sample name of the sample
#' @param dat data.table with likelihood information per locus
#' @param per_chr data.table with per chromosome metrics
#' @param ext_chr_table extension for per chr output table
#' @param ext_loh_table extension for loh output table
#' @param ext_loh_plot extension for loh plot
#' @param min_maf minimum allele frequency used to filter SNPs
#'
#' @return none
#'
#' @export
plot_lr <- function(save_dir, sample, dat, per_chr,
                    ext_chr_table = "per_chr.tsv",
                    ext_loh_table = "per_bin.tsv",
                    ext_loh_plot = "bin.lr.png",
                    min_maf = 0.25) {

  # Write per_chr results to file
  per_chr_file <- file.path(save_dir, paste(sample, ext_chr_table, sep = "."))
  write.table(format(per_chr, digits = 4), file = paste(per_chr_file),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Plot sorted log-likelihoods to see whether there are one or two platoes
  plot_vfn_cp(dat, save_dir, sample, min_maf)

  # Generate plots per bin
  plot_lr_per_bin(dat, save_dir, sample,
                  ext_table = ext_loh_table,
                  ext_plot = ext_loh_plot)
}


#' Plot sorted log likelihoods
#'
#' @param dat data.table with likelihood information per locus
#' @param save_dir folder to write out the results
#' @param sample sample name to put in filename and plot title
#' @param min_maf minimum allele frequency used to filter SNPs
#' @return none
#'
#' @export
plot_vfn_cp <- function(dat, save_dir, sample, min_maf) {
  png(file.path(save_dir, paste(sample, "vfn.cp.png", sep = ".")),
      width = 1280, height = 720)
  p <- ggplot(dat[ gt != "0/1" & vfn != 0, ], aes(y = vfn, x = cp))
  p <- p + geom_point(shape = ".") + geom_density2d(aes(y = vfn, x = cp))
  p <- p + scale_x_log10( breaks = c(0, 0.001, 0.01, 0.1, 0.4, 1),
                          limits = c(min_maf ^ 2, 1))
  p <- p + scale_y_log10( breaks = c(0, 0.001, 0.01, 0.1, 0.4, 1),
                          limits = c(0.00001, 1))
  p <- p + xlab("Contamination Probability")
  p <- p + ylab("Contaminant Allele Frequency")
  p <- p + scale_fill_gradient("Count", low = "white", high = "blue")
  p <- p + theme(axis.text = element_text(size = 14),
                 axis.text.x = element_text(size = 14),
                 axis.title = element_text(size = 14),
                 legend.text = element_text( size = 14))
  print(p)
  msg.trap <- capture.output(suppressWarnings(dev.off()))
}

#' Plot per bin likelihood across chromosomes
#'
#' The goal is to find localized inconcistencies
#' @param dat data.table with likelihood information per locus
#' @param save_dir folder to write out the results
#' @param sample sample name to put in filename and plot title
#' @param ext_table extension for the table output
#' @param ext_plot extension for the figure output
#' @return none
#'
#' @export
plot_lr_per_bin <- function(dat, save_dir, sample,
                            ext_table = "per_bin.tsv",
                            ext_plot = "bin.lr.png") {

  per_bin <- dat[gt != "0/1", .(numSnps = .N,
                                depth = mean(depth, na.rm = TRUE),
                                avg_lr = mean(lr, na.rm = TRUE),
                                lr = sum(lr, na.rm = TRUE),
                                vfn = mean(vfn, na.rm = TRUE),
                                cp = mean(cp, na.rm = TRUE),
                                vr = mean(vr, na.rm = TRUE)),
                 by = .(chrom, chunk)]

  per_bin_file <- file.path(save_dir, paste(sample, ext_table, sep = "."))
  write.table(format(per_bin, digits = 4), file = paste(per_bin_file),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  png(file.path(save_dir, paste(sample, ext_plot, sep = ".")),
      width = 1280, height = 720)

  p <- ggplot(per_bin, aes(chunk, avg_lr, colour = chrom))
  p <- p + geom_bar(stat = "identity") + facet_wrap(~chrom, ncol = 6)
  p <- p + ylab("Average log likelihood ratio")
  p <- p + theme(axis.text = element_text(size = 14),
                 axis.text.x = element_blank(),
                 axis.title = element_text(size = 14),
                 strip.text = element_text(size = 14),
                 legend.position = "none")
  print(p)
  msg.trap <- capture.output(suppressMessages(dev.off()))
}

#' Plot depth by chromosome
#'
#' The goal is to find depth ups and downs
#' @param dat data.table with likelihood information per locus
#' @param save_dir folder to write out the results
#' @param sample sample name to put in filename and plot title
#' @param min_depth minimum depth to visualize
#' @param ext_plot extension for the figure output
#' @return none
#'
#' @export
plot_depth_by_chr <- function(dat, save_dir, sample, min_depth,
                              ext_plot = "depth.png") {

  dat <- dat[depth >= min_depth, ]
  png(file.path(save_dir, paste(sample, ext_plot, sep = ".")),
      width = 1280, height = 720)
  p <- ggplot(dat, aes(pos, depth, colour = chrom))
  p <- p + geom_line(linetype = "dashed") + facet_wrap(~chrom, nrow = 6)
  p <- p + ylab("Depth")
  p <- p + theme(axis.text = element_text(size = 14),
                 axis.text.x = element_blank(),
                 axis.title = element_text(size = 14),
                 strip.text = element_text(size = 14),
                 legend.position = "none")
  print(p)
  msg.trap <- capture.output(suppressMessages(dev.off()))
}

#' Plot minor allele ratio sorted by chromosomes
#'
#' @param dat data.table with snp allele ratio information per locus
#' @param save_dir folder to write out the results
#' @param sample sample name to put in filename and plot title
#' @param ext_plot extension for the figure output
#' @param maxp maximum number of points to plot
#' @param seed random seed
#' @return none
#'
#' @export
plot_minor_ratio <- function(dat, save_dir, sample, ext_plot = "vr.png",
                             maxp = 20000, seed = 1) {

  # Set plot_dat as a subset of dat if it exceeds a pre-determined size
  plot_dat <- dat
  if (nrow(dat) > maxp) {
    set.seed(seed)
    plot_dat <- dat[sort(sample(1:nrow(dat), maxp)), ]
  }

  png(file.path(save_dir, paste(sample, ext_plot, sep = ".")),
      width = 1280, height = 720)
  plot_cols <- ifelse(plot_dat$gt == "0/1", "green",
                ifelse(plot_dat$gt == "1/1", "red", "blue"))
  plot(plot_dat[, minor_ratio], col = plot_cols,
       ylab = "Observed Minor Allele Frequency", xlab = "Sorted positions",
       cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.25)
  msg.trap <- capture.output(suppressMessages(dev.off()))
}
