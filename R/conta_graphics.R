
#' Plot log likelihoods and mark the top score
#'
#' Each contamination fraction (cf) value that was tested will be plotted
#' and the top scoring cf highlighted with a red cross. This plot enables
#' us to see the distribution of log likelihoods.
#'
#' @param cf_range numeric range of contamination fraction
#' @param grid_lr numeric initial grid that was searched
#' @param opt_val list optimized result
#' @param save_dir character location to save the results
#' @param sample character sample name
#'
#' @return none
#'
#' @export
plot_max_likelihood <- function(cf_range, grid_lr, opt_val, save_dir, sample) {

  out_file <- file.path(save_dir, paste(sample, "likelihood.png", sep = "."))
  png(filename = out_file)
  cf_range <- c(cf_range, opt_val$maximum)
  grid_lr <- c(grid_lr, opt_val$objective)[order(cf_range)]
  cf_range <- cf_range[order(cf_range)]
  plot(cf_range, grid_lr, log = "x", las = 2, type = "b",
       ylab = "average log-likelihood ratio",
       xlab = "contamination level",
       main = paste(sample, "cf_est =", round(opt_val$maximum, 5)))
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
#' @return none
#'
#' @export
plot_lr <- function(save_dir, sample, dat, per_chr) {

  # Write per_chr results to file
  per_chr_file <- file.path(save_dir, paste(sample, "per_chr.tsv", sep = "."))
  write.table(format(per_chr, digits = 4), file = paste(per_chr_file),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Plot sorted log-likelihoods to see whether there are one or two platoes
  plot_sorted_lr(dat, save_dir, sample)

  # Generate plots per bin
  plot_lr_per_bin(dat, save_dir, sample)
}


#' Plot sorted log likelihoods
#'
#' @param dat data.table with likelihood information per locus
#' @param save_dir character folder to write out the results
#' @param sample character sample name to put in filename and plot title
#' @return none
#'
#' @export
plot_sorted_lr <- function(dat, save_dir, sample) {
  png(file.path(save_dir, paste(sample, "sorted.lr.png", sep = ".")))
  plot(sort(dat[lr > 0, lr], decreasing = TRUE),
       pch = ".", main = sample, ylab = "log likelihood ratio", xlab = "SNPs")
  msg.trap <- capture.output(suppressMessages(dev.off()))
}

#' Plot per bin likelihood across chromosomes
#'
#' The goal is to find localized inconcistencies
#' @param dat data.table with likelihood information per locus
#' @param save_dir character folder to write out the results
#' @param sample character sample name to put in filename and plot title
#' @return none
#'
#' @export
plot_lr_per_bin <- function(dat, save_dir, sample) {

  portions <- 10 # partions per chr
  dat[, chunk := 0]
  for (j in dat[, unique(chrom)]) {
    bin_size <- ceiling(as.numeric(dat[chrom == j, .(.N / portions)]))
    dat[chrom == j, chunk := .(ceiling(.I / bin_size))]
  }
  dat$chrom <- factor(dat$chrom, levels = unique(dat$chrom))

  per_bin <- dat[, .(numSnps = .N,
                        depth = mean(depth, na.rm = TRUE),
                        avg_lr = mean(lr, na.rm = TRUE),
                        lr = sum(lr, na.rm = TRUE),
                        vfn = mean(vfn, na.rm = TRUE),
                        cp = mean(cp, na.rm = TRUE),
                        vr = mean(vr, na.rm = TRUE)),
                     by = .(chrom, chunk)]

  per_bin_file <- file.path(save_dir, paste(sample, "per_bin.tsv", sep = "."))
  write.table(format(per_bin, digits = 4), file = paste(per_bin_file),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  png(file.path(save_dir, paste(sample, "bin.lr.png", sep = ".")))
  p <- ggplot(per_bin, aes(chunk, avg_lr, colour = chrom))
  p <- p + geom_bar(stat = "identity") + facet_wrap(~chrom)
  p <- p + ylab("Average log likelihood ratio")
  print(p)
  msg.trap <- capture.output(suppressMessages(dev.off()))
}

#' Plot per bin het ratio ratio across chromosomes
#'
#' The goal is to find localized inconcistencies
#' @param dat data.table with het rates
#' @param save_dir character folder to write out the results
#' @param sample character sample name to put in filename and plot title
#' @return distortion metric cv of het divergence from 0.5
#'
#' @export
plot_het_dist <- function(dat, save_dir, sample) {

  portions <- 10 # partions per chr
  dat[, chunk := 0]
  for (j in dat[, unique(chrom)]) {
    bin_size <- ceiling(as.numeric(dat[chrom == j, .(.N / portions)]))
    dat[chrom == j, chunk := .(ceiling(.I / bin_size))]
  }
  dat$chrom <- factor(dat$chrom, levels = unique(dat$chrom))

  per_bin <- dat[, .(numSnps = .N,
                     depth = mean(depth),
                     het_dist = mean(abs(0.5 - major_ratio))),
                 by = .(chrom, chunk)]

  png(file.path(save_dir, paste(sample, "bin_het_rate.png", sep = ".")))
  p <- ggplot(per_bin, aes(chunk, het_dist - mean(het_dist), colour = chrom))
  p <- p + geom_bar(stat = "identity") + facet_wrap(~chrom)
  p <- p + ylab("Mean het distortion")
  plot(p)
  msg.trap <- capture.output(suppressMessages(dev.off()))

  png(file.path(save_dir, paste(sample, "bin_depth.png", sep = ".")))
  p <- ggplot(per_bin, aes(chunk, depth - mean(depth)), colour = chrom)
  p <- p + geom_bar(stat = "identity") + facet_wrap(~chrom)
  plot(p)
  msg.trap <- capture.output(suppressMessages(dev.off()))

  return(mad(per_bin[, het_dist]) / mean(per_bin[, het_dist]))
}
