#' Plot presence-absence matrix in expanded coordinate system labeled with original reference coordinates
#'
#' @param coord_map Map linking expanded coordinate system to original reference coordinate system, written by build_coordsystem()
#' @param presence_matrix Matrix of presence/absence proportions, calculated by calculate_bins()
#' @param bin_info Data frame with bin coordinate information, written by calculate_bins()
#' @param output_prefix Output plot file prefix. If not set, will default to '[roi]_[ref_hap]_[chrom]_[region_start]_[region_end]_[bin_size]-INS_PAV'
#' @param output_fmt File format for output plot. Options are 'pdf' (default), 'tiff', or 'both'
#' @param roi Name of region to be plotted
#' @param ref_hap Name of reference haplotype. Should be the name of the haplotype against which variants are
#' described (i.e., that alignments were made against to build the input VCF) and should also be included as
#' a sample in the input VCF.
#' @param chrom Name of chromosome containing region to be plotted
#' @param region_start First position of region in ref_hap coordinate space to be plotted (1-based)
#' @param region_end Last position of region in ref_hap coordinate space to be plotted
#' @param bin_size Bin size (in bases) for summarizing PAV
#' @param color_low Color for absence (bin presence = 0)
#' @param color_high Color for presence (bin presence = 1)
#' @param ins_thresh Threshold (in bases) for the number of additional bases required in a bin, relative to
#' the reference haplotype, to flag the bin as containing ≥1 insertion. Default = 1
#' @param ins_color_low Color for absence in bin where haplotype has more sequence than reference (bin presence = 0)
#' @param ins_color_high Color for presence in bin where haplotype has more sequence than reference (bin presence = 1)
#' @param width Plot width (inches)
#' @param height Plot height (inches)
#' @param gene_bounds Optional BED file of gene regions to be overlaid on PAV plot. Follows typical
#' tab-delimited BED format with 4 columns: (i) chromosome, (ii) 0-based start coordinate, (iii) end
#' coordinate, and (iv) gene name [somewhat functional]
#' @param gene_color If supplying gene regions, color of highlighting polygon [somewhat functional]
#' @param hap_order Method for determining vertical order of haplotypes in plot. Options are: 'refdist'
#' (default) = haplotypes are ordered by distance to the reference haplotype, with the reference
#' haplotype appearing at the top of the plot and haplotype divergence increasing as y decreases; 'clust'
#'  = haplotypes are clustered by using dist() and hclust(), order of clusters is arbitrary
#' @import dplyr
#' @import ggnewscale
#' @import ggplot2
#' @import ggtext
#' @import reshape2
#' @importFrom stats hclust dist
#' @export
plot_pav_hiliteInsertions <- function(coord_map, presence_matrix, bin_info,
                                      output_prefix = NULL, output_fmt = 'pdf', roi = NULL,
                                      ref_hap = NULL, chrom = NULL, region_start = NULL, region_end = NULL,
                                      bin_size = 100, color_low = 'white', color_high = '#0F85A0FF',
                                      ins_thresh = 1,
                                      ins_color_low = 'white', ins_color_high = 'seagreen4',
                                      width = 14, height = 8,
                                      gene_bounds = NULL, gene_color = '#EDD746FF',
                                      hap_order = 'refdist') {

  if(missing(coord_map)){
    stop("Must specify 'coord_map'")
  }
  if(missing(presence_matrix)){
    stop("Must specify 'presence_matrix'")
  }
  if(missing(bin_info)){
    stop("Must specify 'bin_info'")
  }
  if(missing(roi)){
    stop("Must specify 'roi'")
  }
  if(missing(ref_hap)){
    stop("Must specify 'ref_hap'")
  }
  if(missing(chrom)){
    stop("Must specify 'chrom'")
  }
  if(missing(region_start)){
    stop("Must specify 'region_start'")
  }
  if(missing(region_end)){
    stop("Must specify 'region_end'")
  }

  if(is.null(output_prefix)){
    output_prefix <- paste0(roi,'_',ref_hap,'_',chrom,'_',region_start,'_',region_end,'_',bin_size)
  }

  if(hap_order == 'refdist'){
    ## cluster rows of matrix by similarity to reference sample reorder matrix to match that clustering
    row.dists <- dist(presence_matrix)
    mat.dists <- as.matrix(row.dists)
    ref.dists <- mat.dists[ref_hap,]
    ord.ind <- rev(order(ref.dists))
    presence_matrix <- presence_matrix[ord.ind,]
  } else if(hap_order == 'clust'){
    ## cluster rows of matrix by similarity and reorder matrix to match that clustering
    row.dists <- dist(presence_matrix)
    row.clust <- hclust(row.dists, method = 'average')
    ord.ind <- row.clust$order
    presence_matrix <- presence_matrix[ord.ind,]
  }

  ## Convert to long format for ggplot
  df <- melt(presence_matrix)
  colnames(df) <- c("Sample", "Bin", "Presence")

  ## Extract bin number for proper ordering
  df$BinNum <- as.numeric(gsub("Bin_", "", df$Bin))

  ## Merge with bin info to get reference coordinates
  df <- merge(df, bin_info, by.x = "BinNum", by.y = "bin_num")
  if(hap_order == 'refdist'){ ## set levels so that ref_hap is always plotted at the top
    df$Sample <- factor(df$Sample, levels = c(setdiff(levels(df$Sample), ref_hap), ref_hap))
  }

  ## Determine reasonable number of x-axis breaks
  n_breaks <- min(10, nrow(bin_info)) ## minimum of 10 or the number of bins
  break_indices <- round(seq(1, nrow(bin_info), length.out = n_breaks)) ## number of bins / 10 (or min)

  ## Create custom breaks and labels with round numbers
  region_span <- region_end - region_start

  ## Determine appropriate rounding interval for ~10 total labels
  log_span <- log10(region_span)
  interval_magnitude <- floor(log_span - 1)
  base_interval <- 10^interval_magnitude

  ## Choose a "nice" interval (1, 2, 5, or 10 times the base)
  nice_intervals <- c(1, 2, 5, 10) * base_interval
  target_n_labels <- 10
  n_labels_each <- region_span / nice_intervals
  interval <- nice_intervals[which.min(abs(n_labels_each - target_n_labels))]

  ## Generate round number positions
  interior_positions <- seq(
    from = ceiling(region_start / interval) * interval,
    to = floor(region_end / interval) * interval,
    by = interval
  )

  ## Combine with start and end, ensure uniqueness and sort
  ref_positions <- unique(sort(c(region_start, interior_positions, region_end)))

  ## Remove zeros if present
  ref_positions <- ref_positions[ref_positions != 0]

  ## Map reference positions to expanded coordinates
  x_breaks <- numeric(length(ref_positions))
  x_labels <- character(length(ref_positions))

  for (i in seq_along(ref_positions)) {
    ref_pos <- ref_positions[i]

    ## Find expanded position by looking up in coord_map
    map_idx <- which(coord_map$ref_pos == ref_pos)

    if (length(map_idx) > 0) {
      ## Exact match found
      x_breaks[i] <- coord_map$expanded_pos[map_idx[1]]
    } else {
      ## print(paste0('possible issue at ',i))
      ## Interpolate between nearest positions
      lower_idx <- max(which(coord_map$ref_pos < ref_pos))
      upper_idx <- min(which(coord_map$ref_pos > ref_pos))

      if (length(lower_idx) > 0 && length(upper_idx) > 0) {
        ## Linear interpolation
        ref_lower <- coord_map$ref_pos[lower_idx]
        ref_upper <- coord_map$ref_pos[upper_idx]
        exp_lower <- coord_map$expanded_pos[lower_idx]
        exp_upper <- coord_map$expanded_pos[upper_idx]

        frac <- (ref_pos - ref_lower) / (ref_upper - ref_lower)
        x_breaks[i] <- exp_lower + frac * (exp_upper - exp_lower)
      } else if (length(lower_idx) > 0) {
        ## Use lower bound
        x_breaks[i] <- coord_map$expanded_pos[lower_idx]
      } else {
        ## Use upper bound
        x_breaks[i] <- coord_map$expanded_pos[upper_idx]
      }
    }

    x_labels[i] <- formatC(ref_pos, format = "f", digits = 0, big.mark = ",")
  }

  ## categorize hap x bin combos by whether they contain insertions relative to reference
  ## that exceed a threshold, if set
  df$seq.type <- 'ref'
  prop.thresh <- ins_thresh/bin_size
  for(bin in unique(df$BinNum)){
    sub.bin <- df[df$BinNum == bin,]
    for(hap in unique(df$Sample[df$Sample != ref_hap])){
      if(sub.bin[which(sub.bin$Sample == hap), 'Presence'] >=
         (sub.bin[which(sub.bin$Sample == ref_hap), 'Presence'] + prop.thresh)){
        df[which(df$BinNum == bin & df$Sample == hap), 'seq.type'] <- 'ins'
      }
    }
  }

  ## If supplied, convert gene bounds into expanded coordinate space
  n1 <- length(unique(df$Sample))+1
  if(!is.null(gene_bounds)){
    colnames(gene_bounds) <- c('chr','gene.start','gene.end','gene.name')
    gene_bounds$gene.start <- as.numeric(gene_bounds$gene.start)
    gene_bounds$gene.end <- as.numeric(gene_bounds$gene.end)
    gene_bounds$gene.start <- gene_bounds$gene.start+1 ## convert 0-based BED to 1-based coord system used here
    gene_bounds$midpt <- (gene_bounds$gene.start + gene_bounds$gene.end)/2
    filt_gene_bounds <- gene_bounds[which(gene_bounds$gene.start >= min(bin_info$ref_start[bin_info$ref_start != 0]) &
                                            gene_bounds$gene.end <= max(bin_info$ref_end[bin_info$ref_end != 0]) ),]
    if(nrow(filt_gene_bounds) < nrow(gene_bounds)){ print("Some genes not within processed ROI bounds")}
    filt_gene_bounds$y1 <- n1-0.3
    filt_gene_bounds$y2 <- n1
    filt_gene_bounds.x <- NULL
    filt_gene_bounds.y <- NULL
    filt_gene_bounds.group <- NULL
    filt_gene_bounds.geneside <- NULL
    filt_gene_bounds$exp.midpt <- NA
    for(r in 1:nrow(filt_gene_bounds)){
      filt_gene_bounds.x <- c(filt_gene_bounds.x, filt_gene_bounds$gene.start[r], filt_gene_bounds$gene.end[r],
                              filt_gene_bounds$gene.end[r], filt_gene_bounds$gene.start[r])
      filt_gene_bounds.y <- c(filt_gene_bounds.y, filt_gene_bounds$y1[r], filt_gene_bounds$y1[r],
                              filt_gene_bounds$y2[r], filt_gene_bounds$y2[r])
      filt_gene_bounds.group <- c(filt_gene_bounds.group, filt_gene_bounds$gene.name[r],
                                  filt_gene_bounds$gene.name[r], filt_gene_bounds$gene.name[r], filt_gene_bounds$gene.name[r])
      filt_gene_bounds.geneside <- c(filt_gene_bounds.geneside, 'left','right','right','left')
      filt_gene_bounds$exp.midpt[r] <- bin_info[which(bin_info$ref_start <= filt_gene_bounds$midpt[r] &
                                                        bin_info$ref_end >= filt_gene_bounds$midpt[r]), 'expanded_start']+(bin_size/2)
    }
    filt_gene_bounds.poly <- data.frame(x = filt_gene_bounds.x, y = filt_gene_bounds.y,
                                        p = filt_gene_bounds.group, side = filt_gene_bounds.geneside)
    ## convert x coords to appx bin starts/ends in expanded coord space
    filt_gene_bounds.poly$exp.x <- NA
    for(r in 1:nrow(filt_gene_bounds.poly)){
      if(filt_gene_bounds.poly$side[r] == 'left'){
        filt_gene_bounds.poly$exp.x[r] <- bin_info[which(bin_info$ref_start <= filt_gene_bounds.poly$x[r] &
                                                           bin_info$ref_end >= filt_gene_bounds.poly$x[r]), 'expanded_start']
      } else if(filt_gene_bounds.poly$side[r] == 'right'){
        filt_gene_bounds.poly$exp.x[r] <- bin_info[which(bin_info$ref_start <= filt_gene_bounds.poly$x[r] &
                                                           bin_info$ref_end >= filt_gene_bounds.poly$x[r]), 'expanded_end']
      }

    }
  }

  axis_faces <- ifelse(levels(factor(df$Sample)) == ref_hap, "bold", "plain")

  if(exists('filt_gene_bounds.poly')){
    ins_p <- ggplot(df, aes(x = .data$expanded_start, y = .data$Sample, fill = .data$Presence)) +
      geom_polygon(data = filt_gene_bounds.poly, mapping = aes(x = .data$exp.x, y = .data$y, group = .data$p),
                   inherit.aes = FALSE,
                   fill = gene_color, alpha = 1) +
      geom_tile(data = df %>% filter(.data$seq.type == "ref"),
                aes(fill = .data$Presence),
                color = NA, width = bin_size) +
      scale_fill_continuous(name = "Presence (ref)", limits = c(0,1),
                            low = color_low, high = color_high) +
      new_scale_fill() +
      geom_tile(data = df %>% filter(.data$seq.type == "ins"),
                aes(fill = .data$Presence),
                color = NA, width = bin_size) +
      scale_fill_continuous(name = "Presence (ins)", limits = c(0,1),
                            low = ins_color_low, high = ins_color_high) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0)) +
      labs(title = sprintf("%s PAV", formatC(roi)),
           subtitle = sprintf("%s %s: %s - %s",
                              formatC(ref_hap, format = "f"),
                              formatC(chrom, format = "f"),
                              formatC(region_start, format = "f",
                                      digits = 0, big.mark = ","),
                              formatC(region_end, format = "f",
                                      digits = 0, big.mark = ",")),
           x = "Reference Coordinate Position",
           y = "Sample") +
      theme_minimal() +
      ## annotate("text", x = filt_gene_bounds$exp.midpt,
      ##          y = rep(c(n1+0.05, n1+0.1, n1+0.15), nrow(filt_gene_bounds))[1:nrow(filt_gene_bounds)],
      ##          label = filt_gene_bounds$gene.name, size = 1) +
      coord_cartesian(xlim = range(df$expanded_start), ylim = c(1,n1+0.1)) +
      geom_hline(yintercept = rep(1:length(unique(df$Sample)), each = 2) - 0.5, linewidth = 0.25) +
      ## geom_vline(xintercept = filt_gene_bounds.poly$exp.x, linewidth = 0.5, color = gene_color) +
      theme(axis.text.y = element_markdown(face = axis_faces, color = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_markdown(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right",
            legend.title = element_text(hjust = 0.5),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
            axis.ticks.x = element_line(color = 'black', linewidth = 0.25),
            text = element_text(color = 'black'))
  } else{
    ins_p <- ggplot(df, aes(x = .data$expanded_start, y = .data$Sample, fill = .data$Presence)) +
      geom_tile(data = df %>% filter(.data$seq.type == "ref"),
                aes(fill = .data$Presence),
                color = NA, width = bin_size) +
      scale_fill_continuous(name = "Presence (ref)", limits = c(0,1),
                            low = color_low, high = color_high) +
      new_scale_fill() +
      geom_tile(data = df %>% filter(.data$seq.type == "ins"),
                aes(fill = .data$Presence),
                color = NA, width = bin_size) +
      scale_fill_continuous(name = "Presence (ins)", limits = c(0,1),
                            low = ins_color_low, high = ins_color_high) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0)) +
      labs(title = sprintf("%s PAV", formatC(roi)),
           subtitle = sprintf("%s %s: %s - %s",
                              formatC(ref_hap, format = "f"),
                              formatC(chrom, format = "f"),
                              formatC(region_start, format = "f",
                                      digits = 0, big.mark = ","),
                              formatC(region_end, format = "f",
                                      digits = 0, big.mark = ",")),
           x = "Reference Coordinate Position",
           y = "Sample") +
      theme_minimal() +
      geom_hline(yintercept = rep(1:length(unique(df$Sample)), each = 2) - 0.5, linewidth = 0.25) +
      theme(axis.text.y = element_markdown(face = axis_faces, color = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_markdown(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right",
            legend.title = element_text(hjust = 0.5),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
            axis.ticks.x = element_line(color = 'black', linewidth = 0.25),
            text = element_text(color = 'black'))
  }
  if(output_fmt %in% c('pdf', 'both')){
    ggsave(paste0(output_prefix,'.pdf'), ins_p, device = 'pdf',
           width = width, height = height, units = 'in')
  }
  if(output_fmt %in% c('tiff', 'both')){
    ggsave(paste0(output_prefix,'.tiff'), ins_p, device = 'tiff',
           width = width, height = height, units = 'in', dpi = 400)
  }

  return(ins_p)
}
