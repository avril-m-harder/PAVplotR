#' Plot presence-absence matrix with reference coordinate labels
#'
#' @param presence_matrix Matrix of presence/absence proportions, calculated by calculate_bins.R
#' @param bin_info Data frame with bin coordinate information, written by calculate_bins.R
#' @param ins_output_prefix Output plot file prefix. Default is to substitute '-PAV_INS.pdf' (and/or '-PAV_INS.tiff') for '.vcf.gz' in input VCF file name.
#' @param output_fmt File format for output plot. Options are 'pdf' (default), 'tiff', or 'both'
#' @param roi Name of region to be plotted
#' @param ref_hap Name of reference haplotype. Should be the name of the haplotype against which variants are described (i.e., that alignments were made against to build the input VCF) and should also be included as a sample in the input VCF.
#' @param chrom Name of chromosome containing region to be plotted
#' @param region_start First position of region in ref_hap coordinate space to be plotted (1-based)
#' @param region_end Last position of region in ref_hap coordinate space to be plotted
#' @param bin_size Bin size (in bases) for summarizing PAV
#' @param ins_thresh Threshold (in bases) for the number of additional bases required in a bin, relative to the reference haplotype, to flag the bin as containing ≥1 insertion. Default = 1
#' @param color_low Color for absence (bin presence = 0)
#' @param color_high Color for presence (bin presence = 1)
#' @param ins_color_low Color for absence (insertion bin presence = 0)
#' @param ins_color_high Color for presence (insertion bin presence = 1)
#' @param width Plot width (inches)
#' @param height Plot height (inches)
#' @param gene_bounds Optional BED file of gene regions to be overlaid on PAV plot. Follows typical tab-delimited BED format with 4 columns: (i) chromosome, (ii) 0-based start coordinate, (iii) end coordinate, and (iv) gene name
#' @param gene_color If supplying gene regions, color of highlighting polygon
#' @param hap_order Method for determining vertical order of haplotypes in plot. Options are: 'refdist' (default) = haplotypes are ordered by distance to the reference haplotype, with the reference haplotype appearing at the top of the plot and haplotype divergence increases as y decreases; 'clust' = haplotypes are clustered by using dist() and hclust(), order of clusters is arbitrary
#' @import dplyr
#' @import ggnewscale
#' @import ggplot2
#' @import ggtext
#' @import reshape2
#' @importFrom stats hclust dist
#' @export
plot_pav_hiliteInsertions <- function(presence_matrix, bin_info,
                                      ins_output_prefix, output_fmt = 'pdf', roi = NULL,
                                      ref_hap = NULL, chrom = NULL, region_start = NULL, region_end = NULL,
                                      bin_size = 100, ins_thresh = 1,
                                      color_low = 'white', color_high = '#0F85A0FF',
                                      ins_color_low = 'white', ins_color_high = '#DD4124FF',
                                      width = 14, height = 8,
                                      gene_bounds = NULL, gene_color = '#EDD746FF',
                                      hap_order = 'refdist') {
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
  if(hap_order == 'refdist'){
    ## cluster rows of matrix by similarity to reference sample reorder matrix to match that clustering
    row.dists <- dist(presence_matrix)
    mat.dists <- as.matrix(row.dists)
    ref.dists <- mat.dists[ref_hap,]
    ord.ind <- order(ref.dists)
    presence_matrix <- presence_matrix[ord.ind,]
  } else if(hap_order == 'clust'){
    ## cluster rows of matrix by similarity and reorder matrix to match that clustering
    row.dists <- dist(presence_matrix)
    row.clust <- hclust(row.dists, method = 'average')
    ord.ind <- row.clust$order
    presence_matrix <- presence_matrix[ord.ind,]
  }

  # Convert to long format for ggplot
  df <- melt(presence_matrix)
  colnames(df) <- c("Sample", "Bin", "Presence")

  # Extract bin number for proper ordering
  df$BinNum <- as.numeric(gsub("Bin_", "", df$Bin))

  # Merge with bin info to get reference coordinates
  df <- merge(df, bin_info, by.x = "BinNum", by.y = "bin_num")

  # Determine reasonable number of x-axis breaks
  n_breaks <- min(10, nrow(bin_info)) ## minimum of 10 or the number of bins
  break_indices <- round(seq(1, nrow(bin_info), length.out = n_breaks)) ## number of bins / 10 (or min)

  # Create custom breaks and labels with round numbers
  region_span <- region_end - region_start

  # Determine appropriate rounding interval for ~10 total labels
  log_span <- log10(region_span)
  interval_magnitude <- floor(log_span - 1)
  base_interval <- 10^interval_magnitude

  # Choose a "nice" interval (1, 2, 5, or 10 times the base)
  nice_intervals <- c(1, 2, 5, 10) * base_interval
  target_n_labels <- 10
  n_labels_each <- region_span / nice_intervals
  interval <- nice_intervals[which.min(abs(n_labels_each - target_n_labels))]

  # Generate round number positions
  interior_positions <- seq(
    from = ceiling(region_start / interval) * interval,
    to = floor(region_end / interval) * interval,
    by = interval
  )

  # Combine with start and end, ensure uniqueness and sort
  ref_positions <- unique(sort(c(region_start, interior_positions, region_end)))

  # Remove zeros if present
  ref_positions <- ref_positions[ref_positions != 0]

  # Map reference positions to expanded coordinates
  x_breaks <- numeric(length(ref_positions))
  x_labels <- character(length(ref_positions))

  for (i in seq_along(ref_positions)) {
    ref_pos <- ref_positions[i]

    # Find expanded position by looking up in coord_map
    map_idx <- which(coord_map$ref_pos == ref_pos)

    if (length(map_idx) > 0) {
      # Exact match found
      x_breaks[i] <- coord_map$expanded_pos[map_idx[1]]
    } else {
      # print(paste0('possible issue at ',i))
      # Interpolate between nearest positions
      lower_idx <- max(which(coord_map$ref_pos < ref_pos))
      upper_idx <- min(which(coord_map$ref_pos > ref_pos))

      if (length(lower_idx) > 0 && length(upper_idx) > 0) {
        # Linear interpolation
        ref_lower <- coord_map$ref_pos[lower_idx]
        ref_upper <- coord_map$ref_pos[upper_idx]
        exp_lower <- coord_map$expanded_pos[lower_idx]
        exp_upper <- coord_map$expanded_pos[upper_idx]

        frac <- (ref_pos - ref_lower) / (ref_upper - ref_lower)
        x_breaks[i] <- exp_lower + frac * (exp_upper - exp_lower)
      } else if (length(lower_idx) > 0) {
        # Use lower bound
        x_breaks[i] <- coord_map$expanded_pos[lower_idx]
      } else {
        # Use upper bound
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

  # Transform data to scale each column relative to reference sample
  # df <- df %>%
  #   group_by(expanded_start) %>%
  #   mutate(
  #     ref_value = Presence[Sample == ref][1],
  #     # Scale so that: 0 -> 0, ref_value -> 0.5, 1 -> 1
  #     Presence_scaled = case_when(
  #       Presence == 0 ~ 0,  # Keep zeros as zero
  #       is.na(ref_value) | ref_value == 0 ~ Presence,  # If ref is 0, use original
  #       Presence <= ref_value ~ Presence * (0.5 / ref_value),  # Scale 0-ref to 0-0.5
  #       Presence > ref_value ~ 0.5 + (Presence - ref_value) * (0.5 / (1 - ref_value))  # Scale ref-1 to 0.5-1
  #     )
  #   ) %>%
  #   ungroup()

  if(exists('gene_bounds')){
    n1 <- length(unique(df$Sample))+1
    gene_bounds$gene.start <- gene_bounds$gene.start+1 ## convert 0-based BED to 1-based coord system used here
    filt_gene_bounds <- gene_bounds[which(gene_bounds$gene.start >= min(bin_info$ref_start[bin_info$ref_start != 0]) &
                                            gene_bounds$gene.end <= max(bin_info$ref_end[bin_info$ref_end != 0]) ),]
    if(nrow(filt_gene_bounds) < nrow(gene_bounds)){ print("Some genes not within processed ROI bounds")}
    filt_gene_bounds$y1 <- -1
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
  if(exists('gene_bounds.poly')){
    ## newer -- works, but meh
    # ins_p <- ggplot(df, aes(x = expanded_start, y = Sample, fill = Presence_scaled)) +
    #   geom_tile(color = NA) +
    #   # scale_fill_gradient2(low = color_low, mid = color_high, high = ins_color_high,
    #   #                     name = "Proportion\nPresent",
    #   #                     midpoint = 0.5) +
    #   scale_fill_gradientn(
    #     colors = c("white", color_low, color_high, ins_color_high),
    #     values = scales::rescale(c(0, 0.001, 0.5, 1)),  # 0=white, just above 0=low_color, 0.5=mid, 1=high
    #     name = "Proportion\nPresent",
    #     limits = c(0, 1)
    #   ) +
    #   scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0,0)) +
    #   scale_y_discrete(expand = c(0,0)) +
    #   labs(title = sprintf("%s PAV", formatC(roi)),
    #        subtitle = sprintf("%s %s: %s - %s",
    #                           formatC(ref, format = "f"),
    #                           formatC(chrom, format = "f"),
    #                           formatC(region_start, format = "f",
    #                                   digits = 0, big.mark = ","),
    #                           formatC(region_end, format = "f",
    #                                   digits = 0, big.mark = ",")),
    #        x = "Reference Coordinate Position",
    #        y = "Sample") +
    #   theme_minimal() +
    #   geom_hline(yintercept = rep(1:length(unique(df$Sample)), each = 2) - 0.5, linewidth = 0.25) +
    #   geom_vline(xintercept = exp.gene.start, linewidth = 1, color = gene_color) +
    #   geom_vline(xintercept = exp.gene.end, linewidth = 1, color = gene_color) +
    #   theme(axis.text.y = element_markdown(face = axis_faces, color = 'black'),
    #         axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         plot.title = element_markdown(hjust = 0.5, face = "bold"),
    #         plot.subtitle = element_text(hjust = 0.5),
    #         legend.position = "right",
    #         # panel.background = element_rect(fill = plot_bg),
    #         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    #         axis.ticks.x = element_line(color = 'black', linewidth = 0.25),
    #         text = element_text(color = 'black'))

    ## old == works
    ins_p <- ggplot(df, aes(x = .data$expanded_start, y = .data$Sample)) +
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
            # panel.background = element_rect(fill = plot_bg),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
            axis.ticks.x = element_line(color = 'black', linewidth = 0.25),
            text = element_text(color = 'black'))
  } else{
    ins_p <- ggplot(df, aes(x = .data$expanded_start, y = .data$Sample)) +
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
            # panel.background = element_rect(fill = plot_bg),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
            axis.ticks.x = element_line(color = 'black', linewidth = 0.25),
            text = element_text(color = 'black'))
  }

  if(output_fmt %in% c('pdf', 'both')){
    ggsave(paste0(ins_output_prefix,'.pdf'), ins_p, device = 'pdf', width = width, height = height, units = 'in')
  }
  if(output_fmt %in% c('tiff', 'both')){
    ggsave(paste0(ins_output_prefix,'.tiff'), ins_p, device = 'tiff', width = width, height = height, units = 'in', dpi = 400)
  }

  return(ins_p)
}
