#' Calculate presence/absence for each sample in bins
#'
#' @param vcf A vcfR object
#' @param coord_map Coordinate mapping data frame
#' @param bin_size Size of bins in expanded coordinates
#' @param start_pos Start position in reference coordinates
#' @param end_pos End position in reference coordinates
#' @return A list with presence matrix and bin information
#' @export
calculate_bins <- function(vcf, coord_map, bin_size, start_pos, end_pos) {

  # Get genotype data
  if(any(duplicated(vcf@fix[,'ID']))){
    vcf@fix[,'ID'] <- paste0(vcf@fix[,'ID'],'.',seq(1, nrow(vcf@fix)))
  }
  gt <- extract.gt(vcf, element = "GT")

  # Get sample names
  samples <- colnames(gt)

  # Get positions and alleles
  fix <- getFIX(vcf)

  ## combine gts and site data to sort, then separate
  fix <- cbind(fix, gt)
  rm(gt)
  fix <- fix[order(as.numeric(fix[,2])),]
  positions <- as.numeric(fix[, "POS"])
  refs <- fix[, "REF"]
  alts <- fix[, "ALT"]
  gts <- fix[,c(8:ncol(fix))]
  fix <- fix[,c(1:7)]

  # Filter to region
  in_region <- positions >= start_pos & positions <= end_pos
  positions_region <- positions[in_region]
  refs_region <- refs[in_region]
  alts_region <- alts[in_region]
  gt_region <- gts[in_region, , drop = FALSE]

  # Determine expanded coordinate range
  max_expanded <- max(coord_map$expanded_pos + coord_map$max_length - 1)

  # Create bins
  bin_starts <- seq(1, max_expanded, by = bin_size)
  bin_ends <- c(bin_starts[-1] - 1, max_expanded)
  n_bins <- length(bin_starts)

  # Calculate reference position for each bin (for labeling)
  bin_ref_starts <- numeric(n_bins)
  bin_ref_ends <- numeric(n_bins)

  for (b in 1:n_bins) {
    # Find reference positions corresponding to this expanded bin
    in_bin <- coord_map$expanded_pos >= bin_starts[b] &
      coord_map$expanded_pos <= bin_ends[b]
    if (any(in_bin)) {
      bin_ref_starts[b] <- min(coord_map$ref_pos[in_bin])
      bin_ref_ends[b] <- max(coord_map$ref_pos[in_bin])
    }
  }

  # Initialize presence matrix
  presence_matrix <- matrix(0, nrow = length(samples), ncol = n_bins)
  rownames(presence_matrix) <- samples
  colnames(presence_matrix) <- paste0("Bin_", 1:n_bins)

  # For each sample, calculate presence in each bin
  for (s in seq_along(samples)) {
    s.del.pos <- NULL
    sample_name <- samples[s]

    # Track which expanded positions are present
    present_positions <- rep(FALSE, max_expanded)

    # Initially, all reference positions are present
    for (i in 1:nrow(coord_map)) {
      if (coord_map$variant_type[i] == "ref") {
        present_positions[coord_map$expanded_pos[i]] <- TRUE
      }
    }

    # Process each variant in the region (positions_region == positions of all variants in the region)
    for (v in seq_along(positions_region)){
      genotype <- gt_region[v, s]

      # Skip missing genotypes
      if (is.na(genotype) || genotype == "." || genotype == "./.") {
        next
      }

      # Parse genotype
      alleles <- as.numeric(strsplit(gsub("[|/]", " ", genotype), " ")[[1]])

      ref_pos <- positions_region[v]
      ref_allele <- refs_region[v]
      alt_alleles <- strsplit(alts_region[v], ",")[[1]]
      all_alleles <- c(ref_allele, alt_alleles)

      # Find corresponding position in coord_map
      map_idx <- which(coord_map$ref_pos == ref_pos)[1]
      if (is.na(map_idx)) next

      expanded_start <- coord_map$expanded_pos[map_idx]
      max_len <- coord_map$max_length[map_idx]

      # For each allele in genotype, mark presence
      for (allele_idx in alleles) {
        if (is.na(allele_idx)) next

        # Get the actual allele sequence (0 = ref, 1+ = alt)
        if (allele_idx == 0) {
          allele_seq <- ref_allele
          var.type <- 'ref'
        } else {
          allele_seq <- alt_alleles[allele_idx]
          if(nchar(allele_seq) > nchar(ref_allele)){
            var.type <- 'insertion'
          } else if(nchar(allele_seq) < nchar(ref_allele)){
            var.type <- 'deletion'
          } else{
            var.type <- 'snp'
          }
        }

        allele_len <- nchar(allele_seq)

        #### added accounting for deletions == absence of locus
        # Mark positions as present for this allele
        if(var.type %in% c('insertion','snp','ref')){
          for (offset in 0:(allele_len - 1)) {
            pos <- expanded_start + offset
            if (pos <= max_expanded) {
              present_positions[pos] <- TRUE
            }
          }
        } else if(var.type == 'deletion'){
          for (offset in 0:(nchar(ref_allele) - 1)) {
            pos <- expanded_start + offset
            if (pos <= max_expanded) {
              s.del.pos <- c(s.del.pos, pos) ## save it instead and overlay at the end
              # present_positions[pos] <- FALSE
            }
          }
        }
      }
    }

    ## overlay deletions (even in VCF exported from MC graph w/ vg deconstruct, '0' alleles show up
    ## where a haplotype should have a missing '.' call, i.e., at sites that should be covered by an
    ## upstream deletion in that haplotype. this approach allows those '0' calls to make it through,
    ## then be negated by documented deletion alleles in the below step)
    present_positions[s.del.pos] <- FALSE

    # Calculate proportion present in each bin
    for (b in 1:n_bins) {
      bin_positions <- bin_starts[b]:bin_ends[b]
      bin_positions <- bin_positions[bin_positions <= max_expanded]

      if (length(bin_positions) > 0) {
        presence_matrix[s, b] <- sum(present_positions[bin_positions]) / length(bin_positions)
      }
    }
  }

  # Create bin info dataframe
  bin_info <- data.frame(
    bin_num = 1:n_bins,
    expanded_start = bin_starts,
    expanded_end = bin_ends,
    ref_start = bin_ref_starts,
    ref_end = bin_ref_ends,
    ref_midpoint = (bin_ref_starts + bin_ref_ends) / 2
  )

  return(list(
    matrix = presence_matrix,
    bin_info = bin_info
  ))
}
