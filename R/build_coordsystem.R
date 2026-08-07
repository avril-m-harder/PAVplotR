#' Build expanded coordinate system accounting for all alleles
#'
#' @param vcf_fn A VCF file containing haploid calls, including for the reference haplotype
#' @param start_pos Start position in reference coordinates (NULL for beginning)
#' @param end_pos End position in reference coordinates (NULL for end)
#' @return A data frame with mapping between reference and expanded coordinates
#' @import vcfR
#' @export
build_coordsystem <- function(vcf_fn, start_pos = NULL, end_pos = NULL) {

  if(missing(vcf_fn)){
    stop("Must specify 'vcf_fn': VCF path + file name")
  }
  if(missing(start_pos)){
    stop("Must specify 'start_pos'")
  }
  if(missing(end_pos)){
    stop("Must specify 'end_pos'")
  }

  vcf <- read.vcfR(vcf_fn)

  # Extract positions and alleles
  fix <- getFIX(vcf)
  fix <- fix[order(as.numeric(fix[,2])),]
  positions <- as.numeric(fix[, 'POS'])
  refs <- fix[, "REF"]
  alts <- fix[, "ALT"]

  # Determine region boundaries
  if (is.null(start_pos)) {
    start_pos <- min(positions)
  }
  if (is.null(end_pos)) {
    end_pos <- max(positions)
  }

  # Filter variants to region
  in_region <- positions >= start_pos & positions <= end_pos
  positions <- positions[in_region]
  refs <- refs[in_region]
  alts <- alts[in_region]

  # Initialize coordinate mapping
  coord_map <- data.frame(
    ref_pos = integer(),
    expanded_pos = integer(),
    variant_type = character(),
    max_length = integer(),
    stringsAsFactors = FALSE
  )

  expanded_pos <- 1
  last_ref_pos <- start_pos - 1

  for (i in seq_along(positions)) {
    ref_pos <- positions[i]
    ref_allele <- refs[i]
    alt_alleles <- strsplit(alts[i], ",")[[1]]

    # Add intervening reference positions
    if (ref_pos > last_ref_pos + 1) {
      tmp.ref_pos <- (last_ref_pos + 1):(ref_pos - 1)
      tmp.expanded_pos <- seq(from = expanded_pos, by = 1, length.out = length(tmp.ref_pos))
      tmp.variant_type <- rep('ref', length(tmp.ref_pos))
      tmp.max_length <- rep(1, length(tmp.ref_pos))
      interven.map <- data.frame(ref_pos = tmp.ref_pos,
                                 expanded_pos = tmp.expanded_pos,
                                 variant_type = tmp.variant_type,
                                 max_length = tmp.max_length)
      coord_map <- rbind(coord_map, interven.map)
      expanded_pos <- max(coord_map$expanded_pos) + 1
    }

    # Calculate max allele length at this position
    all_alleles <- c(ref_allele, alt_alleles)
    allele_lengths <- nchar(all_alleles)
    max_len <- max(allele_lengths)

    # Determine variant type
    ref_len <- nchar(ref_allele)
    has_insertion <- any(allele_lengths > ref_len)
    has_deletion <- any(allele_lengths < ref_len)

    if (has_insertion && has_deletion) {
      var_type <- "indel"
    } else if (has_insertion) {
      var_type <- "insertion"
    } else if (has_deletion) {
      var_type <- "deletion"
    } else {
      var_type <- "snp"
    }

    # Add mapping for this variant position
    # Use max_length to accommodate longest allele
    coord_map <- rbind(coord_map, data.frame(
      ref_pos = ref_pos,
      expanded_pos = expanded_pos,
      variant_type = var_type,
      max_length = max_len
    ))

    if(var_type == 'deletion'){ ## edited - deletions shouldn't expand the coord system
      expanded_pos <- expanded_pos + 1
    } else{
      expanded_pos <- expanded_pos + max_len
    }

    last_ref_pos <- ref_pos
  }

  # Add any remaining reference positions to end_pos
  if (end_pos > last_ref_pos) {
    tmp.ref_pos <- (last_ref_pos + 1):end_pos
    tmp.expanded_pos <- seq(from = expanded_pos, by = 1, length.out = length(tmp.ref_pos))
    tmp.variant_type <- rep('ref', length(tmp.ref_pos))
    tmp.max_length <- rep(1, length(tmp.ref_pos))
    interven.map <- data.frame(ref_pos = tmp.ref_pos,
                               expanded_pos = tmp.expanded_pos,
                               variant_type = tmp.variant_type,
                               max_length = tmp.max_length)
    coord_map <- rbind(coord_map, interven.map)
    expanded_pos <- max(coord_map$expanded_pos) + 1
  }
  return(coord_map)
}
