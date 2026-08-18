#' Build expanded coordinate system accounting for all alleles
#'
#' @param vcf A vcfR object containing haploid calls, including for the reference haplotype
#' @param start_pos Start position in reference coordinates (NULL for beginning)
#' @param end_pos End position in reference coordinates (NULL for end)
#' @param var_handling How to handle expanded coordinate system when there are >1 distinct alleles at a single locus.
#' Options are: 'collapsed' (default) = all SV alleles begin at the same recorded REF position *and* expanded
#' coordinate position (e.g., all deletion and insertion alleles will be vertically stacked, regardless of allelic
#' type/identity); 'expanded' = each SV allele (i.e., VCF record) receives its own expanded coordinate start position
#' (still not sure how best to implement this, so not currently available)
#' @return A data frame with mapping between reference and expanded coordinates
#' @import vcfR
#' @export
build_coordsystem <- function(vcf, start_pos = NULL, end_pos = NULL, var_handling = 'collapsed') {

  if(missing(vcf)){
    stop("Must specify 'vcf': a vcfR object")
  }
  if(missing(start_pos)){
    stop("Must specify 'start_pos' for region to be plotted")
  }
  if(missing(end_pos)){
    stop("Must specify 'end_pos' for region to be plotted")
  }

  ## Extract positions and alleles
  fix <- getFIX(vcf)
  fix <- fix[order(as.numeric(fix[,2])),]
  positions <- as.numeric(fix[, 'POS'])
  refs <- fix[, 'REF']
  alts <- fix[, 'ALT']

  ## Determine region boundaries
  if (is.null(start_pos)) {
    start_pos <- min(positions)
  }
  if (is.null(end_pos)) {
    end_pos <- max(positions)
  }

  ## Filter variants to region
  in_region <- positions >= start_pos & positions <= end_pos
  positions <- positions[in_region]
  refs <- refs[in_region]
  alts <- alts[in_region]
  ## make sure all lengths are equal
  if(length(positions) != length(refs) | length(positions) != length(alts)){
    stop("Check VCF formatting - issue with POS / REF / ALT counts")
  }
  ## make sure no multiallelic variants
  if(any(grepl(',', alts)) | any(grepl(',', refs))){
    stop("Check VCF formatting - must not contain multiallelic records")
  }

  ## Initialize coordinate mapping
  coord_map <- data.frame(
    ref_pos = integer(),
    expanded_pos = integer(),
    variant_type = character(),
    max_length = integer(),
    stringsAsFactors = FALSE
  )

  expanded_pos <- 1   ## expanded coordinate system starts at 1 (exp 1 == ref start_pos)
  last_ref_pos <- start_pos - 1

  for(i in unique(positions)){ ## instead of processing variants sequentially, process within each REF position sequentially
    ref_pos <- i
    var.idx <- which(positions == i) ## get indices (== VCF row nums) for variants at this reference position
    ref_alleles <- refs[var.idx]
    alt_alleles <- alts[var.idx]

    # probably don't need this check? ope, we do.
    if(substr(ref_alleles[1], 1, 1) != substr(alt_alleles[1], 1, 1)){
      warning(paste0('Check VCF formatting - first nucleotide not identical between REF / ALT allele at position ',i))
    }

    ## Add intervening reference positions if moving more than 1 base position (i.e., not processing
    ## another allele at the same locus just processed or at the subsequent position)
    if(ref_pos > last_ref_pos + 1){
      tmp.ref_pos <- (last_ref_pos + 1):(ref_pos - 1) ## the intervening positions in ref coords
      tmp.expanded_pos <- seq(from = expanded_pos, by = 1, length.out = length(tmp.ref_pos)) ## " expanded coords
      tmp.variant_type <- rep('ref', length(tmp.ref_pos))
      tmp.max_length <- rep(1, length(tmp.ref_pos))
      interven.map <- data.frame(ref_pos = tmp.ref_pos, ## coord_map rows for intervening coords
                                 expanded_pos = tmp.expanded_pos,
                                 variant_type = tmp.variant_type,
                                 max_length = tmp.max_length)
      coord_map <- rbind(coord_map, interven.map)
      expanded_pos <- max(coord_map$expanded_pos) + 1
    }

    ## Don't need to record variant types at each ref position, actually,
    ## if just collapsing all variants at each locus. may need to deal with this
    ## differently if keeping distinct alleles expanded. So keep for now for
    ## matrix structure consistency between the 2 approaches.

    ##### Collapsed variant handling approach #####
    ## Summarize variant types at locus (will separate them for expanded variant approach)
    ## e.g., here: if(var_handling == 'collapsed'){

    if(var_handling == 'collapsed'){
    ## Determine maximum ALT allele length to get amount expanded coordinates need to accommodate
    ## If only variant = deletion, coordinates will not expand, will just advance by 1 like
    ## the ref coord system
      max_alt_len <- max(nchar(alt_alleles))
      max_ref_len <- max(nchar(ref_alleles))

      if(max_alt_len > 1 & max_ref_len > 1){
        var_type <- 'indel'
      } else if(max_alt_len > 1 & max_ref_len == 1){
        var_type <- 'ins'
      } else if(max_alt_len == 1 & max_ref_len > 1){
        var_type <- 'del'
      } else if(max_alt_len == 1 & max_alt_len == max_ref_len){
        var_type <- 'snp'
      } else if(max_alt_len > 1 & max_alt_len == max_ref_len &
                length(alt_alleles) == 1 & length(ref_alleles) == 1){
        var_type <- 'mnp'
      }

      if(var_type == 'del'){
        ## deletions don't expand the coord system
        coord_map <- rbind(coord_map, data.frame(
          ref_pos = ref_pos,
          expanded_pos = expanded_pos,
          variant_type = var_type,
          max_length = max_ref_len
        ))
        expanded_pos <- expanded_pos + 1
      } else if(max_ref_len > max_alt_len){
        ## indels where REF > ALT don't expand the coord system
        coord_map <- rbind(coord_map, data.frame(
          ref_pos = ref_pos,
          expanded_pos = expanded_pos,
          variant_type = var_type,
          max_length = max_ref_len
        ))
        expanded_pos <- expanded_pos + 1
      } else{
        ## everything else (insertions, SNPs, MNPs, indels where ALT > REF) does by max allele size
        coord_map <- rbind(coord_map, data.frame(
          ref_pos = ref_pos,
          expanded_pos = expanded_pos,
          variant_type = var_type,
          max_length = max_alt_len
        ))
        expanded_pos <- expanded_pos + max_alt_len
      }

      last_ref_pos <- ref_pos
    }

    ##### Expanded variant approach #####
    ##### keep these bones for expanded var approach v v v
    # ref_len <- nchar(ref_allele)
    # alt_len <- nchar(alt_allele)
    # max_len <- max(ref_len, alt_len)
    # if(alt_len > ref_len){
    #   var_type <- 'ins'
    # } else if(ref_len > alt_len){
    #   var_type <- 'del'
    # } else if(ref_len == alt_len & ref_len == 1){
    #   var_type <- 'snp'
    # } else if(ref_len > alt_len & ref_len > 1){
    #   var_type <- 'mnp'
    # }

    #   ## Add mapping for this variant position
    #   ## Use max_length to accommodate longest allele
    #   coord_map <- rbind(coord_map, data.frame(
    #     ref_pos = ref_pos,
    #     expanded_pos = expanded_pos,
    #     variant_type = var_type,
    #     max_length = max_len
    #   ))
    #
    #   if(var_type == 'del'){ ## deletions don't expand the coord system
    #     expanded_pos <- expanded_pos + 1
    #   } else{
    #     expanded_pos <- expanded_pos + max_len ## everything else does by max allele size
    #   }
    #
    #   last_ref_pos <- ref_pos
    # }

  }
  ##### Applies to both approaches #####
  ## Add any remaining reference positions to end_pos
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
  # plot(coord_map$ref_pos, coord_map$expanded_pos) ## QC check plot
}
