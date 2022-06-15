annotation_from_barcode_epicc = function (barcodes, extract = FALSE) {

  epicc_regex =
    paste0(
      "(?P<read_barcode>",
      "(?P<lane_barcode>",
      "(?P<iteration_barcode>",
      "(?P<sample_barcode>",
      "(?P<project>EPICC)_",
      "(?P<patient>C[0-9]+)_",
      "(?P<region>[A-Z])(?P<region_number>[0-9]+)_",
      "(?P<sample_type>[BGL])(?P<sample_number>[0-9]+)_",
      "(?P<analyte>[DRCBL])(?P<analyte_number>[0-9]+))",
      "(?:_I(?P<iteration>[0-9]+))?)",
      "(?(iteration)_(?P<lane_id>L[0-9A-Za-z]+)|)?)",
      "(?(lane_id)_R(?P<read_number>[012]))?)"
    )

  numeric_columns =    # list of columns/elements that
    c("region_number", # should be converted to numeric values
      "sample_number",
      "analyte_number",
      "iteration",
      "read_number")


  # mapping of ids to annotations:
  analyte_id_names = c(
    "D" = "WGS",
    "L" = "LP-WGS",
    "C" = "ATAC-seq",
    "R" = "RNA-seq",
    "B" = "Bisulfit WGS"
  )

  type_id_names = c(
    "B" = "bulk",
    "G" = "gland",
    "L" = "interglandular",
    "Z" = "blood"
  )

  tt_order = c(
    "normal",
    "adenoma",
    "cancer"
  )

  if (!exists("msi_positiv")) { # default value for msi positive cases if variable does not exist globally
    msi_positiv = c("C536","C548","C516","C518","C552","C562")
  }

  if (!exists("msi_positiv_adenoma")) { # default value for msi positive cases if variable does not exist globally
    msi_positiv_adenoma = c("C516")
  }

  # check if the input is valid:
  if (!is.vector(barcodes)) {
    stop(paste("Argument", sQuote("barcodes"), "has to be a vector.\n"))
  } else {
    barcodes = as.character(barcodes)
  }

  if (!is.logical(extract)) {
    stop(paste("Argument", sQuote("extract"), "has to be a boolean.\n"))
  }

  # check for non matching barcodes:
  regexpr_result = regexpr(epicc_regex, barcodes, perl = TRUE)
  nerr = sum(attr(regexpr_result, "match.length") == -1, na.rm=TRUE)
  if (nerr) {
    stop(sprintf("Error: %d barcode(s) do not meet EPICC specification.\n", nerr))
  }


  # check if a valid barcode can be extracted from the input:
  barcodes_extracted = regmatches(barcodes, regexpr_result)
  n_extr = sum(barcodes_extracted != barcodes)
  if (n_extr) {
    if (extract) {
      msg = sprintf("Extracted %d barcode(s) from supplied strings.\n", n_extr)
      warning(msg)
    }
    else {
      msg = sprintf("Error: %d barcode(s) do not meet EPICC specification.\n", n_extr)
      stop(msg)
    }
    regexpr_result = regexpr(epicc_regex, barcodes_extracted,  perl=TRUE)
  }


  # get the annotation elements:
  annotation =
    regcapturedmatches(barcodes_extracted, regexpr_result) %>%
    data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate(lane_barcode=ifelse(lane_id == "", NA, lane_barcode)) %>%
    dplyr::mutate(iteration_barcode=ifelse(iteration == "", NA, iteration_barcode)) %>%
    dplyr::mutate(read_barcode=ifelse(read_number == "", NA, read_barcode)) %>%
    dplyr::mutate(tissue_barcode=gsub("_[DRCBL][0-9]+$", "", sample_barcode))

  if (sum(duplicated(barcodes)) == 0) {
    rownames(annotation) = barcodes
  }


  # insert tissue type:
  annotation$tissue_type =
    with(annotation, {
      dplyr::case_when( # some exceptions from the rule ...
        region %in% c("F") & patient == "C542" ~ "cancer",
        region %in% c("C", "D") & patient == "C516" ~ "adenoma",
        region %in% c("E", "Z", "W") ~ "normal",
        region %in% c("A", "B", "C", "D") ~ "cancer",
        region %in% c("F", "G", "H", "I") ~ "adenoma"
      ) %>% factor(tt_order, ordered=TRUE)
    })


  # insert long name for analytes:
  annotation$analyte_name =
    with(annotation, {
      analyte_id_names[as.character(analyte)] %>%
        factor(analyte_id_names, ordered=TRUE)
    })


  # insert long name for sample type:
  annotation$sample_type_name =
    with(annotation, {
      dplyr::case_when(
        region == "Z" ~ "blood",
        TRUE ~ type_id_names[as.character(sample_type)]
      ) %>% factor(c(type_id_names), ordered=TRUE)
    })


  # convert some cols to numeric:
  for (col in numeric_columns) {
    annotation[, col] = as.numeric(as.character(annotation[, col]))
  }


  # insert a label for each tumour (e.g. independed adenomas):
  group = paste0(annotation$patient, ".", annotation$tissue_type)
  wh_adenoma = grepl("adenoma", group) # add adenoma number to labels
  adenoma_regions = annotation$region[wh_adenoma]
  adenoma_regions = gsub("[CD]", "C+D", adenoma_regions)

  adenoma_region_label_list =
    split(adenoma_regions, annotation$patient[wh_adenoma]) %>%
    lapply(function(x) {
      xu = unique(x)
      if (length(xu) > 1) { l = paste0(" (", xu, ")") } else { l = "" }
      names(l) = xu
      return(l)
    }) %>% unlist()

  key_label = paste0(annotation$patient[wh_adenoma], ".", adenoma_regions)
  adenoma_labels = adenoma_region_label_list[key_label]
  group[wh_adenoma] = paste0(group[wh_adenoma], adenoma_labels)
  annotation$tumour_id = group


  # add msi status
  annotation$msi_status =
    with(annotation, {
      dplyr::case_when(
        tissue_type == "normal" ~ as.character(NA),
        tissue_type == "cancer" & patient %in% msi_positiv ~ "MSI",
        tissue_type == "adenoma" & patient %in% msi_positiv_adenoma ~ "MSI",
        TRUE ~ "MSS"
      )
    })

  return(annotation)
}

regcapturedmatches = function (x, m)
{
  if (length(x) != length(m))
    stop(gettextf("%s and %s must have the same length",
                  sQuote("x"), sQuote("m")), domain = NA)
  ili = is.list(m)
  useBytes = if (ili) {
    any(unlist(lapply(m, attr, "useBytes")))
  }
  else {
    any(attr(m, "useBytes"))
  }
  if (useBytes) {
    asc = iconv(x, "latin1", "ASCII")
    ind = is.na(asc) | (asc != x)
    if (any(ind))
      Encoding(x[ind]) = "bytes"
  }
  if (ili) {
    if (any(sapply(m, function(x) {
      is.null(attr(x, "capture.start"))
    }) == T)) {
      stop("No capture data found (did you use perl=T?)")
    }
    starts = lapply(m, function(x) {
      attr(x, "capture.start")
    })
    lengths = lapply(m, function(x) {
      attr(x, "capture.length")
    })
  }
  else {
    if (is.null(attr(m, "capture.start"))) {
      stop("No capture data found (did you use perl=T?)")
    }
    x = list(x)
    starts = list(attr(m, "capture.start"))
    lengths = list(attr(m, "capture.length"))
  }
  cleannames = function(x) {
    if (!is.null(colnames(x))) {
      colnames(x) = make.unique(make.names(colnames(x)))
    }
    x
  }
  starts = lapply(starts, cleannames)
  lengths = lapply(lengths, cleannames)
  Substring = function(x, starts, lens) {
    if (all(starts < 0)) {
      return(character())
    }
    else {
      x = t(mapply(function(x, st, ln) substring(x, st, st + ln - 1), x, data.frame(t(starts)), data.frame(t(lens)),
                   USE.NAMES = F))
      if (!is.null(colnames(starts))) {
        colnames(x) = colnames(starts)
      }
      x
    }
  }
  y = Map(function(x, sos, mls) {
    Substring(x, sos, mls)
  }, x, starts, lengths, USE.NAMES = FALSE)
  if (ili) {
    y
  }
  else {
    y[[1]]
  }
}

annotation_from_barcode = function(barcodes, extract=FALSE) {

  # quick check of input
  if (is.factor(barcodes)) barcodes = as.character(barcodes)
  checkmate::assertCharacter(barcodes, null.ok = FALSE)
  checkmate::assertFlag(extract)

  if (length(barcodes) == 0)
    return(NULL)

  # if any barcodes are duplicated only annotate unique barcodes
  if (any(duplicated(barcodes))) {
    annot = annotation_from_barcode(unique(barcodes))[barcodes,]
    rownames(annot) = NULL
    return(annot)
  }

  # try annotation with all known barcode fucntions
  annot_fun_to_test =
    list(
      annotation_from_barcode_epicc
    )

  for (annot_fun in annot_fun_to_test) {
    try({
      annot = annotation_from_barcode_epicc(barcodes, extract)
      rownames(annot) = barcodes
      return(annot)
    }, silent = TRUE)
  }

  stop("Couldn't annotate barcodes.\n")
}


get_cnas = function(mutation_data, cna_segments, sample_ids=NULL) {

  if (inherits(mutation_data, "CollapsedVCF")) {
    sample_ids = colnames(mutation_data)
    chrs = as.character(GenomeInfoDb::seqnames(mutation_data))
    pos = BiocGenerics::start(mutation_data)
    mut_gr = GenomicRanges::GRanges(chrs, IRanges::IRanges(pos, pos))
    mut_ids = rownames(mutation_data)
  } else if (inherits(mutation_data, "GRanges")) {
    mut_gr = mutation_data
    mut_ids = as.character(mut_gr)
  } else if (is.character(mutation_data)) {
    mut_gr = as(gsub("_.*", "", mutation_data), "GRanges")
    mut_ids = mutation_data
  } else {
    stop("Invalid 'mutation_data' input.\n")
  }

  if (is.null(sample_ids)) {
    stop("Missing 'sample_ids' please pass these as input.\n")
  }

  # check for missing segment data:
  sample_with_cnas = sample_ids %in% names(cna_segments)
  if (any(!sample_with_cnas)) {
    missing = paste0(sample_ids[!sample_with_cnas],collapse=", ")
    warning("Missing segment data: ", missing, ".\n", sep="")
  }

  # fill cna matrix:
  cna_status = matrix(NA, NROW(mut_gr), length(sample_ids))
  dimnames(cna_status) = list(mut_ids, sample_ids)

  for (i in which(sample_with_cnas)) { # each sample

    # get samples cna data convert to granges
    sample_cna = cna_segments[[sample_ids[i]]] # samples segment data
    wh_cols = c(chr="chromosome",start="start.pos",end="end.pos")

    sample_cna_gr =
      magrittr::set_colnames(sample_cna[, wh_cols], names(wh_cols)) %>%
      as.data.frame() %>%
      as("GRanges")

    # find overlaps add cn to cna_status matrix
    ol = findOverlaps(mut_gr, sample_cna_gr, type="within", select="first")
    cna_status[,i] = c(sample_cna$CNt[ol])

    # try to calculate average of segments for missing data
    wh_missing = which(is.na(ol))
    if (length(wh_missing) == 0) next()

    ol = findOverlaps(mut_gr[wh_missing,], sample_cna_gr, type="any")
    if (length(ol) == 0) next()

    cna_status[wh_missing,i] =
      sapply(wh_missing, function(i) {
        wh = queryHits(ol) == i
        if (sum(wh) == 0) return(NA)
        sw = pintersect(mut_gr[wh,], sample_cna_gr[subjectHits(ol[wh])])
        weighted.mean(sample_cna$CNt[[subjectHits(ol[wh])]], width(sw))
      })
  }

  return(cna_status)
}
