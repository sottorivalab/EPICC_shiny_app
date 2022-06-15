library(dplyr)
library(GenomicRanges)
source("functions/atac.R")
source("functions/annotation.R")

data_source = "~/Dropbox (ICR)/EPICC_MENDELEY/atac"
extend_by = 10000
out_dir = "data/bed_files"

min_q_value_peaks = 0.001
min_enrichment_peaks = 4.0

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# datasets ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

final_rec_count = readRDS("data/other/final_reccurence_data.rds")
genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

atac_purity_data = readRDS("data/other/genotyping_estimates_per_sample.rds")
genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

pltA = readRDS("data/other/pltA.rds")
pltB = readRDS("data/other/pltB.rds")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(out_dir)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load insertion data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

peaks_of_interest = 
  as.character(unique(c(pltA$data$peak, pltB$data$peak)))

regions_of_interest =
  peaks_of_interest %>% 
  unique() %>% 
  as.character() %>%
  as("GRanges") %>%
  "+"(extend_by)

input_beds = 
  file.path(data_source, "reads") %>%
  list.files("EPICC", full.names = TRUE) %>% 
  magrittr::set_names(gsub("_nuc.*", "", basename(.))) 

insertion_data_annot = 
  basename(input_beds) %>% 
  annotation_from_barcode(TRUE) %>% 
  dplyr::mutate(purity = atac_purity_data[sample_barcode, "estimated_purity"])

all_peaks =  as(peaks_of_interest, "GRanges")
all_genome = as(genome@seqinfo, "GRanges")[1:22, ]
all_non_peaks = setdiff(all_genome, all_peaks)

grouping = 
  with(insertion_data_annot,{
    tissue_type = as.character(tissue_type)
    ifelse(
      region == "E", 
      "Normals", 
      paste0(patient, ".", region, " (", Hmisc::capitalize(tissue_type), ")")
    )})

bed_grouped = split(input_beds, grouping)

group_annot = NULL
for (i in seq_along(bed_grouped)) {
  
  print(names(bed_grouped)[i])
  
  d = bed_grouped[[i]] %>% 
    pbapply::pblapply(subset_bed_file) %>% 
    as("GRangesList") %>% 
    unlist()
  
  c_group_annot = 
    data.frame(
      group = names(bed_grouped)[i], 
      n_peak = sum(overlapsAny(d, all_peaks)), 
      n_bkgr = sum(overlapsAny(d, all_non_peaks)), 
      size_peak = sum(width(reduce(all_peaks))),
      sizebkgr = sum(width(reduce(all_non_peaks)))
    )
  
  group_annot = rbind(group_annot, c_group_annot)
  d_export = d[overlapsAny(d, regions_of_interest),]
  
  out_file = file.path(out_dir, paste0(names(bed_grouped)[i], ".rds"))
  dir.create(dirname(out_file))
  saveRDS(d_export, file=out_file)
}

# save other infos as rds files
saveRDS(insertion_data_annot, "data/other/insertion_data_annot.rds")
saveRDS(group_annot, "data/other/group_annot.rds")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load peak data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Loading peak sets\n")

load_peak_set = function(x, i) {
  readr::read_tsv(x, comment = "#", progress=FALSE, col_types=readr::cols()) %>% 
    dplyr::filter(`-log10(pvalue)` >= -log10(min_q_value_peaks)) %>% 
    dplyr::filter(fold_enrichment >= min_enrichment_peaks) %>%
    dplyr::mutate(start=abs_summit-250, end=abs_summit+250) %>% 
    dplyr::select(chr, start, end) %>% 
    as("GRanges") %>% 
    (function(d) d[overlapsAny(d, i)])
}

get_name = function(x) {
  gsub("-GRCh38.*", "", gsub("_nucleosome_free-region_", ".", basename(x)))
}

input_peaks = 
  file.path(data_source, "peaks") %>%
  list.files(full.names = TRUE) %>% 
  magrittr::set_names(get_name(.)) %>% 
  pbapply::pblapply(load_peak_set, regions_of_interest)

mt = match(names(input_peaks), gsub(" .*", "", group_annot$group))
mt[names(input_peaks) == "pan_patient.E"] = which(group_annot$group == "Normals")
names(input_peaks) = group_annot$group[mt]
input_peaks = input_peaks[!is.na(mt)]

saveRDS(input_peaks, "data/other/input_peaks.rds")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

final_rec_count$summary %>% 
  dplyr::filter(peak %in% peaks_of_interest) %>% 
  saveRDS("data/other/recurrence_summary.rds")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cna_data = c( 
  readRDS("./data/other/sequenza_cna_data.rds"),
  readRDS("./data/other/lowpass_cn_calls.rds") 
) %>% lapply(select, chromosome, start.pos, end.pos, CNt)

cna_table = 
  get_cnas(
    mutation_data=as(peaks_of_interest, "GRanges"), 
    cna_segments=cna_data, 
    sample_ids=names(cna_data)
  )

saveRDS(cna_table, "data/other/cna_data.rds")
