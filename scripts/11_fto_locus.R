# 11_fto_locus.R
# Session 10: is the "FTO-high" signal in NPM1-mutated AML FTO itself, or a
# hidden transcript inside the FTO gene (IFEX9, Arza-Apalategi et al.,
# Leukemia 2023)? Beat AML only gives one count per gene, so we look at read
# coverage along the FTO locus in TCGA-LAML, re-processed by recount3.
#
#  1. Find TCGA-LAML in recount3 and check that per-base coverage (bigWig) and
#     splice-junction files exist for our samples. Record the genome build,
#     annotation and read lengths (results/fto_recount3_availability.csv).
#  2. Groups: the Session 7 TCGA cohort, NPM1-mutated split by the three-gene
#     NEK10/FTO/USP44 score (6 high, 25 low), and all 80 NPM1 wild-type as the
#     reference group.
#  3. Read coverage along FTO +/- 20 kb (GRCh38, GENCODE v26 as used by
#     recount3), scaled to 40 million 100-base reads. Group curves, each
#     sample on its own, FTO exons and IFEX9 exon 1 drawn underneath.
#  4. Split the locus into FTO exons and 1 kb intron / flank windows. Test
#     FTO-high vs FTO-low on (a) depth-adjusted coverage and (b) coverage
#     relative to the sample's own FTO exon coverage (log2 ratio).
#  5. Splice junctions inside the locus: which ones show up mostly in FTO-high
#     samples, and are they in the annotation?
#  6. IFEX9 coordinates from the Arza-Apalategi et al. supplement, checked
#     against the genome sequence with their qPCR primer.
#  7. Beat AML: annotated genes inside FTO in our counts file (Ensembl 75)
#     against the Beat AML NEK10/FTO/USP44 high/low split.
#
# Run from the repo root:  Rscript scripts/11_fto_locus.R
# Needs results/tcga_cohort.csv and results/tcga_block_scores.csv from
# 08_tcga_replication.R (Session 7), and for step 7 data/ from 01_download.R
# plus results/block_nek10_fto_scores.csv and results/fto_locus_genes_grch37.csv.
# Downloads are cached in data/recount3/ (about 230 MB: junction tables,
# annotation and the FTO slice of each bigWig, not the whole 43 GB).
# First run about 15 minutes (reading 111 bigWig slices over the network);
# later runs about 3 minutes.

suppressPackageStartupMessages({
  library(recount3)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(DESeq2)
})

set.seed(2026)
options(timeout = 3600)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
r3_dir <- "data/recount3"
dir.create(r3_dir, recursive = TRUE, showWarnings = FALSE)
# recount3 caches its own downloads with BiocFileCache; keep them in data/ too
bfc <- BiocFileCache::BiocFileCache(file.path(r3_dir, "cache"), ask = FALSE)

# ---- 1. What does recount3 have for TCGA-LAML? ---------------------------------
# recount3 (Wilks et al., Genome Biol 2021) re-aligned every TCGA RNA-seq run
# with one pipeline (STAR to GRCh38 / hg38). For each sample it stores a bigWig
# file: the number of reads covering every single base of the genome. For each
# project it stores a junction table: every place where a read was split
# across an intron (one end on one exon, the other on the next).

hp   <- available_projects(organism = "human", bfc = bfc)
proj <- hp %>% filter(project == "LAML", file_source == "tcga")
if (nrow(proj) != 1) stop("TCGA-LAML not found (or found twice) in recount3")
cat("recount3 project:", proj$project, "| source:", proj$file_source,
    "| samples:", proj$n_samples, "\n")

md <- read_metadata(file_retrieve(locate_url("LAML", proj$project_home,
                                             type = "metadata"), bfc = bfc))
md$patient <- substr(md$tcga.tcga_barcode, 1, 12)
cat("Metadata rows:", nrow(md), "| patients:", n_distinct(md$patient), "\n")

# Match to the Session 7 cohort by RNA-seq aliquot barcode (exact), not only by
# patient, so we are sure it is the same sequencing library.
cohort <- read.csv("results/tcga_cohort.csv")
blocks <- read.csv("results/tcga_block_scores.csv")
main   <- cohort %>% filter(main_cohort)
main$in_recount3 <- main$aliquot %in% md$tcga.tcga_barcode
cat("Main cohort (111) aliquots found in recount3:", sum(main$in_recount3), "of",
    nrow(main), "\n")
if (!all(main$in_recount3)) stop("Some main-cohort aliquots are missing in recount3")
cat("NPM1-mutated:", sum(main$NPM1 == "mut"), "| wild-type:", sum(main$NPM1 == "wt"),
    "| core split (Session 7):", paste(names(table(blocks$core_group)),
                                       table(blocks$core_group), collapse = ", "), "\n")

# Check that every bigWig we need is really there (HTTP HEAD: status and size,
# no download). A small file would mean an error page, not coverage. httr
# follows the server's redirect from http to https, so these are the real files.
main <- main %>% left_join(md %>% select(aliquot = tcga.tcga_barcode, external_id,
                                         rail_id), by = "aliquot")
main$bw_url <- vapply(main$external_id, function(s)
  unname(locate_url("LAML", proj$project_home, type = "bw", sample = s)), "")
head_one <- function(u) {
  h <- httr::HEAD(u)
  c(status = httr::status_code(h), bytes = as.numeric(httr::headers(h)[["content-length"]]))
}
hd <- t(vapply(main$bw_url, head_one, c(status = 0, bytes = 0)))
main$bw_status <- hd[, "status"]; main$bw_mb <- round(hd[, "bytes"] / 1e6)
cat("bigWig files: HTTP 200 for", sum(main$bw_status == 200), "of", nrow(main),
    "| size", min(main$bw_mb), "to", max(main$bw_mb), "MB each,",
    round(sum(main$bw_mb) / 1e3, 1), "GB in total\n")

jx_urls <- locate_url("LAML", proj$project_home, type = "jxn")
jx <- t(vapply(jx_urls, head_one, c(status = 0, bytes = 0)))
cat("Junction files (project level, all samples in one table):\n")
print(data.frame(file = basename(jx_urls), status = jx[, "status"],
                 MB = round(jx[, "bytes"] / 1e6, 1)), row.names = FALSE)

# Annotation recount3 uses for TCGA: GENCODE v26 (also v29) on GRCh38 / hg38.
gtf_url <- paste0(getOption("recount3_url", "http://duffel.rail.bio/recount3"),
                  "/human/annotations/exon_sums/human.exon_sums.G026.gtf.gz")
gtf_h <- head_one(gtf_url)
cat("GENCODE v26 exon annotation (recount3):", gtf_h[["status"]],
    round(gtf_h[["bytes"]] / 1e6, 1), "MB\n")

# Sequencing depth and read length for the depth adjustment later. bc_auc is
# the total coverage of the sample (sum of reads x bases over the genome).
qc <- md %>% transmute(aliquot = tcga.tcga_barcode,
                       read_length = `recount_qc.star.average_input_read_length`,
                       mapped_reads = `recount_qc.star.all_mapped_reads`,
                       auc_all = `recount_qc.bc_auc.all_reads_all_bases`,
                       junction_count = recount_qc.junction_count)
main <- main %>% left_join(qc, by = "aliquot")
cat("Average input read length:", paste(names(table(main$read_length)),
                                        table(main$read_length), sep = ": n=", collapse = "; "), "\n")
cat("Mapped reads (millions): median", round(median(main$mapped_reads) / 1e6, 1),
    "range", round(min(main$mapped_reads) / 1e6, 1), "to",
    round(max(main$mapped_reads) / 1e6, 1), "\n")

avail <- main %>%
  left_join(blocks %>% select(patient, core_group), by = "patient") %>%
  mutate(fto_group = case_when(NPM1 == "wt" ~ "NPM1_wt",
                               core_group == "high" ~ "FTO_high",
                               core_group == "low" ~ "FTO_low")) %>%
  select(patient, aliquot, NPM1, fto_group, core_group, FLT3_ITD, sex, plate,
         external_id, rail_id, bw_url, bw_status, bw_mb, read_length, mapped_reads,
         auc_all, junction_count)
write.csv(avail, "results/fto_recount3_availability.csv", row.names = FALSE)
cat("Groups:", paste(names(table(avail$fto_group)), table(avail$fto_group),
                     collapse = ", "), "\n")
cat("Saved results/fto_recount3_availability.csv\n")

# Download one recount3 file over https into data/recount3/raw and check it is
# real data: a redirect or error page is a few hundred bytes of HTML, and a
# real .gz file starts with the gzip magic bytes 1f 8b.
raw_dir <- file.path(r3_dir, "raw"); dir.create(raw_dir, showWarnings = FALSE)
# Use the system CA bundle when one is set (needed behind some HTTPS proxies).
ca_cfg <- function() {
  ca <- Sys.getenv("CURL_CA_BUNDLE", Sys.getenv("SSL_CERT_FILE"))
  if (nzchar(ca)) httr::config(cainfo = ca) else httr::config()
}
fetch_r3 <- function(url) {
  url  <- sub("^http://", "https://", url)
  dest <- file.path(raw_dir, basename(url))
  is_gz <- function(f) file.exists(f) && file.size(f) > 100 &&
    identical(readBin(f, "raw", 2), as.raw(c(0x1f, 0x8b)))
  if (!is_gz(dest)) {
    for (attempt in 1:4) {
      ok <- tryCatch({ download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) FALSE)
      if (ok && is_gz(dest)) break
      Sys.sleep(2^attempt)
    }
  }
  if (!is_gz(dest)) stop("Download failed or not a gzip file: ", url)
  dest
}

# ---- 2. Groups ----------------------------------------------------------------
# Same cohort and NPM1 calls as Session 7. FTO-high / FTO-low = the k-means
# split of the three-gene NEK10/FTO/USP44 core score inside the 31 NPM1-mutated
# samples (results/tcga_block_scores.csv, core_group). All 80 wild-type
# samples are the reference group (approved).

grp_levels <- c("FTO_high", "FTO_low", "NPM1_wt")
grp_labels <- c(FTO_high = "NPM1-mut, FTO-high", FTO_low = "NPM1-mut, FTO-low",
                NPM1_wt = "NPM1 wild-type")
grp_cols   <- c(FTO_high = "#D55E00", FTO_low = "#0072B2", NPM1_wt = "grey55")
avail <- avail %>%
  left_join(blocks %>% select(patient, core_score, FTO_vst = FTO), by = "patient") %>%
  mutate(fto_group = factor(fto_group, levels = grp_levels))
stopifnot(identical(as.vector(table(avail$fto_group)), c(6L, 25L, 80L)))
hi <- avail$fto_group == "FTO_high"; lo <- avail$fto_group == "FTO_low"
wt <- avail$fto_group == "NPM1_wt"

# ---- 3. Coverage along the FTO locus ------------------------------------------
# Annotation: GENCODE v26 (Ensembl 88) on GRCh38, the exon annotation recount3
# itself uses. Canonical FTO = ENST00000471389 (9 exons), the transcript the
# IFEX9 paper used. Exon 9 is FTO's last exon and holds the long 3' UTR.

gtf <- import(fetch_r3(paste0("http://duffel.rail.bio/recount3/human/annotations/",
                              "exon_sums/human.exon_sums.G026.gtf.gz")))
fto_ex  <- gtf[gtf$gene_name == "FTO"]
canon   <- sort(fto_ex[fto_ex$transcript_id == "ENST00000471389.5" |
                       sub("\\..*", "", fto_ex$transcript_id) == "ENST00000471389"])
stopifnot(length(canon) == 9)
canon$exon <- paste0("exon ", seq_along(canon))
# Coverage is unstranded (both strands added up), so drop strand from the
# ranges; otherwise setdiff() treats "+" and "*" ranges as unrelated.
fto_ex <- unstrand(fto_ex); canon <- unstrand(canon)
fto_span <- range(fto_ex)
region   <- GRanges("chr16", IRanges(start(fto_span) - 20000, end(fto_span) + 20000))
reg_pos  <- start(region):end(region)
cat("\nFTO (GENCODE v26, GRCh38): chr16:", start(fto_span), "-", end(fto_span),
    " (", round(width(fto_span) / 1e3), " kb, ", n_distinct(fto_ex$transcript_id),
    " transcripts); plotted region chr16:", start(region), "-", end(region), "\n", sep = "")
cat("Canonical exons:\n"); print(data.frame(exon = canon$exon, start = start(canon),
                                            end = end(canon), bp = width(canon)), row.names = FALSE)

# IFEX9 (step 6 checks these against the genome): exon 1 sits in FTO intron 8,
# exon 2 is FTO exon 9. Coordinates from the paper's supplement, GRCh38.
ifex9_e1 <- GRanges("chr16", IRanges(54063753, 54063943))

# Read just the locus from each bigWig (a range request; the file has an index,
# so only the needed blocks travel). One small .rds per sample in data/.
cov_dir <- file.path(r3_dir, "fto_cov"); dir.create(cov_dir, showWarnings = FALSE)
get_cov <- function(i) {
  f <- file.path(cov_dir, paste0(avail$rail_id[i], ".rds"))
  if (file.exists(f)) return("cached")
  url <- sub("^http://", "https://", avail$bw_url[i])
  for (attempt in 1:4) {
    x <- tryCatch(import(BigWigFile(url), which = region, as = "RleList")[["chr16"]][reg_pos],
                  error = function(e) NULL)
    if (!is.null(x) && length(x) == length(reg_pos)) { saveRDS(x, f); return("downloaded") }
    Sys.sleep(2^attempt)
  }
  "failed"
}
todo <- sum(!file.exists(file.path(cov_dir, paste0(avail$rail_id, ".rds"))))
cat("Reading the FTO slice of", todo, "bigWig files (cached:", nrow(avail) - todo, ")\n")
invisible(parallel::mclapply(seq_len(nrow(avail)), get_cov, mc.cores = 8,
                              mc.preschedule = FALSE))
# A parallel worker can die without returning anything, so check the files
# themselves and redo any missing sample one at a time.
cov_files <- file.path(cov_dir, paste0(avail$rail_id, ".rds"))
for (i in which(!file.exists(cov_files))) get_cov(i)
if (!all(file.exists(cov_files))) stop("Coverage missing for ", sum(!file.exists(cov_files)),
                                       " samples")

# Depth adjustment (approved): recount3's standard, scale each sample so its
# total coverage (AUC) equals 40 million 100-base reads.
avail$scale <- 40e6 * 100 / avail$auc_all
cm <- vapply(seq_len(nrow(avail)), function(i)
  as.numeric(readRDS(file.path(cov_dir, paste0(avail$rail_id[i], ".rds")))) * avail$scale[i],
  numeric(length(reg_pos)))
colnames(cm) <- avail$patient
cat("Coverage matrix:", nrow(cm), "bases x", ncol(cm), "samples\n")
idx <- function(s, e) (s:e) - start(region) + 1
mean_cov <- function(gr) {      # mean coverage over all bases of a GRanges
  ii <- unlist(lapply(seq_along(gr), function(k) idx(start(gr)[k], end(gr)[k])))
  colMeans(cm[ii, , drop = FALSE])
}

pc <- 0.05   # pseudocount for log2: about one 100-base read per 2 kb at this depth
bin_cov <- function(width_bp, from = start(region), to = end(region)) {
  br <- seq(from, to + width_bp, by = width_bp)
  b  <- findInterval(reg_pos, br)
  keep <- reg_pos >= from & reg_pos <= to
  m  <- rowsum(cm[keep, ], b[keep]) / as.vector(table(b[keep]))
  data.frame(mid = br[as.integer(rownames(m))] + width_bp / 2, m, check.names = FALSE)
}
long_groups <- function(bins) {
  bins %>% tidyr::pivot_longer(-mid, names_to = "patient", values_to = "cov") %>%
    left_join(avail %>% select(patient, fto_group, core_score), by = "patient")
}

# Gene model track: canonical exons tall, other FTO exons short, other genes
# grey, IFEX9 exon 1 in red.
others <- gtf[seqnames(gtf) == "chr16" & overlapsAny(gtf, region) & gtf$gene_name != "FTO"]
gene_track <- function(xlim, label_exons = TRUE) {
  ex <- as.data.frame(reduce(fto_ex))
  ot <- as.data.frame(others) %>% filter(end >= xlim[1], start <= xlim[2]) %>%
    group_by(gene_name) %>% mutate(gs = min(start), ge = max(end)) %>% ungroup()
  ot_lab <- ot %>% distinct(gene_name, gs, ge, strand)
  cn <- as.data.frame(canon) %>% mutate(mid = (start + end) / 2)
  p <- ggplot() +
    annotate("segment", x = max(start(fto_span), xlim[1]), xend = min(end(fto_span), xlim[2]),
             y = 1, yend = 1, colour = "black", linewidth = 0.4) +
    geom_rect(data = ex, aes(xmin = start, xmax = end, ymin = 0.85, ymax = 1.15),
              fill = "grey40") +
    geom_rect(data = cn, aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1.3), fill = "black") +
    geom_rect(aes(xmin = start(ifex9_e1) - 300, xmax = end(ifex9_e1) + 300, ymin = 0.7, ymax = 1.3),
              fill = "red") +
    annotate("text", x = start(ifex9_e1), y = 1.55, label = "IFEX9 exon 1", colour = "red",
             size = 2.6) +
    geom_segment(data = ot_lab, aes(x = pmax(gs, xlim[1]), xend = pmin(ge, xlim[2]), y = 0.3, yend = 0.3),
                 colour = "grey60") +
    geom_rect(data = ot, aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.4), fill = "grey60") +
    geom_text(data = ot_lab, aes(x = pmin(pmax((gs + ge) / 2, xlim[1]), xlim[2]), y = 0.05,
                                 label = paste0(gene_name, " (", strand, ")")),
              size = 2, colour = "grey40", check_overlap = TRUE) +
    coord_cartesian(xlim = xlim, ylim = c(-0.1, 1.75), expand = FALSE) +
    labs(x = "chr16 position (GRCh38, Mb)", y = NULL) +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
  if (label_exons) p <- p + geom_text(data = cn %>% filter(mid >= xlim[1], mid <= xlim[2]),
                                      aes(x = mid, y = 1.5, label = sub("exon ", "", exon)),
                                      size = 2.4, check_overlap = TRUE)
  p
}
cov_panel <- function(lg, xlim, title, show_high = TRUE) {
  med <- lg %>% group_by(fto_group, mid) %>% summarise(cov = median(cov), .groups = "drop")
  p <- ggplot() +
    annotate("rect", xmin = start(ifex9_e1) - 300, xmax = end(ifex9_e1) + 300, ymin = pc / 2, ymax = Inf,
             fill = "red", alpha = 0.08)
  if (show_high) p <- p + geom_line(data = lg %>% filter(fto_group == "FTO_high"),
                                    aes(mid, cov + pc, group = patient),
                                    colour = grp_cols[["FTO_high"]], alpha = 0.35, linewidth = 0.3)
  p + geom_line(data = med, aes(mid, cov + pc, colour = fto_group), linewidth = 0.7) +
    scale_colour_manual(values = grp_cols, labels = paste0(grp_labels, " (median, n = ",
                                                           table(avail$fto_group), ")")) +
    scale_y_log10() +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    labs(x = NULL, y = "Coverage per base\n(40 M 100-base reads, log)", colour = NULL,
         title = title) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "top", axis.text.x = element_blank(), panel.grid.minor = element_blank())
}

full_x <- c(start(region), end(region))
zoom_x <- c(54050000, 54125000)
lg_full <- long_groups(bin_cov(500))
lg_zoom <- long_groups(bin_cov(100, zoom_x[1], zoom_x[2]))
p_full <- cov_panel(lg_full, full_x, paste("A. Whole FTO locus +/- 20 kb, 500 bp bins;",
                                           "thin orange lines = each FTO-high sample"))
p_zoom <- cov_panel(lg_zoom, zoom_x, "B. Zoom: end of intron 8 (IFEX9 exon 1) to exon 9, 100 bp bins")
fig <- p_full / gene_track(full_x) / p_zoom / gene_track(zoom_x) +
  plot_layout(heights = c(3, 1, 3, 1))
ggsave("figures/11_fto_coverage.png", fig, width = 11, height = 9, dpi = 300)
cat("Saved figures/11_fto_coverage.png\n")

# Each NPM1-mutated sample on its own (1 kb bins), wild-type median in grey.
lg_1k  <- long_groups(bin_cov(1000))
wt_med <- lg_1k %>% filter(fto_group == "NPM1_wt") %>% group_by(mid) %>%
  summarise(cov = median(cov), .groups = "drop")
ps <- lg_1k %>% filter(fto_group != "NPM1_wt") %>%
  mutate(panel = paste0(patient, ", ", ifelse(fto_group == "FTO_high", "high", "low"),
                        " (score ", sprintf("%+.1f", core_score), ")")) %>%
  mutate(panel = reorder(panel, -core_score))
p_ind <- ggplot(ps, aes(mid, cov + pc)) +
  annotate("rect", xmin = start(ifex9_e1) - 1500, xmax = end(ifex9_e1) + 1500, ymin = pc / 2, ymax = Inf,
           fill = "red", alpha = 0.15) +
  geom_line(data = wt_med, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(colour = fto_group), linewidth = 0.35) +
  scale_colour_manual(values = grp_cols, labels = grp_labels) +
  scale_y_log10() +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6), breaks = c(53.75e6, 53.95e6, 54.15e6)) +
  facet_wrap(~ panel, ncol = 4) +
  labs(x = "chr16 (GRCh38, Mb)", y = "Coverage per base (log)", colour = NULL,
       title = "FTO locus coverage in every NPM1-mutated TCGA sample (1 kb bins)",
       subtitle = paste("Ordered by NEK10/FTO/USP44 score. Grey = median of 80 wild-type.",
                        "Red band = IFEX9 exon 1; FTO exon 9 is the block at the right end.")) +
  theme_minimal(base_size = 8) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/11_fto_coverage_per_sample.png", p_ind, width = 11, height = 14, dpi = 300)
cat("Saved figures/11_fto_coverage_per_sample.png\n")

# ---- 4. Exons and 1 kb windows: where do FTO-high and FTO-low differ? --------
# Features: the 9 canonical FTO exons; other FTO exon pieces (alternative
# transcripts); FTO introns (gene body minus every FTO exon) and the 20 kb
# flanks cut into windows of about 1 kb; and IFEX9 exon 1 on its own.
all_fto_ex <- reduce(fto_ex)
alt_ex  <- setdiff(all_fto_ex, canon)
introns <- setdiff(fto_span, all_fto_ex)
flanks  <- setdiff(region, fto_span)
canon_introns <- GRanges("chr16", IRanges(end(canon)[-9] + 1, start(canon)[-1] - 1))
win <- function(gr, type) {
  t <- unlist(tile(gr, width = 1000))
  t$type <- type; t
}
feat <- c(
  { x <- canon; mcols(x) <- NULL; x$type <- "exon"; x$label <- canon$exon; x },
  { x <- alt_ex; x$type <- "alt_exon"
    o <- findOverlaps(x, fto_ex)
    x$label <- paste0("alt exon ", seq_along(x), " (", vapply(seq_along(x), function(k)
      paste(sort(unique(fto_ex$transcript_name[subjectHits(o)[queryHits(o) == k]])), collapse = ","), ""), ")")
    x },
  { x <- win(introns, "intron"); x$label <- NA; x },
  { x <- win(flanks, "flank"); x$label <- NA; x },
  { x <- ifex9_e1; x$type <- "IFEX9_exon1"; x$label <- "IFEX9 exon 1"; x })
feat$intron <- NA_character_
ov <- findOverlaps(feat, canon_introns, select = "first")
feat$intron[!is.na(ov)] <- paste0("intron ", ov[!is.na(ov)])
feat$label <- ifelse(is.na(feat$label),
                     paste0(ifelse(feat$type == "flank",
                                   ifelse(start(feat) < start(fto_span), "5' flank", "3' flank"),
                                   ifelse(!is.na(feat$intron), feat$intron,
                                          ifelse(start(feat) < start(canon)[1], "FTO gene body before exon 1",
                                                 "FTO gene body after exon 9"))),
                            " ", round(start(feat) / 1e3), " kb"), feat$label)
oth <- findOverlaps(feat, others)
feat$other_gene_exon <- vapply(seq_along(feat), function(k)
  paste(unique(others$gene_name[subjectHits(oth)[queryHits(oth) == k]]), collapse = ";"), "")
feat$contains_ifex9_e1 <- overlapsAny(feat, ifex9_e1) & feat$type == "intron"
cat("\nFeatures: ", paste(names(table(feat$type)), table(feat$type), collapse = ", "), "\n", sep = "")

fc <- t(vapply(seq_along(feat), function(k) mean_cov(feat[k]), numeric(ncol(cm))))
colnames(fc) <- avail$patient
exon_all <- mean_cov(canon)          # base-weighted, so mostly exon 9 (10.2 of 11.4 kb)
exon_1_8 <- mean_cov(canon[1:8])     # the part of FTO that IFEX9 doesn't share
abs_l2  <- log2(fc + pc)
rat_all <- log2(sweep(fc + pc, 2, exon_all + pc, "/"))
rat_18  <- log2(sweep(fc + pc, 2, exon_1_8 + pc, "/"))

test_rows <- function(m, what) {
  r <- data.frame(what = what,
    med_high = apply(m[, hi], 1, median), med_low = apply(m[, lo], 1, median),
    med_wt = apply(m[, wt], 1, median))
  r$diff_high_low <- r$med_high - r$med_low
  r$p <- apply(m, 1, function(v) suppressWarnings(wilcox.test(v[hi], v[lo])$p.value))
  r$padj <- p.adjust(r$p, "BH")
  r
}
fdf <- as.data.frame(feat) %>% select(seqnames, start, end, width, type, label, intron,
                                      other_gene_exon, contains_ifex9_e1)
win_tab <- bind_rows(
  cbind(fdf, test_rows(abs_l2, "log2 coverage")),
  cbind(fdf, test_rows(rat_all, "log2 ratio to mean FTO exon coverage (exons 1-9)")),
  cbind(fdf, test_rows(rat_18, "log2 ratio to FTO exons 1-8 coverage")))
win_tab$mean_cov_high <- rep(rowMeans(fc[, hi]), 3)
win_tab$mean_cov_low  <- rep(rowMeans(fc[, lo]), 3)
write.csv(win_tab, "results/fto_window_tests.csv", row.names = FALSE)
cat("Saved results/fto_window_tests.csv\n")
show_top <- function(w) {
  x <- win_tab %>% filter(what == w) %>% arrange(p) %>% head(12) %>%
    transmute(label, type, start, mean_high = round(mean_cov_high, 2), mean_low = round(mean_cov_low, 2),
              diff = round(diff_high_low, 2), p = signif(p, 2), padj = signif(padj, 2))
  cat("\nTop features,", w, "(diff = median high - median low, log2):\n"); print(x, row.names = FALSE)
  s <- win_tab %>% filter(what == w)
  cat("padj < 0.05:", sum(s$padj < 0.05), "of", nrow(s), "| up in high:",
      sum(s$padj < 0.05 & s$diff_high_low > 0), "| exons among them:",
      sum(s$padj < 0.05 & s$type == "exon"), "\n")
}
for (w in unique(win_tab$what)) show_top(w)
cat("\nCanonical exons (log2 coverage and log2 ratio to exons 1-8):\n")
print(win_tab %>% filter(type %in% c("exon", "IFEX9_exon1"), what != "log2 ratio to mean FTO exon coverage (exons 1-9)") %>%
        transmute(label, what = ifelse(grepl("1-8", what), "ratio_1_8", "abs"),
                  diff = round(diff_high_low, 2), p = signif(p, 2), padj = signif(padj, 2)),
      row.names = FALSE)

# Where along the locus do the groups differ? Both versions, window by window.
pw <- win_tab %>% filter(what != "log2 ratio to FTO exons 1-8 coverage") %>%
  mutate(mid = (start + end) / 2, sig = ifelse(padj < 0.05, "padj < 0.05", "not significant"),
         what = factor(what, levels = unique(what)))
p_w <- ggplot(pw, aes(mid, diff_high_low)) +
  annotate("rect", xmin = start(ifex9_e1) - 1000, xmax = end(ifex9_e1) + 1000, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.15) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(aes(colour = sig, shape = type == "exon" | type == "IFEX9_exon1"), size = 1.1) +
  scale_colour_manual(values = c("padj < 0.05" = "#D55E00", "not significant" = "grey60")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("window", "exon")) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6)) +
  facet_wrap(~ what, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "Median FTO-high minus median FTO-low (log2)", colour = NULL, shape = NULL,
       title = "Where FTO-high and FTO-low differ along the FTO locus (TCGA, 6 vs 25)") +
  theme_minimal(base_size = 9) + theme(legend.position = "top")
ggsave("figures/11_fto_window_tests.png", p_w / gene_track(full_x) + plot_layout(heights = c(5, 1)),
       width = 11, height = 7, dpi = 300)
cat("Saved figures/11_fto_window_tests.png\n")

# ---- 5. Splice junctions inside the FTO locus ---------------------------------
# recount3 stores junctions for the whole project in three files: RR (one row
# per junction: chromosome, intron start and end, strand, and whether it is in
# recount3's annotation list), MM (a sparse matrix: junction row, sample column,
# read count) and ID (the rail_id of each sample column).
rr_file <- fetch_r3(locate_url("LAML", proj$project_home, type = "jxn")[["tcga.junctions.LAML.ALL.RR.gz"]])
mm_file <- fetch_r3(locate_url("LAML", proj$project_home, type = "jxn")[["tcga.junctions.LAML.ALL.MM.gz"]])
id_file <- fetch_r3(locate_url("LAML", proj$project_home, type = "jxn")[["tcga.junctions.LAML.ALL.ID.gz"]])
rr <- fread(rr_file)
rr[, row := .I]
jl <- rr[chromosome == "chr16" & start >= start(region) & end <= end(region)]
ids <- fread(id_file)$rail_id
mm_head <- readLines(gzfile(mm_file), n = 3)
stopifnot(grepl("MatrixMarket", mm_head[1]),
          as.numeric(strsplit(mm_head[3], "\t| ")[[1]][1:2]) == c(nrow(rr), length(ids)))
mm <- fread(mm_file, skip = 3, header = FALSE, col.names = c("row", "col", "count"))
mm <- mm[row %in% jl$row]
mm[, rail_id := ids[col]]
mm <- mm[rail_id %in% avail$rail_id]
cat("\nJunctions in the locus (any sample in recount3):", nrow(jl),
    "| with reads in our 111 samples:", uniqueN(mm$row), "\n")

jm <- matrix(0, nrow(jl), nrow(avail), dimnames = list(jl$row, avail$patient))
jm[cbind(match(mm$row, jl$row), match(mm$rail_id, avail$rail_id))] <- mm$count

# Our own annotation check against every GENCODE v26 transcript in the region:
# an intron = the gap between two neighbouring exons of one transcript.
reg_tx <- gtf[overlapsAny(gtf, region)]
tx_introns <- do.call(rbind, lapply(split(as.data.frame(reg_tx), reg_tx$transcript_id), function(d) {
  d <- d[order(d$start), ]
  if (nrow(d) < 2) return(NULL)
  data.frame(start = d$end[-nrow(d)] + 1, end = d$start[-1] - 1, gene = d$gene_name[1],
             transcript = d$transcript_id[1])
}))
ex_ends <- unique(end(reg_tx)); ex_starts <- unique(start(reg_tx))
canon_intr <- data.frame(start = end(canon)[-9] + 1, end = start(canon)[-1] - 1,
                         name = paste0("FTO exon ", 1:8, " -> ", 2:9))
jt <- data.frame(row = jl$row, start = jl$start, end = jl$end, length = jl$length,
                 strand = jl$strand, recount3_annotated = jl$annotated == 1,
                 motif = paste0(jl$left_motif, "-", jl$right_motif))
key <- function(s, e) paste(s, e)
jt$gencode26 <- key(jt$start, jt$end) %in% key(tx_introns$start, tx_introns$end)
jt$donor_known    <- (jt$start - 1) %in% ex_ends     # + strand view: exon end before the intron
jt$acceptor_known <- (jt$end + 1) %in% ex_starts
jt$class <- case_when(jt$gencode26 ~ "annotated (GENCODE v26)",
                      jt$donor_known & jt$acceptor_known ~ "new pairing of known exon ends",
                      jt$donor_known | jt$acceptor_known ~ "one end new",
                      TRUE ~ "both ends new")
jt$canonical <- canon_intr$name[match(key(jt$start, jt$end), key(canon_intr$start, canon_intr$end))]
# Within 2 bp: the paper's exon 1 end (54063943) is 1 bp past the splice donor
# that GENCODE (FTO-007) and the junction reads put at 54063942, most likely a
# coordinate-convention difference.
jt$touches_ifex9_e1 <- abs(jt$start - 1 - end(ifex9_e1)) <= 2 | abs(jt$end + 1 - start(ifex9_e1)) <= 2
jt$ends_at_fto_exon9 <- jt$end + 1 == start(canon)[9]
jt$genes <- vapply(seq_len(nrow(jt)), function(k) {
  g <- unique(tx_introns$gene[tx_introns$start <= jt$end[k] & tx_introns$end >= jt$start[k]])
  if (length(g)) paste(g, collapse = ";") else "" }, "")
jt$inside_fto <- jt$start >= start(fto_span) & jt$end <= end(fto_span)
cat("Coordinate check, our GENCODE v26 match vs recount3 annotated flag:\n"); print(table(ours = jt$gencode26, recount3 = jt$recount3_annotated))

# Reads per 10 million mapped reads, so deeper samples don't look "higher".
jn <- sweep(jm, 2, avail$mapped_reads / 1e7, "/")
jt$total_reads <- rowSums(jm)
jt$n_high <- rowSums(jm[, hi, drop = FALSE] >= 2)   # samples with at least 2 reads
jt$n_low  <- rowSums(jm[, lo, drop = FALSE] >= 2)
jt$n_wt   <- rowSums(jm[, wt, drop = FALSE] >= 2)
jt$reads_high <- rowSums(jm[, hi]); jt$reads_low <- rowSums(jm[, lo]); jt$reads_wt <- rowSums(jm[, wt])
jt$med_norm_high <- apply(jn[, hi], 1, median); jt$med_norm_low <- apply(jn[, lo], 1, median)
jt$mean_norm_high <- rowMeans(jn[, hi]); jt$mean_norm_low <- rowMeans(jn[, lo])
jt$mean_norm_wt <- rowMeans(jn[, wt])
jt$tested <- tested <- jt$total_reads >= 10
jt$p <- NA_real_
jt$p[tested] <- apply(jn[tested, , drop = FALSE], 1, function(v)
  suppressWarnings(wilcox.test(v[hi], v[lo])$p.value))
jt$padj <- NA_real_; jt$padj[tested] <- p.adjust(jt$p[tested], "BH")
# "Mostly in FTO-high" (my rule, set before looking): seen (>= 2 reads) in at
# least half of the 6 FTO-high samples and in at most 20% of the 25 FTO-low.
jt$mostly_high <- jt$n_high >= 3 & jt$n_low <= 5
jt <- jt %>% arrange(p)
tested <- jt$tested   # keep the flag in the same order as the sorted rows
write.csv(jt %>% select(-row), "results/fto_junctions.csv", row.names = FALSE)
cat("Saved results/fto_junctions.csv\n")
cat("Junctions with >= 10 reads in our samples:", sum(tested), "(inside FTO:",
    sum(tested & jt$inside_fto), ")\n")
cat("By class (>= 10 reads, inside FTO):\n")
print(table(jt$class[tested & jt$inside_fto]))
cat("\nJunctions that show up mostly in FTO-high (seen in >= 3 of 6 high, <= 5 of 25 low):\n")
print(jt %>% filter(mostly_high) %>%
        transmute(start, end, length, strand, class, touches_ifex9_e1, ends_at_fto_exon9,
                  n_high, n_low, n_wt, reads_high, reads_low, reads_wt, p = signif(p, 2),
                  padj = signif(padj, 2)), row.names = FALSE)
cat("\nTop 10 by Wilcoxon p (high vs low):\n")
print(head(jt %>% filter(tested) %>%
             transmute(start, end, length, class, canonical, n_high, n_low, n_wt,
                       mean_high = round(mean_norm_high, 2), mean_low = round(mean_norm_low, 2),
                       p = signif(p, 2), padj = signif(padj, 2)), 10), row.names = FALSE)
cat("\nCanonical FTO junctions (exon to exon), mean reads per 10 M mapped:\n")
print(jt %>% filter(!is.na(canonical)) %>% arrange(start) %>%
        transmute(canonical, mean_high = round(mean_norm_high, 1), mean_low = round(mean_norm_low, 1),
                  mean_wt = round(mean_norm_wt, 1), ratio_high_low = round(mean_norm_high / mean_norm_low, 2),
                  p = signif(p, 2)), row.names = FALSE)

# Junction figure: each junction inside FTO with >= 10 reads, drawn as an arc
# from donor to acceptor, height = share of samples in the group where it's seen.
arcs <- jt %>% filter(tested, inside_fto) %>%
  tidyr::pivot_longer(c(n_high, n_low, n_wt), names_to = "grp", values_to = "n") %>%
  mutate(share = n / c(n_high = 6, n_low = 25, n_wt = 80)[grp],
         fto_group = factor(c(n_high = "FTO_high", n_low = "FTO_low", n_wt = "NPM1_wt")[grp], levels = grp_levels),
         kind = ifelse(class == "annotated (GENCODE v26)", "annotated (GENCODE v26)", "not annotated"),
         id = paste(start, end)) %>%
  filter(share > 0)
# Each junction is a half-arch from donor to acceptor; its height is the share
# of samples in that group with at least 2 reads on it.
arc_pts <- arcs %>% rowwise() %>%
  reframe(id, fto_group, kind, share, t = seq(0, 1, length.out = 60)) %>%
  left_join(arcs %>% distinct(id, start, end), by = "id") %>%
  mutate(x = start + t * (end - start), y = share * sin(pi * t))
p_j <- ggplot(arc_pts, aes(x, y, group = id, colour = kind)) +
  geom_vline(xintercept = start(ifex9_e1), colour = "red", linewidth = 0.3, alpha = 0.5) +
  geom_path(linewidth = 0.5, alpha = 0.8) +
  scale_colour_manual(values = c(`annotated (GENCODE v26)` = "grey35", `not annotated` = "#CC79A7")) +
  facet_wrap(~ fto_group, ncol = 1, labeller = as_labeller(grp_labels)) +
  coord_cartesian(xlim = full_x, ylim = c(0, 1.05), expand = FALSE) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6)) +
  labs(x = NULL, y = "Share of samples with >= 2 reads", colour = NULL,
       title = "Splice junctions inside FTO (>= 10 reads in total over the 111 samples)",
       subtitle = "Each arch joins a donor to an acceptor; red line = IFEX9 exon 1") +
  theme_minimal(base_size = 9) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/11_fto_junctions.png", p_j / gene_track(full_x) + plot_layout(heights = c(6, 1)),
       width = 11, height = 7, dpi = 300)
cat("Saved figures/11_fto_junctions.png\n")

# ---- 6. IFEX9 coordinates ----------------------------------------------------
# The supplement (Arza-Apalategi et al., Leukemia 2023) gives IFEX9 exon 1 as
# chr16:54063753-54063943 (GRCh38, used for TCGA) and 54097665-54097855
# (GRCh37, Beat AML), and says IFEX9 exon 2 = FTO exon 9. Check this against
# the genome: their qPCR forward primer should sit inside exon 1 and their
# reverse primer (read backwards, on the other strand) inside FTO exon 9.
genome_seq <- function(s, e) {   # Ensembl REST API, GRCh38, 1-based, + strand
  r <- httr::GET(sprintf("https://rest.ensembl.org/sequence/region/human/16:%d..%d:1?content-type=text/plain",
                         s, e), ca_cfg())
  httr::stop_for_status(r)
  toupper(httr::content(r, as = "text", encoding = "UTF-8"))
}
revcomp <- function(x) chartr("ACGT", "TGCA", paste(rev(strsplit(x, "")[[1]]), collapse = ""))
fwd <- "CTCATCTGTCAAGCCCTTCTTCTG"; rev_p <- "GCCGACATTCTGGCTTCTGATC"; race <- "CTGAAACGATGTCTGTGAGGTC"
chk <- tryCatch({
  e1 <- genome_seq(start(ifex9_e1), end(ifex9_e1))
  e9 <- genome_seq(start(canon)[9], start(canon)[9] + 1999)
  data.frame(check = c("forward primer inside IFEX9 exon 1 (+ strand)",
                       "reverse qPCR primer inside FTO exon 9 (first 2 kb, - strand)",
                       "5' RACE primer inside FTO exon 9 (first 2 kb, - strand)"),
             found = c(grepl(fwd, e1), grepl(revcomp(rev_p), e9), grepl(revcomp(race), e9)),
             position = c(start(ifex9_e1) + regexpr(fwd, e1)[1] - 1,
                          start(canon)[9] + regexpr(revcomp(rev_p), e9)[1] - 1,
                          start(canon)[9] + regexpr(revcomp(race), e9)[1] - 1))
}, error = function(e) data.frame(check = "Ensembl sequence lookup", found = NA,
                                  position = NA_real_, note = conditionMessage(e)))
cat("\nIFEX9 coordinate check against the hg38 sequence:\n"); print(chk, row.names = FALSE)
ifex9_tab <- data.frame(
  feature = c("IFEX9 exon 1", "IFEX9 exon 1", "IFEX9 exon 2 = FTO exon 9 (GENCODE v26)"),
  build = c("GRCh38", "GRCh37", "GRCh38"), chrom = "chr16",
  start = c(54063753, 54097665, start(canon)[9]), end = c(54063943, 54097855, end(canon)[9]),
  source = c("Arza-Apalategi et al. 2023, Supplementary methods", "same", "GENCODE v26 ENST00000471389"))
write.csv(ifex9_tab, "results/fto_ifex9_coordinates.csv", row.names = FALSE)
write.csv(chk, "results/fto_ifex9_primer_check.csv", row.names = FALSE)
cat("Saved results/fto_ifex9_coordinates.csv and results/fto_ifex9_primer_check.csv\n")

# The junction IFEX9 needs: from the end of its exon 1 straight to FTO exon 9.
ifx_j <- jt %>% filter(abs(start - 1 - end(ifex9_e1)) <= 2, end == start(canon)[9] - 1)
stopifnot(nrow(ifx_j) <= 1)
cat("\nIFEX9 exon 1 -> FTO exon 9 junction: ", if (nrow(ifx_j)) paste0("found at intron ", ifx_j$start,
                                     "-", ifx_j$end, "; class '", ifx_j$class, "'; seen in ", ifx_j$n_high,
                                     " of 6 high, ", ifx_j$n_low, " of 25 low, ", ifx_j$n_wt,
                                     " of 80 wild-type") else "not found", "\n", sep = "")

# Per-sample summary: the numbers behind the conclusion, one row per sample.
ifx_reads <- if (nrow(ifx_j)) jm[as.character(ifx_j$row), ] else rep(0, ncol(jm))
canon89 <- jt %>% filter(canonical == "FTO exon 8 -> 9")
per_sample <- avail %>% select(patient, fto_group, core_score, FTO_vst, FLT3_ITD, mapped_reads) %>%
  mutate(cov_exons_1_8 = exon_1_8, cov_exon_9 = mean_cov(canon[9]),
         cov_ifex9_exon1 = mean_cov(ifex9_e1), cov_intron8 = mean_cov(introns[overlapsAny(introns, canon_introns[8])]),
         cov_other_introns = mean_cov(introns[!overlapsAny(introns, canon_introns[8])]),
         log2_exon9_vs_1_8 = log2((cov_exon_9 + pc) / (cov_exons_1_8 + pc)),
         log2_ifex9e1_vs_1_8 = log2((cov_ifex9_exon1 + pc) / (cov_exons_1_8 + pc)),
         jx_ifex9_to_exon9 = as.numeric(ifx_reads),
         jx_exon8_to_exon9 = if (nrow(canon89)) as.numeric(jm[as.character(canon89$row), ]) else NA)
write.csv(per_sample, "results/fto_locus_per_sample.csv", row.names = FALSE)
cat("Saved results/fto_locus_per_sample.csv\n")
cat("\nPer-sample key numbers, NPM1-mutated, highest score first:\n")
print(per_sample %>% filter(fto_group != "NPM1_wt") %>% arrange(desc(core_score)) %>%
        transmute(patient, group = fto_group, score = round(core_score, 2),
                  ex1_8 = round(cov_exons_1_8, 1), ex9 = round(cov_exon_9, 1),
                  ifex9_e1 = round(cov_ifex9_exon1, 1), intron8 = round(cov_intron8, 2),
                  other_introns = round(cov_other_introns, 2),
                  jx_ifex9 = jx_ifex9_to_exon9, jx_8_9 = jx_exon8_to_exon9), row.names = FALSE)
cat("\nWild-type medians: exons 1-8", round(median(per_sample$cov_exons_1_8[wt]), 1),
    "| exon 9", round(median(per_sample$cov_exon_9[wt]), 1),
    "| IFEX9 exon 1", round(median(per_sample$cov_ifex9_exon1[wt]), 2),
    "| IFEX9 junction reads > 0 in", sum(per_sample$jx_ifex9_to_exon9[wt] > 0), "of 80\n")
key_tests <- bind_rows(lapply(c("cov_exons_1_8", "cov_exon_9", "cov_ifex9_exon1", "cov_intron8",
                                "cov_other_introns", "log2_exon9_vs_1_8", "log2_ifex9e1_vs_1_8",
                                "jx_ifex9_to_exon9"), function(v) {
  x <- per_sample[[v]]; if (!grepl("^log2", v)) x <- log2(x + pc)
  data.frame(measure = v, med_high = median(x[hi]), med_low = median(x[lo]), med_wt = median(x[wt]),
             p_high_vs_low = suppressWarnings(wilcox.test(x[hi], x[lo])$p.value),
             p_low_vs_wt = suppressWarnings(wilcox.test(x[lo], x[wt])$p.value))
}))
write.csv(key_tests, "results/fto_key_feature_tests.csv", row.names = FALSE)
cat("\nKey features (log2 scale, coverage + pseudocount", pc, "):\n"); print(key_tests, digits = 3, row.names = FALSE)

kf <- per_sample %>%
  transmute(patient, fto_group,
            `FTO exons 1-8` = log2(cov_exons_1_8 + pc), `FTO exon 9 (= IFEX9 exon 2)` = log2(cov_exon_9 + pc),
            `IFEX9 exon 1` = log2(cov_ifex9_exon1 + pc), `Rest of FTO introns` = log2(cov_other_introns + pc),
            `Exon 9 / exons 1-8` = log2_exon9_vs_1_8, `IFEX9 exon 1 / exons 1-8` = log2_ifex9e1_vs_1_8,
            `IFEX9 -> exon 9 junction reads` = log2(jx_ifex9_to_exon9 + 1)) %>%
  tidyr::pivot_longer(-c(patient, fto_group)) %>%
  mutate(name = factor(name, levels = unique(name)))
p_k <- ggplot(kf, aes(fto_group, value, colour = fto_group)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50", width = 0.5) +
  geom_jitter(width = 0.15, height = 0, size = 1.1, alpha = 0.8) +
  scale_colour_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c(FTO_high = "high\n(n=6)", FTO_low = "low\n(n=25)", NPM1_wt = "wt\n(n=80)")) +
  facet_wrap(~ name, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "log2 (coverage, ratio, or junction reads + 1)",
       title = "FTO locus, one dot per TCGA sample",
       subtitle = "Coverage scaled to 40 M 100-base reads; ratios are log2 of mean coverage") +
  theme_minimal(base_size = 9)
ggsave("figures/11_fto_key_features.png", p_k, width = 11, height = 6, dpi = 300)
cat("Saved figures/11_fto_key_features.png\n")

# ---- 7. Beat AML: annotated genes inside FTO ---------------------------------
# Beat AML only has gene counts (Ensembl 75, GRCh37). Which annotated genes in
# the counts file sit inside the FTO locus, does any of them overlap IFEX9
# exon 1 (GRCh37 chr16:54097665-54097855), and do they track the Beat AML
# NEK10/FTO/USP44 high/low split (37 high, 63 low)?
fto_locus <- read.csv("results/fto_locus_genes_grch37.csv")
# Does Beat AML's annotation (Ensembl 75, GRCh37) already have the two FTO
# transcripts that start in intron 8 (FTO-007 = IFEX9, FTO-021)? If yes, their
# reads count towards the FTO gene in our counts file.
ens75 <- bind_rows(lapply(c("ENST00000472835", "ENST00000636091"), function(tx) {
  r <- tryCatch(httr::GET(paste0("https://grch37.rest.ensembl.org/lookup/id/", tx,
                                 "?expand=1;content-type=application/json"), ca_cfg()),
                error = function(e) NULL)
  if (is.null(r) || httr::status_code(r) != 200)
    return(data.frame(transcript = tx, in_grch37_ensembl = FALSE, gene = NA, exons = NA))
  j <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))
  data.frame(transcript = tx, in_grch37_ensembl = TRUE, gene = j$Parent,
             exons = paste(paste0(j$Exon$start, "-", j$Exon$end), collapse = "; "))
}))
cat("\nIntron-8 FTO transcripts in the GRCh37 Ensembl annotation (grch37.rest.ensembl.org):\n")
print(ens75, row.names = FALSE)
write.csv(ens75, "results/fto_ifex9_in_grch37_annotation.csv", row.names = FALSE)
counts <- fread("data/beataml_waves1to4_counts_dbgap.txt")
bblock <- read.csv("results/block_nek10_fto_scores.csv")
stopifnot(all(bblock$dbgap_rnaseq_sample %in% names(counts)), nrow(bblock) == 100)
cmat <- as.matrix(counts[, bblock$dbgap_rnaseq_sample, with = FALSE])
rownames(cmat) <- counts$stable_id
sf <- estimateSizeFactors(DESeqDataSetFromMatrix(cmat, data.frame(row.names = colnames(cmat),
                                                                  g = factor(bblock$block_group)), ~ g))
ncnt <- counts(sf, normalized = TRUE)
bh <- bblock$block_group == "high"; bl <- bblock$block_group == "low"
ba <- bind_rows(lapply(seq_len(nrow(fto_locus)), function(k) {
  id <- fto_locus$id[k]
  v <- if (id %in% rownames(ncnt)) ncnt[id, ] else rep(NA, ncol(ncnt))
  r <- if (id %in% rownames(cmat)) cmat[id, ] else rep(NA, ncol(cmat))
  data.frame(id = id, gene = fto_locus$external_name[k], biotype = fto_locus$biotype[k],
             start = fto_locus$start[k], end = fto_locus$end[k], strand = fto_locus$strand[k],
             overlaps_ifex9_exon1 = fto_locus$start[k] <= 54097855 & fto_locus$end[k] >= 54097665 &
               fto_locus$external_name[k] != "FTO",
             median_raw_high = median(r[bh]), median_raw_low = median(r[bl]),
             n_high_ge5 = sum(r[bh] >= 5), n_low_ge5 = sum(r[bl] >= 5),
             median_log2_high = median(log2(v[bh] + 1)), median_log2_low = median(log2(v[bl] + 1)),
             p = if (all(is.na(v)) || sd(v) == 0) NA else suppressWarnings(wilcox.test(v[bh], v[bl])$p.value),
             spearman_vs_core_score = if (all(is.na(v)) || sd(v) == 0) NA else
               cor(v, bblock$core_score, method = "spearman"))
}))
ba$padj <- p.adjust(ba$p, "BH")
write.csv(ba, "results/fto_locus_beataml_tests.csv", row.names = FALSE)
cat("\nBeat AML, genes in the FTO locus vs the high/low split (37 vs 63):\n")
print(ba %>% transmute(gene, biotype, overlaps_ifex9_exon1, median_raw_high, median_raw_low,
                       n_high_ge5, n_low_ge5, diff_log2 = round(median_log2_high - median_log2_low, 2),
                       p = signif(p, 2), padj = signif(padj, 2), rho = round(spearman_vs_core_score, 2)),
      row.names = FALSE)
cat("Saved results/fto_locus_beataml_tests.csv\n")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_11_fto_locus.txt")
cat("\nSaved results/sessionInfo_11_fto_locus.txt\n")
