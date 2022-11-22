# ggtranscript: transcript visualization
library(ggplot2)
library(ggtranscript)
library(rtracklayer)

gtf <- rtracklayer::import("/data1/wangy/data/ref/Gh/Ghirsutum_HAU_gene_model.gtf")
# class(gtf)
gtfTable <- gtf %>% as_tibble()
# head(gtfTable)

# Ghir_D04G000830_gtf_exon <- gtfTable %>% filter(gene_id == "Ghir_D04G000830", type == "exon")
# Ghir_D04G000830_gtf_cds <- gtfTable %>% filter(gene_id == "Ghir_D04G000830", type == "CDS")
# head(Ghir_D04G000830_gtf)
# Ghir_D04G000830_gtf_exon %>% 
# ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
# geom_range(fill = "white", height = 0.25) +
# geom_range(data = Ghir_D04G000830_gtf_cds) +
# geom_intron(data = to_intron(Ghir_D04G000830_gtf_exon, "transcript_id"),
# aes(strand = strand), arrow.min.intron.length = 500)

# ggsave("figure/alter/Ghir_D04G000830_gtf.pdf", width = 6, height = 4)

transcriptStructure  <- function(transcript) {
    exon <- gtfTable %>% filter(gene_id == transcript, type == "exon")
    cds <- gtfTable %>% filter(gene_id == transcript, type == "CDS")
    exon %>% 
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(fill = "white", height = 0.25) +
    geom_range(data = cds) +
    geom_intron(data = to_intron(exon, "transcript_id"),
    aes(strand = strand), arrow.min.intron.length = 500)
}
transcriptStructure("Ghir_D04G000830")
ggsave("figure/alter/Ghir_D04G000830_gtf.pdf", width = 6, height = 4)
transcriptStructure("Ghir_A01G000210")
ggsave("figure/alter/Ghir_A01G000210_gtf.pdf", width = 6, height = 4)
