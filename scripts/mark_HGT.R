suppressMessages(library(optparse))

option_list = list(
    make_option(c("-d", "--dl"), type="character", default=NULL,
                help="generax DL rates file", dest = "dl_file"),
    make_option(c("-t", "--dtl"), type="character", default=NULL,
                help="generax DTL rates file", dest = "dtl_file"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="output filename", dest = "output"),
    make_option(c("-r", "--reco"), type="character", default=NULL,
                help="output DTL GeneRax reconciliation directory", dest = "reco"),
    make_option(c("-f", "--files"), type="character", default=NULL,
                help="files to plot with thirdkind", dest = "tk"),
    make_option(c("-s", "--signif"), type="double", default=0.05, 
                help="p-value threshold for significance", dest = "alpha")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressMessages(library(tidyverse))

theme_set(theme_bw(base_family = "Helvetica"))

df_DL <- 2
df_DTL <- 3

dl <- read_delim(opt$dl_file, show_col_types = F)
dtl <- read_delim(opt$dtl_file, show_col_types = F)

df <- left_join(select(dl, gene, RecLL, D_rate, L_rate),
                select(dtl, gene, RecLL, contains("rate")), 
                by="gene", suffix = c("_DL", "_DTL")) %>% 
    mutate(diff = RecLL_DL-RecLL_DTL,
           LR_stat = -2*(RecLL_DL-RecLL_DTL),
           LR_pvalue = pchisq(LR_stat, df_DTL-df_DL, lower.tail = FALSE),
           # LR_test = -2*(RecLL_DTL-RecLL_DL) > qchisq(alpha, df_DTL-df_DL, lower.tail = FALSE)
           # AIC_test = diff < -(df_DTL-df_DL),
           LR_test = LR_pvalue < opt$alpha)

write_delim(df, opt$output, delim = "\t")

hgt_ids <- pull(filter(df, LR_test), "gene")

if (length(hgt_ids)>0){
    writeLines(paste0(opt$reco, "/", hgt_ids, "_reconciliated.xml"), opt$tk)
} else {
    writeLines("", opt$tk)
}
