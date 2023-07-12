suppressMessages(library(optparse))

option_list = list(
    make_option(c("-d", "--dl"), type="character", default=NULL,
                help="generax DL rates file", dest = "dl_file"),
    make_option(c("-t", "--dtl"), type="character", default=NULL,
                help="generax DTL rates file", dest = "dtl_file"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="output filename", dest = "output")
    # make_option(c("-a","--approximate"), type = "logical", default=FALSE, action = "store_true",
    #             help="if active, reduce matrix by approximate method", dest = "approx")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressMessages(library(tidyverse))

theme_set(theme_bw(base_family = "Helvetica"))


df_DL <- 2
df_DTL <- 3
alpha <- .05

dl <- read_delim(opt$dl_file, show_col_types = F)
dtl <- read_delim(opt$dtl_file, show_col_types = F)


df <- left_join(select(dl, gene, RecLL, D_rate, L_rate),
                select(dtl, gene, RecLL, contains("rate")), 
                by="gene", suffix = c("_DL", "_DTL")) %>% 
    mutate(diff = RecLL_DTL-RecLL_DL,
           AIC_test = diff < -1,
           LR_test = 2*(RecLL_DL-RecLL_DTL) > qchisq(alpha, df_DTL-df_DL, lower.tail = FALSE))


print(table(df$LR_test, df$AIC_test))

# pivot_longer(df, cols = contains("rate_D")) %>% 
#     separate(name, c("Event","x","model"), "_") %>% 
#     # pivot_wider(names_from = Event, values_from = value) %>% 
#     ggplot(aes(model, value)) +
#     facet_grid(LR_test~Event) +
#     geom_line(aes(group = gene), alpha=.1, linewidth=.2) +
#     geom_boxplot()

hgt_ids <- pull(filter(df, LR_test), "gene")

if (length(hgt_ids)>0){
    writeLines(hgt_ids, opt$output)
} else {
    writeLines("", opt$output)
}
