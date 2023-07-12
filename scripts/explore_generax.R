suppressMessages(library(optparse))

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="generax result directory", dest = "input"),
    make_option(c("-c", "--code"), type="character", default=NULL,
                help="code phylomedb", dest = "code"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="output prefix", dest = "output")
    # make_option(c("-a","--approximate"), type = "logical", default=FALSE, action = "store_true",
    #             help="if active, reduce matrix by approximate method", dest = "approx")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# opt <- NULL
# opt$input <- "output/reconciliation/phylome_0845/generax_DTL"
# opt$code <- "0845"
# opt$output <- "test/test_results.png"

suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(ggtree))

col_events <- c("D" = "#F2B705", "S" = "#04ADBF", "SL" = "#F2380F", "T" = "#146551")
col_rates <- c("D" = "#F2B705", "L" = "#F2380F", "T" = "#146551")
cnt_pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
is_DTL <- FALSE


dir.create(dirname(opt$output), showWarnings = FALSE, recursive = T)
# prefix <- paste0(opt$output,"/generax_", opt$code, "_", opt$model)
prefix <- gsub("_rates.tsv", "", opt$output)


grax_files <- list.files(opt$input, full.names = T)

# whole reconciliation stats
title_sptree <- paste0(readLines(grax_files[grep("stats.txt", grax_files)], warn = F), collapse = " ")

# count of events per species
spcounts <- read_delim(grax_files[grep("per_species_event_counts.txt", grax_files)], 
           col_names = c("species", "S", "SL", "D", "T", "TL"), skip = 1, show_col_types = FALSE) %>% 
    mutate_at(vars(-species), ~as.numeric(gsub(".*=","", .)))

# missing data per species
sp_presence <- read_delim(grax_files[grep("perSpeciesCoverage.txt", grax_files)], skip=1, 
           col_names = c("species", "family_coverage"), show_col_types = FALSE) %>% 
    mutate(species=gsub(":", "", species)) %>% 
    left_join(read_delim(grax_files[grep("fractionMissing.txt", grax_files)], 
                         col_names = c("species", "missing"), show_col_types = FALSE))

# missing data influence in the dataset
plot_presence <- ggplot(sp_presence, aes(missing, family_coverage)) + 
    geom_point() +
    ggrepel::geom_text_repel(aes(label=species)) +
    theme_bw()

# results mapped to species tree
get_internal_id <- function(tree, name) {
    nodeid <- which(tree$node.label==name) + ape::Ntip(tree)
    if (name %in% tree$tip.label){
        nodeid <- name
    }
    return(nodeid)
}

sptree <- read.tree(paste0(grax_files[grep("species_trees", grax_files)], "/inferred_species_tree.newick"))
sptree$edge.length <- NULL

spcounts$node <- sapply(spcounts$species, function(x) get_internal_id(sptree, x), simplify = T)

sptreeplot <- ggtree(sptree) + 
    geom_tiplab() +
    ggtitle(title_sptree) + 
    # theme(legend.position = "right") + 
    xlim(c(0,20))

pies <- nodepie(spcounts[grepl("node_", spcounts$species),], cols=2:5, alpha=.8)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = col_events))

sptreeplot <- inset(sptreeplot, pies, width = .1, height = .1)

# for each gene you have
# results: stats with DTL rates and corrected tree (right now useless but maybe we can recall orthologs)
res_dir <- grax_files[grep("results$", grax_files)]
res_ids <- list.files(res_dir)

fam_sts_nms <- c("LibpllLL", "RecLL", "D_rate", "L_rate", "T_rate")

df_rates <- tibble()

print("starting loop")
for (gene in res_ids) {
    family_stat <- paste(read_lines(paste0(res_dir, "/", gene, "/stats.txt")), collapse = " ")
    family_stat <- c(str_split(gsub("Reconciliation rates = ", "", str_trim(family_stat)), pattern = " ", simplify = T))
    family_stat <- as.numeric(family_stat)
    if (length(family_stat)==4) { 
        family_stat <- c(family_stat, NA)
    }
    names(family_stat) <- fam_sts_nms
    df_rates <- bind_rows(df_rates, family_stat)
}
df_rates$gene <- res_ids
print("finished loop")

# reconciliations: eventcounts, trasnfers-> do a tree with connections based on this!!!! 
# specieseventcounts, reconciliated trees, orthogroups are useless

rec_dir <- grax_files[grep("reconciliations", grax_files)]
ecounts_files <- list.files(rec_dir, pattern = "eventCounts.txt", full.names = T)

df_rates <- read_delim(ecounts_files, id="gene", delim = ":", 
                       col_names = c("event", "n"),
                       show_col_types = FALSE) %>% 
    mutate(gene=gsub("_eventCounts.txt", "", basename(gene))) %>% 
    pivot_wider(id_cols = gene, names_from = event, values_from = n) %>% 
    left_join(df_rates)


rates_plot <- select(df_rates, c(gene, contains("rate"))) %>% 
    pivot_longer(!gene) %>% 
    mutate(name=gsub("_rate", "", name)) %>% 
    ggplot(aes(value, name, fill=name)) + 
    ggdist::stat_slabinterval() +
    scale_fill_manual(values = col_rates) +
    theme_bw()

# REACTIVATE ONCE THE EVAL MODE DOES WHAT IT SHOULD DO
# ll_plot <- ggplot(df_rates, aes(LibpllLL, RecLL)) + 
#     geom_point() +
#     theme_bw()


# if DL mode this is useless
transfer_files <- list.files(rec_dir, pattern = "_transfers.txt", full.names = T) 
transfer_files <- transfer_files[which(file.size(transfer_files)>0)]

if (length(transfer_files)>0){
    is_DTL <- TRUE
    transfers <- map_df(transfer_files, ~read_delim(., delim = " ", 
                                                    col_names = c("from", "to"), show_col_types = FALSE)) %>% 
        # read_delim(transfer_files, delim = " ", col_names = c("from", "to"),
        #                     show_col_types = FALSE) %>% 
        group_by(from, to) %>% 
        count() %>% 
        filter(!is.na(from))

    transfer_plot <- ggtree(sptree, layout = "circular") +
        geom_tiplab() +
        geom_taxalink(data = filter(transfers, n>quantile(transfers$n, .75)), 
                      aes(taxa1=from, taxa2=to, color=n, alpha=n)) + 
        scale_color_gradientn(colours = cnt_pal) + 
        scale_alpha(range = c(0, 1))

    events_plot <- ggplot(df_rates, aes(D/Leaf, SL/Leaf, color=T/Leaf)) + 
        geom_point() + 
        scale_color_gradientn(colours = cnt_pal) +
        theme_bw()
} else {
    events_plot <- ggplot(df_rates, aes(D/Leaf, SL/Leaf)) + 
        geom_point(alpha=.6, size=.5) + 
        theme_bw()
}


## FINAL PLOT
if (is_DTL) {
    dash <- (plot_presence/rates_plot/events_plot) | (sptreeplot/transfer_plot)
    write_delim(transfers, file = paste0(prefix, "_transfers.tsv"), delim = "\t")
} else {
    dash <- (plot_presence/rates_plot/events_plot) | (sptreeplot)
}


# final RESULTS
ggsave(paste0(prefix, "_results.png"), dash, dpi = 300, width = 15, height = 15)
write_delim(spcounts, file = paste0(prefix, "_spcounts.tsv"), delim = "\t")
write_delim(df_rates, file = opt$output, delim = "\t")

# TODO
# maybe add info so that species have normal names
# fix xlim mess ggtree
# save results in tsv files

# when comparing DTL vs DL plot boxplot of t_rates 