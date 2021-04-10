#!/usr/bin/env Rscript

if (!require("pacman")) install.packages("pacman")

pacman::p_load(ggplot2, cowplot, svglite, ggridges, RColorBrewer, optparse)

pacman::p_load_gh("laduplessis/bdskytools")

option_list = list(
    make_option(c("-d", "--dates"), type="character", default="",
                help="input dataset file name", metavar="character"),
    make_option(c("-r", "--date_pair"), type="character", default="",
                help="input dataset file name", metavar="character"),
    make_option(c("-l", "--logpath"), type="character", default=getwd(),
                help="input log file path [default= %default]", metavar="character"),
    make_option(c("-p", "--plotpath"), type="character", default=getwd(),
                help="output plot path [default= %default]", metavar="character"),
    make_option(c("-n", "--plotname"), type="character", default="ridge",
                help="output plot name[default= %default]", metavar="character"),
    make_option(c("-c", "--ci"), type="numeric", default=0.95,
                help="hpd or credible interval [default= %default]", metavar="numeric"),
    make_option(c("-b", "--burnin"), type="numeric", default=0.1,
                help="burnin to discard [default= %default]", metavar="numeric"),
    make_option(c("-f", "--dateformat"), type="character", default="min-max",
                help="date format of --dates [default= %default]", metavar="numeric"),
    make_option(c("-e", "--logext"), type="character", default=".log",
                help="log file extension of beast outputs [default= %default]", metavar="numeric"),
    make_option(c("-g", "--namesplit"), type="character", default=".log",
                help="file name split taking first of split for id with --dates [default= %default]", metavar="numeric"),
    make_option(c("-s", "--slice"), action="store_true", default=FALSE,
                help="sample proportion is slices pre-sample / sample period [default= %default]"),
    make_option(c("-i", "--infectious"), type="numeric", default=NULL,
                help="set a limit on the infectious period (can have long tail sometimes) [default= %default]")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


source("C:/Users/esteinig/Desktop/Backup/BEAST_LINEAGES/DATA/R/beastling.R")


# Required catches


extract_prefix <- function(fname){
    s <- stringr::str_split(fname, opt$namesplit)[[1]]
    label <- s[1]
}

if (opt$date_pair != ""){

    dates <- stringr::str_split(opt$date_pair, ":")[[1]]

    min_date <- dates[1]
    max_date <- dates[2]

    logfiles = list.files(path=opt$logpath, pattern=paste0("*", opt$logext))
    file_names <- sapply(logfiles, basename)
    fnames <- sapply(file_names, tools::file_path_sans_ext)
    prefixes <- sapply(fnames, extract_prefix)

    dates <- data.frame(
        first_sample_year=rep(as.numeric(min_date), length(prefixes)),
        last_sample_year=rep(as.numeric(max_date), length(prefixes))
    )

    rownames(dates) <- prefixes

} else {
    if (opt$dateformat == "min-max"){
        dates <- read.table(file = opt$dates, row.names = 1, sep="\t")
        colnames(dates) <- c('first_sample_year', 'last_sample_year')
    } else {
        stop("Date format must be min-max for now")
    }
}


print(dates)

posteriors <- read_log_files(
    path=opt$logpath, logpattern=paste0("*", opt$logext), data=dates, sample_slice=opt$slice,
    burnin=opt$burnin, hpd=opt$ci, extract_label_function=extract_prefix
)

plot_posterior_multiple(
    df_posterior=posteriors$posterior,
    plot_path=opt$plotpath,
    plot_name=opt$plotname,
    infectious_limit=opt$infectious
)


