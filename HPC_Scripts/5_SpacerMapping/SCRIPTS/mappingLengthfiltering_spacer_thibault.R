#!/usr/bin/env Rscript
library("optparse")
library("tidyverse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input sam", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
RMK202_rawIllumina_bwamem_sorted_mapped_mapQ_30 <- read.delim2(opt$input, header=FALSE, quote="")
RMK202_rawIllumina_bwamem_sorted_mapped_mapQ_30_out <- RMK202_rawIllumina_bwamem_sorted_mapped_mapQ_30 %>% mutate(add = str_extract(V13, "[0-9]+")) %>% mutate_at('add',as.numeric) %>% filter(.,add>42) %>% select(-"add")
write.table(RMK202_rawIllumina_bwamem_sorted_mapped_mapQ_30_out, opt$out, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = FALSE,quote = FALSE)
