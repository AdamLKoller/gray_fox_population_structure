rm(list=ls())
library(tess3r)
library(maps)
library(Rgraphviz)

do_something <- function(data_path, out_path, threads, myparam) {
    print('hello world')
}

do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config[["myparam"]])