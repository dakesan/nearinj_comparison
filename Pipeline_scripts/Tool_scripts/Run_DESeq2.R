#Run_DeSeq2

#* DESeq2がインストールされていればlibrary(DESeq2)
#* なければインストール
deps = c("DESeq2")
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("DESeq2")
  }
  library(dep, character.only = TRUE)
}

library(DESeq2)

#* 引数が適切に渡されているかを確認している
#! Rscript $TOOL_DIR/Run_DESeq2.R $ASV_table_Path $Groupings_Path $out_file_deseq
args <- commandArgs(trailingOnly = TRUE)
#test if there is an argument supply
if (length(args) <= 2) {
  stop("At least three arguments must be supplied", call.=FALSE)
}

#* ASV table path
con <- file(args[1])
file_1_line1 <- readLines(con,n=1)
close(con)

#* constructed from biome fileという余計な一行が先頭にあるかを確認している
if(grepl("Constructed from biom file", file_1_line1)){
  ASV_table <- read.table(args[1], sep="\t", skip=1, header=T, row.names = 1,
                          comment.char = "", quote="", check.names = F)
}else{
  ASV_table <- read.table(args[1], sep="\t", header=T, row.names = 1,
                          comment.char = "", quote="", check.names = F)
}

#* "サンプル名, グループ名"というテーブルを別途用意、読み込んでいる
groupings <- read.table(args[2], sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)

#number of samples
sample_num <- length(colnames(ASV_table))
grouping_num <- length(rownames(groupings))

#check if the same number of samples are being input.
if(sample_num != grouping_num){
  message("The number of samples in the ASV table and the groupings table are unequal")
  message("Will remove any samples that are not found in either the ASV table or the groupings table")
}

#check if order of samples match up.
if(identical(colnames(ASV_table), rownames(groupings))==T){
  message("Groupings and ASV table are in the same order")
}else{
  #* "groupings"というテーブルがdeseqに使うcolData
  rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
  groupings <- groupings[rows_to_keep,,drop=F]
  ASV_table <- ASV_table[,rows_to_keep]
  if(identical(colnames(ASV_table), rownames(groupings))==T){
    message("Groupings table was re-arrange to be in the same order as the ASV table")
    message("A total of ", sample_num-length(colnames(ASV_table)), " from the ASV_table")
    message("A total of ", grouping_num-length(rownames(groupings)), " from the groupings table")
  }else{
    stop("Unable to match samples between the ASV table and groupings table")
  }
}

colnames(groupings)[1] <- "Groupings"
#Run Deseq2

#* メイン関数
dds <- DESeq2::DESeqDataSetFromMatrix(countData = ASV_table,
                                      colData=groupings,
                                      design = ~ Groupings)
dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")

res <- results(dds_res, tidy=T, format="DataFrame")

rownames(res) <- res$row
res <- res[,-1]

write.table(res, file=args[3], quote=FALSE, sep="\t", col.names = NA)

message("Results written to ", args[3])
