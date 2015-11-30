in.names.all <- list.files("output", pattern="*.rds", full.names=TRUE)
in.names.all <- in.names.all[grep("1000x",in.names.all)]
in.names <- list(in.names.all[grep("_SN-",in.names.all)], 
                 in.names.all[grep("_Th-",in.names.all)],
                 in.names.all[grep("_Ctx-",in.names.all)],
                 in.names.all[grep("_Str-",in.names.all)])
for (i in 1:as.integer(length(in.names))) {
  #i <- 2
         out.range <- lapply(in.names[[i]], readRDS)
         out.range <- do.call(GAlignmentsList,unlist(out.range))
         out.range <- cbind(unlist(out.range))[[1]]
         out.names <- gsub("output/found.","",in.names[[i]])
         out.names <- gsub(".rds","",out.names)
         out.name <- paste("output/found.",paste(out.names, collapse="_"),".rds", sep="")
         saveRDS(out.range, file=out.name)
         #unlink(in.names[[i]], recursive = FALSE, force = FALSE) #remove non-merged files
         }
