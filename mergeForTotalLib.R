in.names.all <- list.files("output", pattern="*.rds", full.names=TRUE)
         out.range <- lapply(in.names.all, readRDS)
         out.range <- do.call(GAlignmentsList,unlist(out.range))
         out.range <- cbind(unlist(out.range))[[1]]
         saveRDS(out.range, file="output/total.infectiveLib.rds")
