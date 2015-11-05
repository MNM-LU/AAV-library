retrieveFASTAQID <- function(sequences,
                           PE=FALSE)
{
if (PE) {
sequencesID <- as.character(ShortRead::id(sequences))
#sequencesID <- substr(ShortRead::id(sequences),1,width(ShortRead::id(sequences))) #Use for PE reads
sequencesID <- unlist(lapply(strsplit(sequencesID, " "),"[",1)) #Use for PE reads
} else {
  sequencesID <- as.character(ShortRead::id(sequences))
}

sequencesID
}