# read scanned Pigs data and check against version in package

PigsWide <- read.table("~/Documents/temp/pigs-wide.txt")
PigsWide[PigsWide > 100] <- PigsWide[PigsWide > 100]/10
PigsWide
PigsLong <- reshape(PigsWide, varying=paste0("V", 1:9), direction="long", idvar="id", v.names="weight")
PigsLong <- PigsLong[order(PigsLong$id), ]
PigsLong <- PigsLong[, c("id", "time", "weight")]
names(PigsLong) <-  c("id", "week", "weight")
rownames(PigsLong) <- 1:nrow(PigsLong)
head(PigsLong, 18)
all.equal(Pigs, PigsLong)
