#
# External arguments
#

options(warn=-1)
suppressMessages(library(igraph))

# args = '/Users/mviana/Desktop/MitoGraph/data/MitoGraphTest001'
args = commandArgs(trailingOnly=TRUE)

FileName <- args[1]

#
# Loading graph
#

G <- read.table(paste(FileName,'.gnet',sep=''), skip=1, col.names=c('Source','Target','Length'))
G <- graph.data.frame(as.data.frame(G), directed=F)

#
# Loading volume information
#

Vol <- data.frame(read.table(paste(FileName,'.cc',sep=''), skip=1))
names(Vol) <- c('node','cc','vol_cc')

ids <- as.numeric(V(G)$name)
V(G)$cc_vol <- Vol$vol_cc[match(ids,Vol$node)]

#
# Connected component analysis
#

List <- decompose(G)

Table <- NULL
for (g in List) {
  Table <- rbind(Table,data.frame(vcount(g), ecount(g), sum(E(g)$Length), max(V(g)$cc_vol)))
}

names(Table) <- c('nodes','edges','length_(um)','vol_from_img_(um3)')

Table <- Table[order(Table[,4], decreasing=T),]

#
# Loading global measurements
#

Global <- data.frame(read.table(paste(FileName,'.mitograph',sep=''), skip=1))

names(Global) <- c('vol_from_voxels_(um)', 'avg_width_(um)', 'std_width_(um)', 'total_length_(um)', 'vol_from_length_(um3)')

Global$nodes <- vcount(G)

Global$edges <- ecount(G)

Global$components <- length(List)

output_file <- paste(FileName,'.mitograph',sep='')

cat('::MitoGraph\n',file=output_file, append=F)
cat(paste(Sys.time(),'\n'),file=output_file, append=T)
cat('\n',file=output_file, append=T)
cat('GLOBAL STATISTICS:\n',file=output_file, append=T)
cat('\n',file=output_file, append=T)
write.table(Global,output_file,
          append=T,
          row.names=F,
          col.names=T,
          sep="\t",
          quote=F)
cat('\n',file=output_file, append=T)
cat('COMPONENTS STATISTICS:\n',file=output_file, append=T)
cat('\n',file=output_file, append=T)
write.table(Table,output_file,
          append=T,
          row.names=F,
          col.names=T,
          sep="\t",
          quote=F)
          
