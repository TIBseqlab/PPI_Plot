#This script is prepread to plot the protein protin interaction for a given protein
#Usage: 
# Arguments:
# 1. The given protein 
# 2. The name of output directory
# 3. Width
# 4. Height
# example: Rscript ppi_plot.R POLD1 ./ 700 600

Args <- commandArgs(TRUE)
# output folder
folder = Args[2]

# Input INDEL list
gene =as.character(Args[1])

# import igraph library
library("igraph")

##get interaction from biogrid database
data<-read.delim(paste("http://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=",gene,"&taxId=9606&includeInteractors=true&includeInteractorInteractions=false&accesskey=997eae2170c0043e4a7d24987f163316",sep=""),header=FALSE)

#get the interactions
#interact<-t(apply(data.frame(from=data$V8,to=data$V9),1,sort))
interact<-data.frame(t(apply(data.frame(from=data$V8,to=data$V9),1,sort)),stringsAsFactors=FALSE)

# sort and remove redundants
int <- interact[interact[,1]!= interact[,2],]
links<-unique(int)

#colnames(links)<-c("P1","P2")
write.table(links,paste(folder,gene,"_ppi.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE)

nodes<-data.frame(unique(c(links[,1],links[, 2])))
net <- graph_from_data_frame(links, directed=FALSE, vertices=nodes)
deg <- degree(net, mode="all")
V(net)$size <- log(deg+1)*3

png(file=paste(folder,gene,"_ppi.png",sep=""),  bg="white",width=as.double(Args[3]), height=as.double(Args[4]))
par(mar = c(0, 0, 0, 0));
plot(net,layout=layout.lgl,edge.color="grey50",vertex.color="orange", vertex.frame.color="#ffffff",
     , vertex.label.color="black")
dev.off()
