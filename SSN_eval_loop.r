#Install packages
packs<-c("igraph")
lapply(packs, require, character.only=T)

#Read in all vs. all BLAST results in tabular format
allvall<-read.table("data/1263_ANL_all_v_all.tsv",stringsAsFactors=F)
head(allvall)

#Convert to a SSN
shrt<-allvall %>%
  dplyr::select(prot1 = 1,
                prot2 = 2,
                eval = 11)

#First remove all identical protein comparisons
noprs<-shrt[!shrt$prot1 == shrt$prot2,]
noprs<-noprs[order(noprs$eval),]

#Example coloring (can be changed based on protein family)
meta<-data.frame(fam = noprs$prot1,stringsAsFactors=F)
colors <- c("red", "dodgerblue") # more colors can be added
metu<-unique(meta)
metu$color <- rep(colors, nrow(metu)/length(colors)) # a

#Set a similarity threshold 
colors<-unique(metu$color)
proteins<-unique(metu$fam)

for(i in seq(from=0, to=100,by=10)){
  thresh<-as.numeric(paste0("1.00e-",i))
  net<-noprs[noprs$eval<(thresh),]
  dim(net)
  
  g<-simplify(graph_from_data_frame(net,vertices=metu,directed = FALSE))
  
  V(g)$color<-V(g)$color
  head(V(g)$color)
  gr<-delete_vertices((g),degree(g)==0)
  l<-layout_components(gr)
  
  png(file = paste0("output/SSN_trim_",thresh,".png"))
  par(mar=c(0.01,0.01,0.01,0.01))
  plot(gr,vertex.label=NA, 
       layout=l,edge.color="gray50",vertex.size=2)
  dev.off()
}

