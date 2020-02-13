#######################################################
#                                                     #
#     Network Analysis With Cosine Similarity         #
#                                                     #
#######################################################


#  put all therequired data under the following directory 
setwd("W:/BOB HAZEN_NETWORK ANALYSIS/NETWORK ANALYSYS/Hoover Res/Test2")

setwd("~/helps/Gerald/")

#  download your own data here
Hoover = read.table("Test2.csv", header=TRUE, sep = ',')

#  Hoover Data
Hoover.active <- Hoover[1:82, c(2:16,19,20,23:27)]      


#  load igraph
library(igraph)
library(reshape2)
library(ggplot2)
library(corrgram)
library(GGally)
library(PerformanceAnalytics)
library(psych)
library(randomForest)

Hoover.active.scaled <- scale(Hoover.active)   # sacle all the values

## examining variable distributions
summary(Hoover.active.scaled)

## boxplots of all variables
# ggplot(melt(Hoover.active), aes(variable, value)) + geom_boxplot() + facet_wrap( ~ variable, scales = "free")
ggplot(melt(Hoover.active.scaled), aes(Var2, value)) + geom_boxplot() + facet_wrap( ~ Var2, scales = "free")

## scatter plot matrix
## saving plot directly to file to avoid error in rendering due to large margins
png('Algal_Bloom_scatterplotmatrix.png',width=2500,height = 2500)
pairs.panels(Hoover.active.scaled, 
                     method = "pearson", # correlation method
                     hist.col = "#00AFBB",
                     density = TRUE,  # show density plots
                     ellipses = TRUE, # show correlation ellipses
                     cex.cor = 1
)
dev.off()


#  same Hoover.active normalized data
write.csv(Hoover.active, "Hoover.active.csv", row.names=TRUE)  # both
Hoover.active <- read.table("Hoover.active.csv", header=T, sep =',', row.names =1)


#  library cosine similarity instead of correlation
library(coop)

mat_both = cosine(t(Hoover.active))

## evaluating cos similarity distribution
library(dplyr)

## get upper triangle of similarity matrix
mat.up.tri <- mat_both[upper.tri(mat_both, diag = FALSE)]

## examine flattened upper triangle
summary(c(mat.up.tri))
boxplot(c(mat.up.tri))

#  correlation below 0.4 ignore
mat_both[mat_both<0.35]=0  # you may change the value r = 0.5, 0.6, ....


#  network analysis using igraph
network_both = graph_from_adjacency_matrix(mat_both, weighted=T, mode="undirected", diag=F)
set.seed(4)
clp <- cluster_louvain(network_both)

#  you can choose different cluster module and pick the highest cluster module that gives you the highest modularity
#  for exmaple,
clp2 <- cluster_walktrap(network_both)  
modularity(clp)
modularity(clp2)

## evaluating community detection
### assessing modularity matrices

#  add communities to network_both data
V(network_both)$community<- clp$membership

#  this is showing a membership and date
Hoover['membership'] <- clp$membership
Hoover[,c('Date','membership')]


# pam cluster incidence of Hoover by numbers
library(cluster)
res.pam  <- pam(Hoover$TotalP, 4)

## Unsupervised clustering
## evaluating cluster numbers using silhouette width

sil_width <- c()
x.dist <- daisy(as.data.frame(na.omit(Hoover)), metric = "gower")

invisible(sapply(3:9,function(i){pam_fit <- pam(x.dist,diss = TRUE,k = i); sil_width[i] <<- pam_fit$silinfo$avg.width; }))

plot(sil_width,xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(sil_width)


################################################################################
####       Create Loop for each variable
################################################################################

file_list <- c("NN","Ammonia","TotalP","OrthoP","Diatoms","G","YG",
               "BG","ZooP","DinoF")

for (i in file_list){
   library(dplyr)
   H_pam<-select(Hoover,c('Date',i))
   print(colnames(H_pam))
   print(H_pam)
   pam_m <- H_pam[order(H_pam[i],decreasing = FALSE),]
   print(pam_m)
   res.pam2  <- pam(pam_m[i], 4)
   print(res.pam2)

   #######################################
   ###  create pam number in the order ###
   #######################################
   library(dplyr)
   #H_pam <-select(Hoover, c('Date','BG'))
   #pam_m <- H_pam[order(H_pam$BG,decreasing = FALSE),]
   #res.pam2  <- pam(pam_m$BG, 4)
   
   #  add pam results
   pam_m$pam <- res.pam2$clustering
   pam_m
   
   # color palette for pam
   library(RColorBrewer)
   
   coul = brewer.pal(nlevels(as.factor(pam_m$pam)), "Set2")
   coul = c("#8DA0CB","#66C2A5", "#FC8D62","#E78AC3")
   # Map the color to cylinders
   my_color=coul[as.numeric(as.factor(pam_m$pam))]
   pam_m$my_color <- my_color
   
   S_pam <-merge(Hoover,pam_m, by=c("Date"),all=T,sort=F)
   S_pam
   
   my_color = S_pam$my_color
   # color palette for membership
   #library(RColorBrewer)
   #coul = brewer.pal(nlevels(as.factor(Hoover$membership)), "Set3")
   # Map the color to cylinders
   #my_color=coul[as.numeric(as.factor(Hoover$membership))]
   
   
   # change plot based on different layout
   l <- layout.fruchterman.reingold(network_both)
   set.seed(4)
   
   #plot(network_both, layout = l, vertex.label = Hoover$Date, mark.group=clp)
   
   set.seed(4)
   plot(network_both, edge.color = 5, ,vertex.color=my_color, layout = l, vertex.size=Hoover$TotalP/0.01, vertex.label = NA, mark.group=clp)
   legend(x=0.5, y=1.4, legend= c("Low","Medium", "High","Extreme"), #format(res.pam$medoids,digits=1), #c("Low","Mid","High"), 
          col= coul, bty = "n", pch=20 , pt.cex = 2, cex = 1, text.col="black" , horiz = F)
   text(x= 0,y = 1.3,i, col="black")
}




# change cluster and the layout
set.seed(4)
clp <- cluster_louvain(network_both)
V(network_both)$community <- clp$membership
plot(network_both, vertex.color=my_color, vertex.label = Hoover$Date, layout = l, mark.group=clp, vertex.size = Hoover$NN/1500)
#cluster_optimal(network_both)


set.seed(4)


# We can also plot the communities without relying on their built-in plot:
V(network_both)$community <- clp$membership

set.seed(4)
plot(network_both, mark.group = clp, edge.color = 5, vertex.color=NA, vertex.label=Hoover$Date)#, vertex.color=coul, vertex.label=Hoover$Date)
legend(x=0.5, y=1.4, legend= c("Low","Extreme", "High","Medium"), #format(res.pam$medoids,digits=1), #c("Low","Mid","High"), 
       col= coul , bty = "n", pch=20 , pt.cex = 2, cex = 1, text.col="black" , horiz = F)

write.csv(Hoover, "test.csv", row.names=TRUE)  # both


