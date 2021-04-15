library(data.table)
library(Rtsne)
library(ggplot2)
library(caret)
library(ggplot2)
library(ClusterR)

data<-fread("cladocopium_msats_use_tsne.csv")

# we are running 100 t-sne analyses to demonstrate that despite the random nature of the spacial embeddings
# the relevent conclusions about the groupings of samples are always supported by the analysis. 

for (i in 1:100){
seed<-i

set.seed(seed)

# load in data data<-fread("~/Downloads/cladocopium_msats_use_tsne copy.csv")


m_data<-melt(data,id.vars=c("sample","site","spp"))
m_data<-m_data[!m_data$value==0]
new_data<-dcast(m_data,sample+site+spp ~ variable + value,length,value.var = "value" )

spp<-new_data$spp
site<-new_data$site
sample<-new_data$sample
#This will then remove those columns so you dont feed them into the numeric operations that happen next:
  new_data$spp<-NULL
new_data$site<-NULL
new_data$sample<-NULL



# do a pca
pca<-prcomp(new_data)

# look at the percent variance explained by each pca
screeplot(pca)

# look at the rotation of the variables on the PCs
pca

# see the values of the scree plot in a table 
summary(pca)

# see a biplot of the first 2 PCs
#biplot(pca)

# use the unclass() function to get the data in PCA space
pca_dt<-data.table(unclass(pca)$x)

# add back the party to prove to ourselves that this works
pca_dt$spp<-spp

# see a plot with the spp data 
ggplot(pca_dt,aes(x=PC1,y=PC2,col=spp))+geom_point()




# run t-sne on the PCAs, note that if you already have PCAs you need to set pca=F or it will run a pca again. 
# pca is built into Rtsne, ive run it seperatly for you to see the internal steps

## updated Aug 2020
# I increased max_iter and lowered theta to insure a stable embedding was reached.
# Lower theta and high max_iter should produce better results at the cost of computation 
# speed. I did not have an issue running these settings on a macbook pro 

tsne<-Rtsne(new_data,pca = T,theta=0.1,perplexity = 10,max_iter = 10000)

# grab out the coordinates
tsne_dt<-data.table(tsne$Y)

# add back in spp and site so we can see what the analysis did with them
tsne_dt$spp<-spp
tsne_dt$site<-site
tsne_dt$sample<-sample

# plotting in black to see if there are groups. 
ggplot(tsne_dt,aes(x=V1,y=V2))+geom_point()


# use a gaussian mixture model to find optimal k and then get probability of membership for each row to each group

# this fits a gmm to the data for all k=1 to k= max_clusters, we then look for a major change in likelihood between k values
k_bic<-Optimal_Clusters_GMM(tsne_dt[,.(V1,V2)],max_clusters = 10,criterion = "BIC",dist_mode = 'eucl_dist')

# now we will look at the change in model fit between successive k values
delta_k<-c(NA,k_bic[-1] - k_bic[-length(k_bic)])

# I'm going to make a plot so you can see the values, this part isnt necessary
del_k_tab<-data.table(delta_k=delta_k,k=1:length(delta_k))

# plot 
ggplot(del_k_tab,aes(x=k,y=-delta_k))+geom_point()+geom_line()+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_text(aes(label=k),hjust=0, vjust=-1)



##########################
# function for picking k #
##########################

getOptK <- function(x){
  x <- abs(x)
  max_delta_k_pos <- which.max(x)
  max_delta_k <- max(na.omit(x))
  n2eval<-(length(x) - max_delta_k_pos) - 2
  for(i in max_delta_k_pos:(max_delta_k_pos+n2eval)){
    if(x[max_delta_k_pos + 1]/max_delta_k < 0.15 & x[max_delta_k_pos + 2]/max_delta_k < 0.15 & x[max_delta_k_pos + 3]/max_delta_k < 0.15){
    }else{
      max_delta_k_pos <- max_delta_k_pos + 1
    }
  }
  max_delta_k_pos
}


# You may visualy inspect the plot to pick the optimal k, I have writen a function that expresses the logic that I use
#opt_k<-getOptK(delta_k)

#opt k was 3 but based on the other data collected in this study we are exploring k=2 to 
#see if the subdivision seen at k=3 consolidates at k=2 such that the clusters match the
#ITS haplotypes. We discus the significance of k=3 as potential 3rd species or genetic 
#subdivision within a species. 

opt_k<-2

## updated Aug 2020
# I increased km_iter and em_iter to insure that the groupings were stable.
# The clusterR package says that only 10 iterations are often neccessary but I observed
# very poor groupings that did not seem consistant with the spacial distribution of
# points. Increasing these parmeters resulted in clusters that followed logical visual groupings.
# Higher iterations should always result in better groupings, the only trade of is computation time,
# this was not an issue for my macbook pro.

# now we run the model with our chosen k value
gmm_data<-GMM(tsne_dt[,.(V1,V2)],opt_k, km_iter = 4000,em_iter=4000,dist_mode = 'eucl_dist')

# the model gives a log-likelihood for each datapoint's membership to each cluster, me need to convert this 
# log-likelihood into a probability

l_clust<-gmm_data$Log_likelihood^10

l_clust<-data.table(l_clust)

net_lh<-apply(l_clust,1,FUN=function(x){sum(1/x)})

cluster_prob<-1/l_clust/net_lh

# we can now plot to see what cluster 1 looks like

tsne_dt$Cluster_1_prob<-cluster_prob$V1

tsne_dt$spp_num<-0
tsne_dt[spp=="c1d"]$spp_num<-15
tsne_dt[spp=="c1c"]$spp_num<-17

ggplot(tsne_dt,aes(x=V1,y=V2,col=Cluster_1_prob,shape=spp))+geom_point()+scale_shape_manual(values=c(15,17))
filename<-paste0("./probability/clusterprobability_",seed,".pdf")
ggsave(filename)

ggplot(tsne_dt,aes(x=V1,y=V2,col=spp))+geom_point()
filename<-paste0("./spp/clusterspp_",seed,".pdf")
ggsave(filename)
##plot by species
#ggplot(tsne_dt,aes(x=V1,y=V2,col=spp))+geom_point()

##plot by site
ggplot(tsne_dt,aes(x=V1,y=V2,col=site))+geom_point()
filename<-paste0("./site/clustersite_",seed,".pdf")
ggsave(filename)
#cluster 2

#tsne_dt$Cluster_2_prob<-cluster_prob$V2

#ggplot(tsne_dt,aes(x=V1,y=V2,col=Cluster_2_prob))+geom_point()

#cluster 3 (for k=3)
 
#tsne_dt$Cluster_3_prob<-cluster_prob$V3

#ggplot(tsne_dt,aes(x=V1,y=V2,col=Cluster_3_prob))+geom_point()

#plot by species
#ggplot(tsne_dt,aes(x=V1,y=V2,col=spp))+geom_point()

#plot by site
#ggplot(tsne_dt,aes(x=V1,y=V2,col=site))+geom_point()



}

