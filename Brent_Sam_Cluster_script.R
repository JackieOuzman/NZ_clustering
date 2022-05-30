
## This is a script from Brent Robs phD student and how he does Clustering

## he said the NbClusters package is quite good.


#install.packages("cluster.datasets")
#install.packages("NbClust")
#install.packages("clValid")
#install.packages("ggfortify")
#install.packages("clustree") #won't install
#install.packages("dendextend")
#install.packages("factoextra")
#install.packages("ggiraphExtra")
#install.packages("kableExtra")

#library(clustree)#won't install

library(stats)
library(tidyr)
library(ggplot2)
library(ggplot)
library(tidyverse)
library(dplyr)
library(magrittr)
library(cluster)
library(cluster.datasets)
library(cowplot)
library(corrplot)
library(NbClust)
library(clValid)
library(ggfortify)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(GGally)
library(ggiraphExtra)
library(knitr)
library(kableExtra)

qual<-read.csv("insert_filename_here.csv")

qual<-as_tibble(qual)

glimpse(qual)

summary(qual) %>% kable() %>% kable_styling()

qual %>%
  gather(attributes, value, 1:6) %>%
  ggplot(aes(x=value)) +
  geom_histogram(fill = "lightblue2", color = "black") +
  facet_wrap(~attributes, scales = "free_x") +
  labs(x = "Value", y = "Frequency")

corrplot(cor(qual), type = "upper", method = "number", tl.cex = 0.9)

qual_scaled<- scale(qual)
res.pca<-PCA(qual_scaled)

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0,50))

var<-get_pca_var(res.pca)

fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE #Avoid text overlapping
) + theme_minimal() + ggtitle("Variables - PCA")


km2<-kmeans(qual_scaled, centers = 2, nstart = 30)

kmean_calc<- function(df, ...) {
  kmeans(df, scaled = ..., nstart = 30)
}

km2<-kmean_calc(qual_scaled, 2)
km3<-kmean_calc(qual_scaled, 3)
km4<-kmean_calc(qual_scaled, 4)
km5<-kmean_calc(qual_scaled, 5)
km6<-kmean_calc(qual_scaled, 6)
km7<-kmean_calc(qual_scaled, 7)
km8<-kmean_calc(qual_scaled, 8)
km9<-kmean_calc(qual_scaled, 9)
km10<-kmean_calc(qual_scaled, 10)

p1<-fviz_cluster(km2, data = qual_scaled, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 2")
p2<-fviz_cluster(km3, data = qual_scaled, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 3")
p3<-fviz_cluster(km4, data = qual_scaled, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 4")
p4<-fviz_cluster(km5, data = qual_scaled, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 5")
p5<-fviz_cluster(km6, data = qual_scaled, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 6")
p6<-fviz_cluster(km7, data = qual_scaled, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 7")

plot_grid(p1, p2, p3, p4, p5, p6, labels = c("k2", "k3", "k4", "k5", "k6", "k7"))

set.seed(31)

fviz_nbclust(qual_scaled, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("elbow")

gap_stat<-clusGap(qual_scaled, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: gap statistic")

fviz_nbclust(qual_scaled, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("Silhouette")

ssc<-data.frame(
  kmeans = c(2,3,4,5,6,7,8),
  within_ss = c(mean(km2$withinss), mean(km3$withinss), mean(km4$withinss), mean(km5$withinss), mean(km6$withinss), 
                mean(km7$withinss),mean(km8$withinss)), 
  between_ss = c(km2$betweenss, km3$betweenss, km4$betweenss, km5$betweenss, km6$betweenss, km7$betweenss, km8$betweenss)
)

ssc %<>% gather(., key = "measurement", value = value, -kmeans)

ssc %>% ggplot(., aes(x=kmeans, y=log10(value), fill = measurement)) +
  geom_bar(stat = "identity", position = "dodge") + ggtitle("Cluster Model Comparison") +
  xlab("Number of Clusters") + ylab("Log10 Total Sum of Squares") +
  scale_x_discrete(name = "Number of Clusters", limits = c("0", "2", "3", "4", "5", "6", "7", "8"))

res.nbclust<-NbClust(qual_scaled, distance = "euclidean",
                     min.nc = 2, max.nc = 9,
                     method = "complete", index = "all")

factoextra::fviz_nbclust(res.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")




ss_kmeans<- t(sapply(2:14,
                     FUN = function(k)
                       kmeans(x = qual_scaled,
                              centers = k,
                              nstart = 20,
                              iter.max = 25)[c('tot.withinss', 'betweenss')]))
plot(2:14, unlist(ss_kmeans[,1]), xlab = 'Clusters', ylab = 'Within Cluster SSE')

tot.ss <- sum(apply(qual_scaled, 2, var)) * (nrow(qual_scaled) - 1)
var_explained<- unlist(ss_kmeans[,2]) / tot.ss
plot(2:14, var_explained, xlab = 'Clusters', ylab = '% of Total Variation Explained')

kmeansAIC<- function(fit){
  m = ncol(fit$centers)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

aic_k<- sapply(2:14, FUN = 
                 function(k)
                   kmeansAIC(kmeans(qual_scaled, centers = k, nstart = 20, iter.max = 25)))
plot(2:14, aic_k, xlab = 'Clusters', ylab = 'AIC from kmeans')


kmeansBIC <- function(fit){
  m = ncol(fit$centers) 
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + log(n) * m * k) 
}

bic_k <- sapply(2:14, FUN = 
                  function(k) 
                    kmeansBIC(kmeans(qual_scaled, centers = k, nstart = 20, iter.max = 25)))

plot(2:14, aic_k, xlab = 'Clusters', ylab = 'BIC from kmeans')

kmeans2 <- function(data, center_range, iter.max, nstart, plot = TRUE){
  
  #fit kmeans for each center
  all_kmeans <- lapply(center_range, 
                       FUN = function(k) 
                         kmeans(data, center = k, iter.max = iter.max, nstart = nstart))
  
  #extract AIC from each
  all_aic <- sapply(all_kmeans, kmeansAIC)
  #extract BIC from each
  all_bic <- sapply(all_kmeans, kmeansBIC)
  #extract tot.withinss
  all_wss <- sapply(all_kmeans, FUN = function(fit) fit$tot.withinss)
  #extract between ss
  btwn_ss <- sapply(all_kmeans, FUN = function(fit) fit$betweenss)
  #extract totall sum of squares
  tot_ss <- all_kmeans[[1]]$totss