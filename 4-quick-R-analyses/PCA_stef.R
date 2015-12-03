library("adegenet")
myFile <- read.structure("Stef_r7_1022-short.stru", onerowperind=FALSE, n.ind=49, n.loc=1022, col.lab=1)
myFile



*********************************************************************
################PCA with my data......

scaling <- scaleGen(myFile, missing="mean", scale=FALSE)
pca1 <- dudi.pca(scaling, center=FALSE, scale=FALSE) ##displays eigenvalues barplot >>>TELL NUMBER OF RETAINED AXES
pca1


##if you want individual labels in graph
s.label(pca1$li)
###this first graph is useful for identifying individua samples that fall within the origin, indicating it's an issue of too much missing data

		
###Plotting with colors		
myCol <-c("orange2", "darkblue", "red")
s.class(pca1$li, fac=pop(myFile),
		col=myCol, axesel=FALSE, cstar=0)



##############################
*********************************************************************
*********************************************************************
##find.clusters transforms the data using PCA, then it runs k-means algorithm with increasing values of K, and computes summary statistics (default BIC). 
grp <- find.clusters(myFile, max.n.clust=40) #the maximum number of clusters here is k=40
##need to pick number of PCs and number of clusters to retain
grp$size
table(pop(myFile), grp$grp)##how well have groups been retrieved by clustering
table.value(table(pop(myFile), grp$grp), col.lab=paste("inf", 1:6), row.lab=paste("ori", 1:4))
###how many clusters are actually in the data?


********************************************
********************************************
#######DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS ###
dapc <- dapc(myFile)
##don't keep that many PCs, since it may result in instability of membership probabilities. Retain as many as possible without sacrificing much information
dapc
summary(dapc)

scatter(dapc, posi.da="bottomright", bg="white", pch=17:20)


myCol <-c("orange2","darkred", "darkblue")
scatter(dapc, posi.da="bottomright", bg="white",cstar=0, col=myCol, scree.pca=TRUE)



###COMPUTE CONTRIBUTIONS OF ALLELES IN DAPC
contrib <- loadingplot(dapc$var.contr, axis=2, thres=.07, lab.jitter=1)



####ASSIGNMENT TO GROUPS AND CLUSTERING

assignplot(dapc)


********************************************
********************************************
#####GENERATING A CLUSTER GRAPH#####

X <- tab(myFile, freq=TRUE, NA.method="mean") ###this is in the tutorial online, but doesn't work for me
D <- dist(myFile)

h1 <- hclust(D, method="complete")
plot(h1, labels=TRUE)


