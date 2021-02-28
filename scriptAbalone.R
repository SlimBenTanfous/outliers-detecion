
set.seed(41593)
##bivariate normal-distribution
x <- rnorm(1000,0,5)
y <- rnorm(1000,0,1)
par(mfrow=c(1,2))
plot(x,y)
x1 <- c(-14,-14,14,14); y1 <- c(-2.5,2.5,-2.5,2.5)
points(x1,y1,col="red",pch=19)
x <- c(x,x1); y <- c(y,y1)
data <- cbind(x,y)

abalone <- read.csv('C:\\Users\\lt32132\\Downloads\\Abalone.csv', header=TRUE, sep=",")
names(abalone) <- c("sex", "length", "diameter", "height", "weight.whole",
                    "weight.shucked", "weight.viscera", "weight.shell", "rings")
aba_out <- abalone[,2:3]
out <- which((aba_out$length<0.25 & aba_out$diameter>0.25)|
               (aba_out$length>0.6 & aba_out$diameter<0.4)|
               (aba_out$length>0.25 & aba_out$diameter<0.15))                
outliers <- rep(1,length(aba_out[,1]))
outliers[out] <- 2
plot(aba_out,col=outliers,pch=18)





# **1. Boxplot Method** #


bp <- function(X,fac){
med <-sapply(X,median) 
q25 <-sapply(X,function(x)quantile(x,prob=0.25)) 
q75 <-sapply(X,function(x)quantile(x,prob=0.75)) 
erg <- t(apply(X, 1, function(x) abs(med-x)-fac*(q75-q25)))
return(as.vector(which(rowSums(erg>0)>0)))}

bp_app <-function(X,a){
outliers <- rep(1,length(X[,1]))
outliers[bp(X,a)] <- 2
outliers<- as.factor(outliers)
levels(outliers) <- c("No Outlier","Outlier")
print(table(outliers))
if(table(outliers)[2] > 0)plot(X,col=outliers,pch=18)
return(outliers)}

#### application bivariate Normal-distribution ####
par(mfrow=c(1,2))                
print("Boxplot methode with factor 1.5 on bivariate normal-distribution")
bp1.5_1=bp_app(as.data.frame(data),1.5)
print("Boxplot methode with factor 3 on bivariate normal-distribution")
bp3_1=bp_app(as.data.frame(data),3)

### abalone Data ###
print("Boxplot methode with factor 1.5 on abalone data")                 
bp1.5_2 <-bp_app(aba_out,1.5)
print("Boxplot methode with factor 3 on abalone data")                 
bp3_2 <- bp_app(aba_out,3)


# **2. Euclid Method** #


euclid <- function(X,fac){
med <-sapply(X,median) 
erg <- t(apply(X, 1, function(x) (med-x)^2))
dist <- sqrt(rowSums(erg))
#   print(plot(dist))
return(dist > fac*median(dist))}

euclid_app <-function(X,a){
outliers <- rep(1,length(X[,1]))
outliers[euclid(X,a)] <- 2
outliers<- as.factor(outliers)
levels(outliers) <- c("No Outlier","Outlier")
print(table(outliers))
if(table(outliers)[2] > 0)plot(X,col=outliers,pch=18)
#return(outliers)
}

#### application on bivariate Normal-Distribution ####
print("Euclid methode on bivariate normal_distribution")                 
par(mfrow=c(1,2))
data_std <- apply(data,2,function(x)x/(max(x)-min(x)))
euclid_app(as.data.frame(data_std),3)

### abalone ###

aba_std <- apply(aba_out,2,function(x)x/(max(x)-min(x)))
print("Euclid methode on abalone data")   
euclid_app(as.data.frame(aba_std),3)
```


# **3.  Outlier detection using PCA ** #

#### bivariate normal-distribution ####
par(mfcol=c(2,2))
print("Outlier detection using PCA with Factor 1.5 on bivariate normal-distribution")                 

PCA1.5_1 <-bp_app(as.data.frame(princomp(data_std)$scores)[,1:2],1.5)
plot(data,col=PCA1.5_1,pch=18)
print("Outlier detection using PCA with Factor 3 on bivariate normal-distribution")                 
PCA3_1 <-bp_app(as.data.frame(princomp(data_std)$scores)[,1:2],3)

### abalone ###
print("Outlier detection using PCA with Factor 3 on abalones data")                 
PCA3_2 <- bp_app(as.data.frame(princomp(aba_std)$scores)[,1:2],3)
plot(aba_out,col=PCA3_2,pch=18)
```




# **4. Outlier detection using K-Means** #


#### bivariate Normal-Distribution ####
par(mfcol= c(3,2))
cl1_1 <- kmeans(data_std,50)
ind <- as.vector(which(table(cl1_1$cluster)< 10))
out <- ifelse(cl1_1$cluster %in% ind,2,1)
plot(as.data.frame(data),pch=18,col=out, main="A1) bivariate normal-distr.:ncl=50 , minpoints=10")

cl2_1 <- kmeans(data_std,100)
ind <- as.vector(which(table(cl2_1$cluster)< 10))
out <- ifelse(cl2_1$cluster %in% ind,2,1)
plot(as.data.frame(data),pch=18,col=out, main="A2) bivariate normal-distr.:ncl=100 , minpoints=10")

cl3_1 <- kmeans(data_std,50)
ind <- as.vector(which(table(cl3_1$cluster)< 5))
out <- ifelse(cl3_1$cluster %in% ind,2,1)
plot(as.data.frame(data),pch=18,col=out, main="A3) bivariate normal-distr.:ncl=50 , minpoints=5")

### abalone ###
cl1_2 <- kmeans(aba_std,200)
ind <- as.vector(which(table(cl1_2$cluster)< 5))
out <- ifelse(cl1_2$cluster %in% ind,2,1)
plot(as.data.frame(aba_out),pch=18,col=out, main="B1) abalones:ncl=200 , minpoints=5")

cl3_2 <- kmeans(aba_std,300)
ind <- as.vector(which(table(cl3_2$cluster)< 5))
out <- ifelse(cl3_2$cluster %in% ind,2,1)
plot(as.data.frame(aba_out),pch=18,col=out, main="B2) abalones:ncl=300 , minpoints=5")

cl2_2 <- kmeans(aba_std,200)
ind <- as.vector(which(table(cl2_2$cluster)< 10))
out <- ifelse(cl2_2$cluster %in% ind,2,1)
plot(as.data.frame(aba_out),pch=18,col=out, main="B3) abalones:ncl=200 , minpoints=10")
```

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 So based on my results k-Means doesn't deliver a confidently outlier detection. The identified candidates need to be evaluated closely.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 # **5. Outlier detection using density-based Clustering** #
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 The last method I wanted to test is the density-based clustering. This method classifies observations as outliers if there are too few observations in their surrounding. The parameter $\text{eps}$ gives the radius of this neighborhood around the point. In this $\epsilon$-neighborhood other observations are searched . If their number is at least $\text{minPts}$ the point is a *core point*. If a points isn't a *core point* but lies in the $\epsilon$-neighborhood of a *core point* the point is  a *border point*. Every observation that neither is a *core point* nor a *border point* is an outlier.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Thus, to perform the density-based clustering we need to set this two parameter. 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Fortunately, there is an standardisized method to determine the optimal $\text{eps}$- value. 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 The idea is to calculate, the average of the distances of every point to its $k$ nearest neighbors. The value of $k$ needs to be set and corresponds to $\text{MinPts}$.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Next, these $k$-distances are plotted in an ascending order. The aim is to determine the elbow, which corresponds to the optimal $\text{eps}$ parameter.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 As I dont want to create the best clusters but to detect clear outliers I choose a value sleightly above of the elbow as $\text{eps}$.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 I tried $\text{MinPts}=10$ and was very much satified by the results. Again I used the standardized data:
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ```{r ,results='hold',fig.show='hold', fig.dim = c(10, 10)}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 #### bivariate Normal-Distribution ####
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 library("dbscan")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 par(mfcol=c(2,2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 kNNdistplot(data_std, k =  10)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 abline(h=0.1,lwd=2)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 dbcl1_1<-dbscan(data_std,eps=0.1,minPts=10)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 plot(as.data.frame(data),pch=18,col=ifelse(dbcl1_1$cluster==0,2,1))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ### abalone ###
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 kNNdistplot(aba_std, k =  10)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 abline(h=0.05,lwd=2)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 dbcl1_2 <- dbscan(aba_std,eps=0.05,minPts=10)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 plot(as.data.frame(aba_out),pch=18,col=ifelse(dbcl1_2$cluster==0,2,1))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ```
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 On the simulated  data the method caught the four constructed outliers as well as some other points. By inspecting the plot we can confirm that all these points are outliers. No point in the middle of the points-cloud are wrongly rated as outlier. Maybe a few more points could have been classified as outliers like the ones around $(x=15, y=0)$ but overall I am quite satisfied with the result. 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 On the abalones data only the five clear oultiers are detected.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 # **Summary** #
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 The clear winner is the density-based clustering. This method is the only one, which works well on both data sets.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 The Euclid method identifies only the four constructed points as outliers in the simulated data. Thus it works quite well to detect only huge outliers in uncorrrelated data. But on the correlated data it classifies points as outliers which are quite central. So this method is only useful in special cases.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 The results of the outlier detection using k-Means is based on chance and also the parameter needs to be set correctly. So this method isn't confidently.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Also on the correlated data even in the best results many observation were classified as outliers which are pretty close to the points-band.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 The boxplot method and the outlier detection using the PCA doesn't supplies satisfying results on the tested datasets.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 