> wssplot <- function(numsub, nc=15, seed=1234) {}
> wssplot <- function(numsub, nc=15, seed=1234) {
+ 	wss <- (nrow(numsub)-1)*sum(apply(numsub,2,var))
+	for (i in 2:nc) {
		set.seed(seed)
		wss[i] <- sum(kmeans(numsub, center=i)$withinss)}
	plot(1:nc, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")}

wssplot(numsub)