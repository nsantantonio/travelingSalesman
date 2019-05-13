library(txtplot)
parDir <- getwd()
usMap <- FALSE
usxy <- TRUE
if(usMap){

	library(maps)
	capitals <- us.cities[us.cities$capital > 0, ]
	capitals48 <- capitals[!capitals$country.etc %in% c("HI", "AK"), ]
	distance <- as.matrix(dist(as.matrix(capitals48[c("long", "lat")])))

	cityCoords <- capitals48
	cityCoords$x <- cityCoords$long
	cityCoords$y <- cityCoords$lat
	par(mar=c(0,0,0,0))
	# map('usa',col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05, mar=rep(0,4), border=0, ylim=c(-80,80))
	usa <- map("usa")
	usa <- map("usa", plot = FALSE)
	points(cityCoords$x, cityCoords$y, pch = 20)
	# thoughts
	#	- split

} else if(usxy) {

	distance <- as.matrix(read.table(paste0(parDir, "/data/att48_d.txt")))
	dimnames(distance) <- list(paste0("city", 1:nrow(distance)), paste0("city", 1:ncol(distance)))
	cityCoords <- as.matrix(read.table(paste0(parDir, "/data/att48_xy.txt")))
	dimnames(cityCoords) <- list(paste0("city", 1:nrow(distance)), c("x", "y"))
	optsol <- read.table(paste0(parDir, "/data/att48_s.txt"))[[1]]

} else {
	set.seed(987)
	nG <- 5
	npg <- 3
	x <- runif(nG)
	y <- runif(nG)
	x <- rep(x, each = npg) + runif(nG * npg, -0.1, 0.1)
	y <- rep(y, each = npg) + runif(nG * npg, -0.1, 0.1)

	plot(x, y)

	cityCoords <- cbind(x, y)
	distance <- as.matrix(dist(cityCoords))
	rownames(cityCoords) <- paste0("city", 1:nrow(cityCoords))
}



subSetProblem <- FALSE
if(subSetProblem) {
	n <- 20
	# set <- 10
	distance <- distance[1:n, 1:n]
	cityCoords <- cityCoords[1:n, ]
}

# range(as.matrix(dist(cityCoords)) - distance)




getIndex <- function(path, sym = FALSE, loop = FALSE) {
	last <- length(path)
	if(!returnToStart) last <- last - 1  
	nextIndex <- c(path[2:length(path)], path[1])
	index <- cbind(path[1:last], nextIndex[1:last])
	if(sym) index <- rbind(index, index[, c(2, 1)])
	index
}

getPath <- function(path, mat, f = identity, start = 1, shouldLoop = FALSE){
	if(path[1] != start) stop("path does not start at:", start)
	f(mat[getIndex(path, sym = FALSE, loop = shouldLoop)])
}


evap <- function(oldP, newP, evapRate) (1 - evapRate) * oldP + newP 
getProb <- function(tao, nu, alpha, beta) {
	# print(mean(tao^alpha) / mean(nu^beta))
	(tao^alpha * nu^beta) / sum(tao^alpha * nu^beta)
}


pherFunc <- function(x, p) (x / max(x))^p
# pherFunc <- function(x, p) x * p



# initFunc <- function(Q, distance) {
# 	distInv <- Q / distance
# 	# distInv <- 1 / distance
# 	diag(distInv) <- NA
# 	distInv
# }

# # sigmoid <- function(x) 1 / (1 + exp(x))


# # pherFunc <- function(x, p) (x / max(x))^p
# # pherFunc <- function(x, p) x * p

initFunc <- function(Q, distance) {
	distInv <- Q / distance
	diag(distInv) <- 0
	distInv
}

pherFuncSum <- function(x, Q, p) x
# pherFuncSum <- function(x, Q, p) (x / max(x))^p
# pherFuncAnt <- function(d, i, Q, p) {Q / d[[i]]}
pherFuncAnt <- function(d, i, Q, p) {
	mind <- d[[which.min(d)]]
	deltaTao <- d[[i]]
	Q / ((deltaTao/mind)^p * mind)
}

subsetAnts <- function(paths, quant = 1) {# return a boolean to indicate which ants are included
	paths <- unlist(paths)
	paths <= quantile(paths, quant)
}

if(usxy) minPath <- getPath(optsol, distance, f = sum)

qthresh <- 0.1 # works pretty well by subsetting to only the fastest ants. SHould try with all ants.
pherPower <- 1

lineBegin <- 400 
nCities <- ncol(distance)
evapRate <- 0.05
nAnts <- 500
startCity <- 1
returnToStart <- TRUE
nPaths <- nCities - 1

shortestPathLen <- sum(distance[lower.tri(distance)])
lastBestPath <- shortestPathLen
maxIter <- 1000
pathCounter <- 0
countThresh <- 10

# Q <- mean(distance) * nCities
Q <- 1

distInv <- initFunc(Q, distance)
pheromone <- matrix(1, nCities, nCities) - diag(nCities)

noPheromone <- matrix(0, nrow(distance), ncol(distance))

alpha <- 0.95 # pheromone
beta <- 1.05 # distance

iter <- 0
plotMap <- TRUE
plotpathLen <- TRUE
imagePheromone <- FALSE
estLine <- TRUE
bestPathLen <- NULL
meanPathLen <- NULL
bestIter <- NULL


# whichInit <- which(lower.tri(distInv), arr.ind = TRUE)[, c(2, 1)]
# barplot(distInv[whichInit])


plot.new()


# updateBeta <- FALSE

while(pathCounter <= countThresh & iter <= maxIter){

	# if(updateBeta & iter > 0) beta <- beta^(1/iter)

	iter <- iter + 1
	ants <- rep(list(startCity), nAnts)
	for (i in 1:nAnts) {
		for (j in 2:nCities) {
			pij <- pheromone[ants[[i]][[j-1]], ]
			pij[ants[[i]]] <- 0 
			distInvij <- distInv[ants[[i]][[j-1]], ]
			distInvij[ants[[i]]] <- 0 
			pathProb <- getProb(tao = pij, nu = distInvij, alpha = alpha, beta = beta) 
			ants[[i]][[j]] <- sample(1:length(pij), 1, prob = pathProb) 
		}
	}
	pathLen <- lapply(ants, getPath, distance, f = sum, shouldLoop = returnToStart)
	cleverAnts <- subsetAnts(pathLen, qthresh)
	pathLen <- pathLen[cleverAnts]
	ants <- ants[cleverAnts]
	# txtdensity(unlist(pathLen))


	shortPath <- which.min(pathLen)
	if (pathLen[[shortPath]] == lastBestPath) {
		pathCounter <- pathCounter + 1
	} else {
		lastBestPath <- pathLen[[shortPath]]
		pathCounter <- 0
	}

	if (pathLen[[shortPath]] < shortestPathLen){
		shortestPath <- shortPath
		shortestPathLen <- pathLen[[shortestPath]]
		bestAnt <- ants[[shortPath]]
		bestIter <- iter
		if (returnToStart) bestAnt <- c(bestAnt, startCity)
	}

	newPheromone <- list()
	for (i in 1:length(ants)) {
		newPheromone[[i]] <- noPheromone
		newPheromone[[i]][getIndex(ants[[i]], sym = TRUE, loop = returnToStart)] <- pherFuncAnt(pathLen, i, Q = Q, p = pherPower)
		# if(!isSymmetric(newPheromone[[i]])) stop("not sym!")
	}
	newPheromone[which.min(pathLen)]
	newPheromone[which.max(pathLen)]
	
	# need to update based on fastest ant, not shortest path!!!!!
	sumNewPheromone <- Reduce("+", newPheromone)
	# newPheromone <- pherFuncSum(sumNewPheromone, Q = Q, p = pherPowerS)
	newPheromoneSum <- pherFuncSum(sumNewPheromone)
	# txtdensity(newPheromone)

	pheromone <- evap(pheromone, newPheromoneSum, evapRate)

	# print(cor(pheromone[lower.tri(pheromone)], distInv[lower.tri(distInv)]))
	# print(pathLen[[shortPath]])

	# this is just to add start to end for plotting...
	whichBest <- apply(getIndex(ants[[shortPath]]), 1, sort)
	if (returnToStart) ants <- lapply(ants, function(x) c(x, startCity))

	# if (plotpathLen & plotMap) par(mfrow = c(2, 1))
	if (plotpathLen & plotMap & imagePheromone) par(mfrow = c(2, 3)) else if(plotpathLen & plotMap) par(mfrow = c(1, 3)) else if(plotpathLen) par(mfrow = c(1, 1)) 

	if (imagePheromone) {
		# par(mfrow = c(1, 3)
		image(distInv, main = "prior")
		image(pheromone,main = "pheromone")
		image(getProb(tao = pheromone, nu = distInv, alpha = alpha, beta = beta), main = "transition probability")
	}

	if (plotpathLen) {
		bestPathLen <- c(bestPathLen, pathLen[[shortPath]])
		meanPathLen <- c(meanPathLen, mean(unlist(pathLen)))
		
		plot(1:iter, meanPathLen, type = "b", pch = 16, ylim = range(c(bestPathLen, meanPathLen)))
		lines(1:iter, bestPathLen, type = "b", pch = 16, col = "firebrick")

		if(estLine & iter > lineBegin + 3) {
			x <- lineBegin:iter
			y <- meanPathLen[x]
			fit <- lm(y ~ x)
			abline(fit)
		}
	}
	if (plotMap) {
		whichrows <- which(lower.tri(pheromone), arr.ind = TRUE)[, c(2, 1)]
		bestIndex <- which(as.data.frame(t(whichrows)) %in% as.data.frame(whichBest))

		# cityCombos2 <- cbind(rownames(distance)[whichrows[, 1]], rownames(distance)[whichrows[, 2]])
		pl <- pheromone[whichrows]
		lw <- {pl - max(pl)} / {min(pl) - max(pl)}
		# lw <- {pl - min(pl)} / {max(pl) - min(pl)}
		barCols <- rep("gray20", choose(nCities, 2))
		barCols[bestIndex] <- "red"
		barplot(pheromone[whichrows], col = barCols)

		cityCombos <- combn(rownames(cityCoords), 2, simplify = FALSE)
		# do.call(rbind, cityCombos), cityCombos2

		pathx <- do.call(rbind, lapply(cityCombos, function(x) cityCoords[x, "x"]))
		pathy <- do.call(rbind, lapply(cityCombos, function(x) cityCoords[x, "y"]))

		if(usMap) map("usa") else plot(NA, xlim = range(cityCoords[, "x"]), ylim = range(cityCoords[, "y"]), xlab = "", ylab = "")
		# segments(pathx[, 1], pathy[, 1], pathx[, 2], pathy[, 2], col = rgb(lw, lw, lw), lwd = 2)
		lines(cityCoords[ants[[shortPath]], "x"], cityCoords[ants[[shortPath]], "y"], lwd = 2)
		lines(cityCoords[bestAnt, "x"], cityCoords[bestAnt, "y"], lwd = 1, col = "firebrick")
		points(cityCoords[, "x"], cityCoords[, "y"], pch = c(16, rep(16, nCities - 1)), col = c("red", rep("black", nCities - 1)), main = "Best ant")
	}	
	# Sys.sleep(0.5)
}

bestAnt
lastAnt <- ants[[shortPath]]

shortestPathLen - minPath

if(usxy){

	pdf(paste0("evap", evapRate, "_qthresh", qthresh, "_", maxIter, "iter_", nAnts, "Ants_pherPower", pherPower, "_alpha", alpha, "_beta", beta, ".pdf"))
	par(mfrow = c(2, 2))


	plot(1:iter, meanPathLen, type = "b", pch = 1, ylim = range(c(minPath, bestPathLen, meanPathLen)), cex = 0.5)
	lines(1:iter, bestPathLen, type = "b", pch = 1, col = "firebrick", cex = 0.5)
	abline(h = minPath, lty = 3)
	legend("topright", legend = c("Mean path length", "Shortest path length"), fill = c("black", "firebrick"))
	text(iter * 0.65, mean(range(meanPathLen)), label = paste0("Min diff from optimal: ", shortestPathLen - minPath))
	# pointCol <- rep("black", iter)
	# pointCol[c(shortestPath, shortPath)] <- "red"
	# plot(1:iter, bestPathLen, type = "b", col = pointCol, pch = 16)


	plot(cityCoords[, "x"], cityCoords[, "y"], pch = c(16, rep(1, nCities - 1)), main = "Optimal Solution", xlab = "", ylab = "")
	lines(cityCoords[optsol, "x"], cityCoords[optsol, "y"], col = "red")

	plot(cityCoords[, "x"], cityCoords[, "y"], pch = c(16, rep(1, nCities - 1)), , main = "My Best Ant", xlab = "", ylab = "")
	lines(cityCoords[bestAnt, "x"], cityCoords[bestAnt, "y"])

	plot(cityCoords[, "x"], cityCoords[, "y"], pch = c(16, rep(1, nCities - 1)), main = "My Last Best Ant", xlab = "", ylab = "")
	lines(cityCoords[lastAnt, "x"], cityCoords[lastAnt, "y"], col = "blue")
	
	dev.off()
}


















# matrixToLong <- function(X, valName = "val", inclDiag = FALSE, tri = upper.tri, sortX = TRUE) {
# 	Xlong <- data.frame(which(tri(X, diag = inclDiag), arr.ind = TRUE))
# 	Xlong[[valName]] <- X[tri(X, diag = inclDiag)]
# 	if(sortX) Xlong <- Xlong[order(Xlong[,1], Xlong[,2]), ]
# 	Xlong
# }

# d <- matrixToLong(dmat, "distance")
# edges <- as.matrix(d[c("row", "col")])
# Z <- diag()
# nodes <- unique(c(edges))
# distance  <- d[["distance"]] 
# # edgeN <- nrow(edges)
# # edgeN <- 1:nrow(edges)




# evap <- function

# getProb <- function(tao, nu, alpha, beta) tao^alpha * nu^beta / sum(tao^alpha * nu^beta)

# alpha <- 1
# beta <- 1
# # nu is attractiveness, tao is trail level
# antPath <- function(edges, nodes, start = 1, pheromone = NULL){

# 	edgeN <- 1:nrow(edges)
# 	nPaths <- length(nodes) - 1
# 	if(is.null(pheromone)) pheromone <- rep(1/nPaths, nrow(edges)) # could just rep 1
# 	path <- start
# 	i <- 1
# 	while(length(path) < nPaths + 1) {
# 		candidatePath <- rowSums(edges == path[i]) == 1
# 		# candidatePath <- edges[, 1] == curNode
# 		distinv <- 1 / distance[candidatePath]
# 		nu <- distinv / sum(distinv)
# 		tao <- pheromone[candidatePath] / sum(pheromone[candidatePath])
# 		p <- getProb(tao, nu, alpha, beta)
# 		samplePath <- sample(edgeN[candidatePath], 1, prob = p)
# 		path <- c(path, samplePath)
# 		i <- i+1
# 	}
# 	tao <-  



# }





# optsol <- c(1,3,2,5,4)

