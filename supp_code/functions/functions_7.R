# Libraries:

# install.packages(c("geoscale", "strap"))
library(strap)
# install.packages("devtools")
library(devtools)
# install_github("graemetlloyd/Claddis")
library(Claddis)

# Functions:

makedeep <- function(x, y, y_label, quartz=TRUE, main="")
{
  y.max=max(sort(y))
  y.min=0
  young.bin <- -8 # Offset for youngest bin label...allows for age < 0 to center label
  plot.min=y.min-.03*(y.max-y.min)
  seg.min=y.min-.07*(y.max-y.min)
  text.min=y.min-.036*(y.max-y.min)
  time.boundaries=c(5.332, 23.03, 33.9, 55.8, 65.5, 99.6, 145.5, 161.2, 175.6, 199.6)
  interval.names=c("Pl", "M", "O", "E", "P", "UK", "LK", "UJ", "MJ", "LJ")
  interval.midpoint=time.boundaries-diff(c(0, time.boundaries))/2
  plot(1, 1, xlim=c(max(time.boundaries), 0), ylim=c(plot.min, y.max), type="n", xlab="Geologic time (Ma)", ylab=y_label, main=main)
  abline(h=y.min)
  segments(c(time.boundaries, 0), y.min, c(time.boundaries, 0), seg.min)
  text(interval.midpoint, text.min, labels=interval.names)
  lines(x, y)
}

makedeeppoly <- function(x, y, y_label, upper95, lower95, quartz=TRUE, main="")
{
  y.max=max(c(sort(y), sort(upper95)))
  y.min=0
  young.bin <- -8 # Offset for youngest bin label...allows for age < 0 to center label
  plot.min=y.min-.03*(y.max-y.min)
  seg.min=y.min-.07*(y.max-y.min)
  text.min=y.min-.036*(y.max-y.min)
  time.boundaries=c(5.332, 23.03, 33.9, 55.8, 65.5, 99.6, 145.5, 161.2, 175.6, 199.6)
  interval.names=c("Pl", "M", "O", "E", "P", "UK", "LK", "UJ", "MJ", "LJ")
  interval.midpoint=time.boundaries-diff(c(0, time.boundaries))/2
  plot(1, 1, xlim=c(max(time.boundaries), 0), ylim=c(plot.min, y.max), type="n", xlab="Geologic time (Ma)", ylab=y_label, main=main)
  abline(h=y.min)
  segments(c(time.boundaries, 0), y.min, c(time.boundaries, 0), seg.min)
  text(interval.midpoint, text.min, labels=interval.names)
  lines(x, y)
  polygon(x=c(x, rev(x)), y=c(upper95, rev(lower95)), col="grey", border=NA)
  points(x, y, type="l")
}

makedeepz <- function(x, y, z, y_label, z_label, quartz=TRUE, main="")
{
	y.max=max(sort(y))
	y.min=0
	z.max=max(sort(z))
	z.min=0
	young.bin <- -8 # Offset for youngest bin label...allows for age < 0 to center label
	yplot.min=y.min-.03*(y.max-y.min)
	yseg.min=y.min-.07*(y.max-y.min)
	ytext.min=y.min-.036*(y.max-y.min)
	zplot.min=z.min-.03*(z.max-z.min)
	time.boundaries=c(5.332, 23.03, 33.9, 55.8, 65.5, 99.6, 145.5, 161.2, 175.6, 199.6)
	interval.names=c("Pl", "M", "O", "E", "P", "UK", "LK", "UJ", "MJ", "LJ")
	interval.midpoint=time.boundaries-diff(c(0, time.boundaries))/2
	par(mar=c(5, 4, 4, 4) + 0.1) # Leave space for z axis
	plot(1, 1, xlim=c(max(time.boundaries), 0), ylim=c(yplot.min, y.max), type="n", xlab="Geologic time (Ma)", ylab=y_label, main=main)
	abline(h=y.min)
	segments(c(time.boundaries, 0), y.min, c(time.boundaries, 0), yseg.min)
	text(interval.midpoint, ytext.min, labels=interval.names)
	lines(x, y)
	par(new=T)
	plot(1, 1, xlim=c(max(time.boundaries), 0), ylim=c(zplot.min, z.max), axes=F, type="n", xlab="", ylab="")
	lines(x, z, col="dark grey")
	axis(4, at=pretty(range(z)))
	mtext(z_label, 4, 3)
}

makephan <- function(x, y, y_label, quartz=TRUE, main="")
{
	y.max=max(sort(y))
	y.min=0
	young.bin <- -8 # Offset for youngest bin label...allows for age < 0 to center label
	plot.min=y.min-.03*(y.max-y.min)
	seg.min=y.min-.07*(y.max-y.min)
	text.min=y.min-.036*(y.max-y.min)
	time.boundaries=c(23.03, 65.5, 145.5, 199.6, 251, 299, 359.2, 416, 443.7, 488.3, 542)
	interval.names=c("Ng", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Cm")
	interval.midpoint=time.boundaries-diff(c(0, time.boundaries))/2
	plot(1, 1, xlim=c(max(time.boundaries), 0), ylim=c(plot.min, y.max), type="n", xlab="Geologic time (Ma)", ylab=y_label, main=main)
	abline(h=y.min)
	segments(c(time.boundaries, 0), y.min, c(time.boundaries, 0), seg.min)
	text(interval.midpoint, text.min, labels=interval.names)
	lines(x, y)
}

AICc <- function(model, n)
{
	require(stats)
	require(nlme)
	p <- length(coef(model))
	lk <- AIC(model)-(2*p)
	lk+(2*p*(n/(n-p-1)))
}

linear.model <- function(x, y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1 <- nls(y~a*x, start=list(a=(max(y)/max(x))), algorithm="port", lower=list(a=0), control=list(warnOnly=TRUE))
	modelout2 <- nls(y~(a*x)+b, start=list(a=(max(y)/max(x)), b=1), algorithm="port", lower=list(a=0, b=0), control=list(warnOnly=TRUE))
	modelout3 <- lm(y~x)
	wts <- akaike.wts(c(AICc(modelout1, n=length(x)), AICc(modelout2, n=length(x)), AICc(modelout3, n=length(x))))
	best <- grep(TRUE, max(wts) == wts)
	if(best == 1) model <- modelout1
	if(best == 2) model <- modelout2
	if(best == 3) model <- modelout3
	sefit.model <- std.error(y-predict(model))*1.96
	sdfit.model <- sd(y-predict(model))*1.96
	result <- list(model, sefit.model, sdfit.model)
	names(result) <- c("model", "sefit.model", "sdfit.model")
	return(result)
}

hyperbolic.model <- function(x, y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1 <- NA; modelout2 <- NA
	try(modelout1 <- nls(y~(a*x)/(b+x), start=list(a=max(y), b=1), algorithm="port", control=list(warnOnly=TRUE)), silent=TRUE)
	try(modelout2 <- nls(y~a+((c*x)/(b+x)), start=list(a=1, b=1, c=max(y)), algorithm="port", lower=list(a=0), control=list(warnOnly=TRUE)), silent=TRUE)
	if (is.na(modelout1[1]) == TRUE) modelout1 <- nls(y~(max(y)*x)/(b+x), start=list(b=1), algorithm="port", control=list(warnOnly=TRUE))
	if (is.na(modelout2[1]) == TRUE) modelout2 <- nls(y~a+((max(y)*x)/(b+x)), start=list(a=1, b=1), algorithm="port", lower=list(a=0), control=list(warnOnly=TRUE))
	wts <- akaike.wts(c(AICc(modelout1, n=length(x)), AICc(modelout2, n=length(x))))
	ifelse(wts[1] >= wts[2], model <- modelout1, model <- modelout2)
	sefit.model <- std.error(y-predict(model))*1.96
	sdfit.model <- sd(y-predict(model))*1.96
	result <- list(model, sefit.model, sdfit.model)
	names(result) <- c("model", "sefit.model", "sdfit.model")
	return(result)
}

logarithmic.model <- function(x, y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	model <- NA
	try(model <- nls(y~a+log(b+x), start=list(a=1, b=1), algorithm="port", lower=list(a=0, b=0), control=list(warnOnly=TRUE)), silent=TRUE)
	if (is.na(model[1]) == TRUE) model <- nls(y~a+log(1+x), start=list(a=1), algorithm="port", lower=list(a=0), control=list(warnOnly=TRUE))
	sefit.model <- std.error(y-predict(model))*1.96
	sdfit.model <- sd(y-predict(model))*1.96
	result <- list(model, sefit.model, sdfit.model)
	names(result) <- c("model", "sefit.model", "sdfit.model")
	return(result)
}

exponential.model <- function(x, y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1 <- NA; modelout2 <- NA
	try(modelout1 <- nls(y~a*(1-exp(-c*x)), start=list(a=max(y), c=0.1), algorithm="port", lower=list(a=0, c=0), control=list(warnOnly=TRUE)), silent=TRUE)
	try(modelout2 <- nls(y~a-(b*exp(-c*x)), start=list(a=max(y), b=max(y), c=0), algorithm="port", control=list(warnOnly=TRUE)), silent=TRUE)
	if (is.na(modelout1[1]) == TRUE) modelout1 <- nls(y~max(y)*(1-exp(-c*x)), start=list(c=0.000001), algorithm="port", lower=list(c=0), control=list(warnOnly=TRUE))
	if (is.na(modelout2[1]) == TRUE) modelout2 <- nls(y~max(y)-(max(y)*exp(-c*x)), start=list(c=0.), algorithm="port", control=list(warnOnly=TRUE))
	wts <- akaike.wts(c(AICc(modelout1, n=length(x)), AICc(modelout2, n=length(x))))
	ifelse(wts[1] >= wts[2], model <- modelout1, model <- modelout2)
	sefit.model <- std.error(y-predict(model))*1.96
	sdfit.model <- sd(y-predict(model))*1.96
	result <- list(model, sefit.model, sdfit.model)
	names(result) <- c("model", "sefit.model", "sdfit.model")
	return(result)
}

sigmoidal.model <- function(x, y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1 <- NA; modelout2 <- NA
	try(modelout1 <- nls(y~b/(c+exp(-x)), start=list(b=max(y)/100, c=max(y)/10000), algorithm="port", control=list(warnOnly=TRUE)), silent=TRUE)
	try(modelout2 <- nls(y~a+(b/(c+exp(-x))), start=list(a=1, b=0.1, c=0.001), algorithm="port", lower=list(a=0), control=list(warnOnly=TRUE)), silent=TRUE)
	if (is.na(modelout1[1]) == TRUE) modelout1 <- nls(y~b/((1/max(y))+exp(-x)), start=list(b=1), algorithm="port", control=list(warnOnly=TRUE))
	if (is.na(modelout2[1]) == TRUE) modelout2 <- nls(y~a+(max(y)/100)/((max(y)/100)+exp(-x)), start=list(a=1), algorithm="port", lower=list(a=0), control=list(warnOnly=TRUE))
	wts <- akaike.wts(c(AICc(modelout1, n=length(x)), AICc(modelout2, n=length(x))))
	ifelse(wts[1] >= wts[2], model <- modelout1, model <- modelout2)
	sefit.model <- std.error(y-predict(model))*1.96
	sdfit.model <- sd(y-predict(model))*1.96
	result <- list(model, sefit.model, sdfit.model)
	names(result) <- c("model", "sefit.model", "sdfit.model")
	return(result)
}

polynomial.model <- function(x, y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout2 <- lm(y~x+I(x^2))
	modelout3 <- lm(y~x+I(x^2)+I(x^3))
	modelout4 <- lm(y~x+I(x^2)+I(x^3)+I(x^4))
	wts <- akaike.wts(c(AICc(modelout2, n=length(x)), AICc(modelout3, n=length(x)), AICc(modelout4, n=length(x))))
	best <- grep(TRUE, max(wts) == wts)
	if(best == 1) model <- modelout2
	if(best == 2) model <- modelout3
	if(best == 3) model <- modelout4
	sefit.model <- std.error(y-predict(model))*1.96
	sdfit.model <- sd(y-predict(model))*1.96
	result <- list(model, sefit.model, sdfit.model)
	names(result) <- c("model", "sefit.model", "sdfit.model")
	return(result)
}

best.model <- function(x, y)
{
	linmod <- linear.model(x, y)$model # Fit linear model
	hypmod <- hyperbolic.model(x, y)$model # Fit hyperbolic model
	logmod <- logarithmic.model(x, y)$model # Fit logarithmic model
	expmod <- exponential.model(x, y)$model # Fit exponential model
	sigmod <- sigmoidal.model(x, y)$model # Fit sigmoidal model
	polmod <- polynomial.model(x, y)$model # Fit polynomial model
	best <- min(c(AICc(linmod, n=length(x)), AICc(hypmod, n=length(x)), AICc(logmod, n=length(x)), AICc(expmod, n=length(x)), AICc(sigmod, n=length(x)), AICc(polmod, n=length(x))))
	if (AICc(linmod, n=length(x)) == best) model <- linmod; sefit.model <- linear.model(x, y)$sefit.model; sdfit.model <- linear.model(x, y)$sdfit.model
	if (AICc(hypmod, n=length(x)) == best) model <- hypmod; sefit.model <- hyperbolic.model(x, y)$sefit.model; sdfit.model <- hyperbolic.model(x, y)$sdfit.model
	if (AICc(logmod, n=length(x)) == best) model <- logmod; sefit.model <- logarithmic.model(x, y)$sefit.model; sdfit.model <- logarithmic.model(x, y)$sdfit.model
	if (AICc(expmod, n=length(x)) == best) model <- expmod; sefit.model <- exponential.model(x, y)$sefit.model; sdfit.model <- exponential.model(x, y)$sdfit.model
	if (AICc(sigmod, n=length(x)) == best) model <- sigmod; sefit.model <- sigmoidal.model(x, y)$sefit.model; sdfit.model <- sigmoidal.model(x, y)$sdfit.model
	if (AICc(polmod, n=length(x)) == best) model <- polmod; sefit.model <- polynomial.model(x, y)$sefit.model; sdfit.model <- polynomial.model(x, y)$sdfit.model
	result <- list(model, sefit.model, sdfit.model)
	names(result) <- c("model", "sefit.model", "sdfit.model")
	return(result)
}

rockmodel.predictCI <- function(rockmeasure, diversitymeasure, CI=0.95)
{
	x <- sort(rockmeasure)
	y <- sort(diversitymeasure)
	model <- best.model(x, y)$model
	sefit.model <- best.model(x, y)$sefit.model
	sdfit.model <- best.model(x, y)$sdfit.model
	predicted <- predict(model, list(x=rockmeasure))
	selowerCI <- predicted-sefit.model
	seupperCI <- predicted+sefit.model
	sdlowerCI <- predicted-sdfit.model
	sdupperCI <- predicted+sdfit.model
	result <- list(predicted, selowerCI, seupperCI, sdlowerCI, sdupperCI, model)
	names(result) <- c("predicted", "selowerCI", "seupperCI", "sdlowerCI", "sdupperCI", "model")
	return(result)
}

sphpolyarea <- function(vlat, vlon)
{
	require(fossil)
	splits <- strsplit(unique(paste(vlon, vlat, sep="_")), "_")
	vlat <- vlon <- vector(mode="numeric")
	for (i in 1:length(splits)) {
		vlon[i] <- splits[[i]][1]
		vlat[i] <- splits[[i]][2]
	}
	vlat <- as.numeric(vlat)
	vlon <- as.numeric(vlon)
	points <- chull(vlon, vlat)
	vlon <- vlon[points]
	vlat <- vlat[points]
	sum <- 0
	nv <- length(vlat)
	R <- 40041.47/(2 * pi)
	while (nv >= 3) {
		lat1 <- vlat[1]
		lat2 <- vlat[2]
		lat3 <- vlat[3]
		long1 <- vlon[1]
		long2 <- vlon[2]
		long3 <- vlon[3]
		cx <- deg.dist(lat1, long1, lat2, long2)/R
		bx <- deg.dist(lat1, long1, lat3, long3)/R
		ax <- deg.dist(lat2, long2, lat3, long3)/R
		A <- acos((cos(ax) - cos(bx) * cos(cx))/(sin(bx) * sin(cx)))
		B <- acos((cos(bx) - cos(cx) * cos(ax))/(sin(cx) * sin(ax)))
		C <- acos((cos(cx) - cos(ax) * cos(bx))/(sin(ax) * sin(bx)))
		SA <- R^2 * ((A + B + C) - pi)
		sum <- sum+SA
		vlat <- vlat[-2]
		vlon <- vlon[-2]
		nv <- length(vlat)
	}
	return(sum)
}

subsample <- function(pool, ntrials, CI = 0.95)
{
  spgout <- gout <- sout <- matrix(nrow=ntrials, ncol=length(pool))
  stopat <- length(unique(pool)) # No point continuing to pull out species if they have all been sampled
  poolholder <- pool
  topCI <- ceiling((1-((1-CI)/2))*ntrials)
  bottomCI <- max(floor(((1-CI)/2)*ntrials), 1)
  for (i in 1:ntrials) {
    pool <- poolholder
    poolsize <- length(pool)
    picked <- vector(mode="character") # Binomial
    pickedgen <- vector(mode="character") # Genus only
    counter <- 1
    while (poolsize > 0 && length(picked) < stopat) {
      pickno <- ceiling(runif(1, 0, length(pool)))
      picked <- unique(c(picked, pool[pickno]))
      pickedgen <- unique(c(pickedgen, strsplit(pool[pickno], " ")[[1]][1]))
      sout[i, counter] <- length(picked)
      gout[i, counter] <- length(pickedgen)
      spgout[i, counter] <- length(picked)/length(pickedgen)
      pool <- pool[-pickno]
      poolsize <- length(pool)
      counter <- counter+1
    }
    sout[i, grep(TRUE, is.na(sout[i, ]))] <- length(picked) # Fill in remaining
    gout[i, grep(TRUE, is.na(gout[i, ]))] <- length(pickedgen) # Fill in remaining
    spgout[i, grep(TRUE, is.na(spgout[i, ]))] <- length(picked)/length(pickedgen) # Fill in remaining
  }
  smeanCI <- supperCI <- slowerCI <- apply(sout, 2, mean)
  gmeanCI <- gupperCI <- glowerCI <- apply(gout, 2, mean)
  spgmeanCI <- spgupperCI <- spglowerCI <- apply(spgout, 2, mean)
  for (i in 1:length(sout[1, ])) {
    supperCI[i] <- sort(sout[, i])[topCI]
    slowerCI[i] <- sort(sout[, i])[bottomCI]
    gupperCI[i] <- sort(gout[, i])[topCI]
    glowerCI[i] <- sort(gout[, i])[bottomCI]
    spgupperCI[i] <- sort(spgout[, i])[topCI]
    spglowerCI[i] <- sort(spgout[, i])[bottomCI]
  }
  sout <- rbind(supperCI, smeanCI, slowerCI)
  gout <- rbind(gupperCI, gmeanCI, glowerCI)
  spgout <- rbind(spgupperCI, spgmeanCI, spglowerCI)
  out <- list(sout, gout, spgout)
  names(out) <- c("sout", "gout", "spgout")
  return(out)
}

five.models <- function(intable, time, log.input=FALSE)
{
  tstable <- matrix(ncol=4, nrow=length(time))
  if (log.input == FALSE) {
    for (i in 1:length(time)) {
      tstable[i, 1] <- ifelse(is.nan(mean(sort(intable[, i]))), 0, mean(sort(intable[, i])))
      tstable[i, 2] <- ifelse(is.nan(var(sort(intable[, i]))), 0, var(sort(intable[, i])))
      tstable[i, 3] <- length(sort(intable[, i]))
      tstable[i, 4] <- time[i]
    }
  }
  if (log.input == TRUE) {
    for (i in 1:length(time)) {
      tstable[i, 1] <- mean(sort(log(intable[, i]))[grep(TRUE, sort(log(intable[, i])) > 0)])
      tstable[i, 2] <- var(sort(log(intable[, i]))[grep(TRUE, sort(log(intable[, i])) > 0)])
      tstable[i, 3] <- length(sort(log(intable[, i]))[grep(TRUE, sort(log(intable[, i])) > 0)])
      tstable[i, 4] <- time[i]
    }
  }
  for (i in length(tstable[, 1]):1) {
    kill <- 0; a <- tstable[i, 1]; b <- tstable[i, 2]; c <- tstable[i, 3]
    if(is.nan(a) || is.na(a) || a == Inf || a == -Inf || a == 0) kill <- 1
    if(is.nan(b) || is.na(b) || b == Inf || b == -Inf || b == 0) kill <- 1
    if(is.nan(c) || is.na(c) || c == Inf || c == -Inf || c == 0) kill <- 1
    if(kill == 1) tstable <- tstable[-i, ]
  }
  tstable <- as.paleoTS(tstable[, 1], tstable[, 2], tstable[, 3], tstable[, 4], start.age=max(tstable[, 4]))
  threemodels <- fit3models(tstable, silent=T)
  twostases <- fitGpunc(tstable, ng=2, oshare=F, silent=T)
  threestases <- fitGpunc(tstable, ng=3, oshare=F, silent=T)
  results <- cbind(c(threemodels$aic, twostases$AIC, threestases$AIC), c(threemodels$aicc, twostases$AICc, threestases$AICc), round(akaike.wts(c(threemodels$aicc, twostases$AICc, threestases$AICc)), digits=5))
  colnames(results) <- c("AIC", "AICc", "Akaike Weights"); rownames(results) <- c("Directional trend", "Random walk", "Stasis", "Punctuation with two stases", "Punctuation with three stases")
  result <- list(tstable, threemodels, twostases, threestases, results, log.input)
  names(result) <- c("tstable", "threemodels", "twostases", "threestases", "results", "log.input")
  return(result)
}

plot.stases <- function(fivemodel, rawtable, y_label)
{
  makedeep(x=fivemodel$tstable$tt, y=fivemodel$tstable$mm, y_label=y_label)
  fivemodel$results
  bestmodel <- rownames(fivemodel$results)[grep(TRUE, fivemodel$results[, "Akaike Weights"] == max(fivemodel$results[, "Akaike Weights"]))]
  text(max(fivemodel$tstable$tt), max(fivemodel$tstable$mm), paste("Best model: ", bestmodel, sep=""), pos=4)
  if (bestmodel == "Punctuation with three stases") {
    firststasis <- 1:(fivemodel$threestases$shift.start[1]-1)
    secondstasis <- fivemodel$threestases$shift.start[1]:(fivemodel$threestases$shift.start[2]-1)
    thirdstasis <- fivemodel$threestases$shift.start[2]:length(fivemodel$tstable$mm)
  }
  if (bestmodel == "Punctuation with two stases") {
    firststasis <- 1:(fivemodel$twostases$shift.start[1]-1)
    secondstasis <- fivemodel$twostases$shift.start[1]:length(fivemodel$tstable$mm)
  }
  SD <- rawtable
  if (fivemodel$log.input == FALSE) {
    for (i in 1:length(SD[1, ])) SD[1, i] <- sd(sort(SD[, i]))
  }
  if (fivemodel$log.input == TRUE) {
    for (i in 1:length(SD[1, ])) SD[1, i] <- sd(log(sort(SD[, i])))
  }
  SD <- 1.96*SD[1, 1:length(fivemodel$tstable$mm)]
  if (bestmodel == "Punctuation with three stases") {
    polygon(c(fivemodel$tstable$tt[firststasis], rev(fivemodel$tstable$tt[firststasis])), c(rep(fivemodel$threestases$par[1]+fivemodel$threestases$par[4], length(firststasis)), rep(fivemodel$threestases$par[1]-fivemodel$threestases$par[4], length(firststasis))), col="grey", border=NA)
    polygon(c(fivemodel$tstable$tt[secondstasis], rev(fivemodel$tstable$tt[secondstasis])), c(rep(fivemodel$threestases$par[2]+fivemodel$threestases$par[5], length(secondstasis)), rep(fivemodel$threestases$par[2]-fivemodel$threestases$par[5], length(secondstasis))), col="grey", border=NA)
    polygon(c(fivemodel$tstable$tt[thirdstasis], rev(fivemodel$tstable$tt[thirdstasis])), c(rep(fivemodel$threestases$par[3]+fivemodel$threestases$par[6], length(thirdstasis)), rep(fivemodel$threestases$par[3]-fivemodel$threestases$par[6], length(thirdstasis))), col="grey", border=NA)
  }
  if (bestmodel == "Punctuation with two stases") {
    polygon(c(fivemodel$tstable$tt[firststasis], rev(fivemodel$tstable$tt[firststasis])), c(rep(fivemodel$twostases$par[1]+fivemodel$twostases$par[4], length(firststasis)), rep(fivemodel$twostases$par[1]-fivemodel$twostases$par[4], length(firststasis))), col="grey", border=NA)
    polygon(c(fivemodel$tstable$tt[secondstasis], rev(fivemodel$tstable$tt[secondstasis])), c(rep(fivemodel$twostases$par[2]+fivemodel$twostases$par[5], length(secondstasis)), rep(fivemodel$twostases$par[2]-fivemodel$twostases$par[5], length(secondstasis))), col="grey", border=NA)
  }
  points(x=fivemodel$tstable$tt, y=fivemodel$tstable$mm, type="l")
  points(fivemodel$tstable$tt, fivemodel$tstable$mm, cex=0.5)
  for (i in 1:length(fivemodel$tstable$mm)) lines(c(fivemodel$tstable$tt[i], fivemodel$tstable$tt[i]), c(fivemodel$tstable$mm[i]+SD[i], fivemodel$tstable$mm[i]-SD[i]))
  if (bestmodel == "Punctuation with three stases") {
    points(fivemodel$tstable$tt[firststasis], rep(fivemodel$threestases$par[1], length(firststasis)), type="l")
    points(fivemodel$tstable$tt[secondstasis], rep(fivemodel$threestases$par[2], length(secondstasis)), type="l")
    points(fivemodel$tstable$tt[thirdstasis], rep(fivemodel$threestases$par[3], length(thirdstasis)), type="l")
  }
  if (bestmodel == "Punctuation with two stases") {
    points(fivemodel$tstable$tt[firststasis], rep(fivemodel$twostases$par[1], length(firststasis)), type="l")
    points(fivemodel$tstable$tt[secondstasis], rep(fivemodel$twostases$par[2], length(secondstasis)), type="l")
  }
}

ts.corr <- function(var1, var2)
{
  shared.time <- intersect(var1$tstable$tt, var2$tstable$tt)
  var1.comp <- var1$tstable$mm[match(shared.time, var1$tstable$tt)]
  var2.comp <- var2$tstable$mm[match(shared.time, var2$tstable$tt)]
  raw.corr <- cor.test(var1.comp, var2.comp, method="spearman")
  firstdiff.corr <- cor.test(diff(var1.comp), diff(var2.comp), method="spearman")
  result <- list(raw.corr, firstdiff.corr)
  names(result) <- c("raw.corr", "firstdiff.corr")
  return(result)
}

multi.ts.3comp <- function(var1, var2, var3, var2name, var3name, cuts=NA) # First variable is dependent, second and third are explanatory
{
	require(paleoTS)
	# Create comparable vectors (i.e. where all variables are sampled in bin)
	shared.time <- intersect(intersect(var1$tstable$tt, var2$tstable$tt), var3$tstable$tt)
	if(!is.na(cuts)[1]) {
	  bottom <- max(cuts)
	  top <- min(cuts)
	  shared.time <- shared.time[intersect(grep(TRUE, shared.time >= top), grep(TRUE, shared.time <= bottom))]
	}
	var1.comp <- var1$tstable$mm[match(shared.time, var1$tstable$tt)]
	var2.comp <- var2$tstable$mm[match(shared.time, var2$tstable$tt)]
	var3.comp <- var3$tstable$mm[match(shared.time, var3$tstable$tt)]
	fd.var1.comp <- diff(var1.comp)
	fd.var2.comp <- diff(var2.comp)
	fd.var3.comp <- diff(var3.comp)
	gd.var1.comp <- gen.diff(var1.comp, shared.time)
	gd.var2.comp <- gen.diff(var2.comp, shared.time)
	gd.var3.comp <- gen.diff(var3.comp, shared.time)
	# Fit linear models:
	lm.var2 <- lm(var1.comp ~ var2.comp)
	lm.var3 <- lm(var1.comp ~ var3.comp)
	lm.var2var3 <- lm(var1.comp ~ var2.comp + var3.comp)
	fd.lm.var2 <- lm(fd.var1.comp ~ fd.var2.comp)
	fd.lm.var3 <- lm(fd.var1.comp ~ fd.var3.comp)
	fd.lm.var2var3 <- lm(fd.var1.comp ~ fd.var2.comp + fd.var3.comp)
	gd.lm.var2 <- lm(gd.var1.comp ~ gd.var2.comp)
	gd.lm.var3 <- lm(gd.var1.comp ~ gd.var3.comp)
	gd.lm.var2var3 <- lm(gd.var1.comp ~ gd.var2.comp + gd.var3.comp)
	# Get proportion explained by each model, plus unexplained:
	var2.rsq <- summary(lm.var2var3)$r.squared-summary(lm.var3)$r.squared
	var3.rsq <- summary(lm.var2var3)$r.squared-summary(lm.var2)$r.squared
	var2var3.rsq <- summary(lm.var2)$r.squared-var2.rsq
	ue.rsq <- 1-(var2.rsq+var3.rsq+var2var3.rsq)
	fd.var2.rsq <- summary(fd.lm.var2var3)$r.squared-summary(fd.lm.var3)$r.squared
	fd.var3.rsq <- summary(fd.lm.var2var3)$r.squared-summary(fd.lm.var2)$r.squared
	fd.var2var3.rsq <- summary(fd.lm.var2)$r.squared-fd.var2.rsq
	fd.ue.rsq <- 1-(fd.var2.rsq+fd.var3.rsq+fd.var2var3.rsq)
	gd.var2.rsq <- summary(gd.lm.var2var3)$r.squared-summary(gd.lm.var3)$r.squared
	gd.var3.rsq <- summary(gd.lm.var2var3)$r.squared-summary(gd.lm.var2)$r.squared
	gd.var2var3.rsq <- summary(gd.lm.var2)$r.squared-gd.var2.rsq
	gd.ue.rsq <- 1-(gd.var2.rsq+gd.var3.rsq+gd.var2var3.rsq)
	# Get AIC and AICc for each model plus Akaike weights:
	AIC.var2 <- AIC(lm.var2)
	AIC.var3 <- AIC(lm.var3)
	AIC.var2var3 <- AIC(lm.var2var3)
	fd.AIC.var2 <- AIC(fd.lm.var2)
	fd.AIC.var3 <- AIC(fd.lm.var3)
	fd.AIC.var2var3 <- AIC(fd.lm.var2var3)
	gd.AIC.var2 <- AIC(gd.lm.var2)
	gd.AIC.var3 <- AIC(gd.lm.var3)
	gd.AIC.var2var3 <- AIC(gd.lm.var2var3)
	AICc.var2 <- AICc(lm.var2, length(var1.comp))
	AICc.var3 <- AICc(lm.var3, length(var1.comp))
	AICc.var2var3 <- AICc(lm.var2var3, length(var1.comp))
	fd.AICc.var2 <- AICc(fd.lm.var2, length(fd.var1.comp))
	fd.AICc.var3 <- AICc(fd.lm.var3, length(fd.var1.comp))
	fd.AICc.var2var3 <- AICc(fd.lm.var2var3, length(fd.var1.comp))
	gd.AICc.var2 <- AICc(gd.lm.var2, length(gd.var1.comp))
	gd.AICc.var3 <- AICc(gd.lm.var3, length(gd.var1.comp))
	gd.AICc.var2var3 <- AICc(gd.lm.var2var3, length(gd.var1.comp))
	wts <- akaike.wts(c(AICc.var2, AICc.var3, AICc.var2var3))
	fd.wts <- akaike.wts(c(fd.AICc.var2, fd.AICc.var3, fd.AICc.var2var3))
	gd.wts <- akaike.wts(c(gd.AICc.var2, gd.AICc.var3, gd.AICc.var2var3))
	# Write formulae for models:
	fm.var2 <- paste("(", round(coef(lm.var2)[2], 2), " x ", var2name, ") + ", round(coef(lm.var2)[1], 2), sep="")
	fm.var3 <- paste("(", round(coef(lm.var3)[2], 2), " x ", var3name, ") + ", round(coef(lm.var3)[1], 2), sep="")
	fm.var2var3 <- paste("(", round(coef(lm.var2var3)[2], 2), " x ", var2name, ") + (", round(coef(lm.var2var3)[3], 2), " x ", var3name, ") + ", round(coef(lm.var2var3)[1], 2), sep="")
	fd.fm.var2 <- paste("(", round(coef(fd.lm.var2)[2], 2), " x ", var2name, ") + ", round(coef(fd.lm.var2)[1], 2), sep="")
	fd.fm.var3 <- paste("(", round(coef(fd.lm.var3)[2], 2), " x ", var3name, ") + ", round(coef(fd.lm.var3)[1], 2), sep="")
	fd.fm.var2var3 <- paste("(", round(coef(fd.lm.var2var3)[2], 2), " x ", var2name, ") + (", round(coef(fd.lm.var2var3)[3], 2), " x ", var3name, ") + ", round(coef(fd.lm.var2var3)[1], 2), sep="")
	gd.fm.var2 <- paste("(", round(coef(gd.lm.var2)[2], 2), " x ", var2name, ") + ", round(coef(gd.lm.var2)[1], 2), sep="")
	gd.fm.var3 <- paste("(", round(coef(gd.lm.var3)[2], 2), " x ", var3name, ") + ", round(coef(gd.lm.var3)[1], 2), sep="")
	gd.fm.var2var3 <- paste("(", round(coef(gd.lm.var2var3)[2], 2), " x ", var2name, ") + (", round(coef(gd.lm.var2var3)[3], 2), " x ", var3name, ") + ", round(coef(gd.lm.var2var3)[1], 2), sep="")
	# Make table:
	models <- c(var2name, var3name, paste(var2name, " + ", var3name, sep=""), "Unexplained")
	results <- cbind(models, c(fm.var2, fm.var3, fm.var2var3, NA), c(fd.fm.var2, fd.fm.var3, fd.fm.var2var3, NA), c(gd.fm.var2, gd.fm.var3, gd.fm.var2var3, NA), c(var2.rsq, var3.rsq, var2var3.rsq, ue.rsq), c(fd.var2.rsq, fd.var3.rsq, fd.var2var3.rsq, fd.ue.rsq), c(gd.var2.rsq, gd.var3.rsq, gd.var2var3.rsq, gd.ue.rsq), c(AIC.var2, AIC.var3, AIC.var2var3, NA), c(fd.AIC.var2, fd.AIC.var3, fd.AIC.var2var3, NA), c(gd.AIC.var2, gd.AIC.var3, gd.AIC.var2var3, NA), c(AICc.var2, AICc.var3, AICc.var2var3, NA), c(fd.AICc.var2, fd.AICc.var3, fd.AICc.var2var3, NA), c(gd.AICc.var2, gd.AICc.var3, gd.AICc.var2var3, NA), c(wts, NA), c(fd.wts, NA), c(gd.wts, NA))
	colnames(results) <- c("Explanatory variable", "Model formula (raw data)", "Model formula (first differences)", "Model formula (generalised differences)", "Proportion variance explained (raw data)", "Proportion variance explained (first differences)", "Proportion variance explained (generalised differences)", "AIC (raw data)", "AIC (first differences)", "AIC (generalised differences)", "AICc (raw data)", "AICc (first differences)", "AICc (generalised differences)", "Akaike weight (raw data)", "Akaike weight (first differences)", "Akaike weight (generalised differences)")
	result <- list(shared.time, var1.comp, var2.comp, var3.comp, fd.var1.comp, fd.var2.comp, fd.var3.comp, gd.var1.comp, gd.var2.comp, gd.var3.comp, lm.var2, lm.var3, lm.var2var3, fd.lm.var2, fd.lm.var3, fd.lm.var2var3, gd.lm.var2, gd.lm.var3, gd.lm.var2var3, results)
	names(result) <- c("shared.time", "var1.comp", "var2.comp", "var3.comp", "fd.var1.comp", "fd.var2.comp", "fd.var3.comp", "gd.var1.comp", "gd.var2.comp", "gd.var3.comp", "lm.var2", "lm.var3", "lm.var2var3", "fd.lm.var2", "fd.lm.var3", "fd.lm.var2var3", "gd.lm.var2", "gd.lm.var3", "gd.lm.var2var3", "results")
	return(result)
}

moving.average <- function(ts.holder, mav)
{
  ts <- vector(mode="numeric")
  for (i in 1:length(ts.holder[1, ])) ts[i] <- mean(sort(ts.holder[, i]))
  trimtb <- floor(mav/2)
  mav.vector <- vector(mode="numeric", length=length(ts)-(2*trimtb))
  for (i in (trimtb+1):(length(ts)-trimtb)) mav.vector[(i-trimtb)] <- mean(sort(ts[(i-trimtb):(i+trimtb)]))
  mav.vector <- as.numeric(gsub(NaN, NA, mav.vector))
  mav.vector <- c(rep(NA, trimtb), mav.vector, rep(NA, trimtb))
  mav.vector
}

moving.averagets <- function(ts, mav)
{
	trimtb <- floor(mav/2)
	mav.vector <- vector(mode="numeric", length=length(ts)-(2*trimtb))
	for (i in (trimtb+1):(length(ts)-trimtb)) mav.vector[(i-trimtb)] <- mean(sort(ts[(i-trimtb):(i+trimtb)]))
	mav.vector <- as.numeric(gsub(NaN, NA, mav.vector))
	mav.vector <- c(rep(NA, trimtb), mav.vector, rep(NA, trimtb))
	mav.vector
}

corr.plot <- function(x.holder, y.holder, time, log.x=FALSE, log.y=FALSE, fd=FALSE, gd=FALSE, do.mav=FALSE, mav, xlab="", ylab="", main="")
{
	y <- x <- vector(mode="numeric")
	for (i in 1:length(x.holder[1, ])) x[i] <- mean(sort(x.holder[, i]))
	for (i in 1:length(y.holder[1, ])) y[i] <- mean(sort(y.holder[, i]))
	if(log.x == FALSE) x1 <- as.numeric(gsub(NaN, NA, x))
	if(log.x == TRUE) x1 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, x)))))
	if(log.y == FALSE) y1 <- as.numeric(gsub(NaN, NA, y))
	if(log.y == TRUE) y1 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, y)))))
	if(fd == TRUE) {
		x1 <- diff(x1)
		y1 <- diff(y1)
	}
	if(gd == TRUE) {
		x1 <- gen.diff(x1, time)
		y1 <- gen.diff(y1, time)
	}	
	if(do.mav == TRUE) {
		trimtb <- floor(mav/2)
		mav.x1 <- vector(mode="numeric", length=length(x1)-(2*trimtb))
		mav.y1 <- vector(mode="numeric", length=length(y1)-(2*trimtb))
		for (i in (trimtb+1):(length(x1)-trimtb)) mav.x1[(i-trimtb)] <- mean(sort(x1[(i-trimtb):(i+trimtb)]))
		for (i in (trimtb+1):(length(y1)-trimtb)) mav.y1[(i-trimtb)] <- mean(sort(y1[(i-trimtb):(i+trimtb)]))
		mav.x1 <- as.numeric(gsub(NaN, NA, mav.x1))
		mav.y1 <- as.numeric(gsub(NaN, NA, mav.y1))
		mav.x1 <- c(rep(NA, trimtb), mav.x1, rep(NA, trimtb))
		mav.y1 <- c(rep(NA, trimtb), mav.y1, rep(NA, trimtb))
		x1 <- x1-mav.x1
		y1 <- y1-mav.y1
	}
	plot(x1, y1, xlab=xlab, ylab=ylab, main=main)
	text(x=min(sort(x1)), y=max(sort(y1)), labels=paste("Spearman rho: ", round(cor.test(x1, y1, method="spearman")$estimate, 2), " (p = ", format(round(cor.test(x1, y1, method="spearman")$p.value, 4), nsmall = 4), ")", sep="", collapse=""), adj=c(0, 1))
	abline(lsfit(x1, y1)$coefficients[1], lsfit(x1, y1)$coefficients[2])
}

ts.corr.plot <- function(var1.holder, var2.holder, var3.holder, log.var1=FALSE, log.var2=FALSE, log.var3=FALSE, var1.name, var2.name, var3.name, time)
{
  var1 <- var2 <- var3 <- vector(mode="numeric")
  for (i in 1:length(var1.holder[1, ])) var1[i] <- mean(sort(var1.holder[, i]))
  for (i in 1:length(var2.holder[1, ])) var2[i] <- mean(sort(var2.holder[, i]))
  for (i in 1:length(var3.holder[1, ])) var3[i] <- mean(sort(var3.holder[, i]))
  if(log.var1 == FALSE) var1 <- as.numeric(gsub(NaN, NA, var1))
  if(log.var1 == TRUE) var1 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, var1)))))
  if(log.var2 == FALSE) var2 <- as.numeric(gsub(NaN, NA, var2))
  if(log.var2 == TRUE) var2 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, var2)))))
  if(log.var3 == FALSE) var3 <- as.numeric(gsub(NaN, NA, var3))
  if(log.var3 == TRUE) var3 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, var3)))))
  par(mfrow=c(2, 2))
  plot(var1, var2, cex=0.5, xlab=var1.name, ylab=var2.name)
  abline(a=lsfit(var1, var2)$coefficients[1], b=lsfit(var1, var2)$coefficients[2])
  text(x=min(sort(var1)), y=max(sort(var2)), paste("rho = ", round(cor.test(var1, var2, method="spearman")$estimate, 2)), adj=c(0, 1))
  plot(gen.diff(var1, time), gen.diff(var2, time), cex=0.5, xlab=var1.name, ylab=var2.name)
  abline(a=lsfit(gen.diff(var1, time), gen.diff(var2, time))$coefficients[1], b=lsfit(gen.diff(var1, time), gen.diff(var2, time))$coefficients[2])
  lines(x=c(0, 0), y=c(-10, 10), col="grey", lty=2)
  lines(x=c(-10, 10), y=c(0, 0), col="grey", lty=2)
  text(x=min(sort(gen.diff(var1, time))), y=max(sort(gen.diff(var2, time))), paste("rho = ", round(cor.test(gen.diff(var1, time), gen.diff(var2, time), method="spearman")$estimate, 2)), adj=c(0, 1))
  plot(var1, var3, cex=0.5, xlab=var1.name, ylab=var3.name)
  abline(a=lsfit(var1, var3)$coefficients[1], b=lsfit(var1, var3)$coefficients[2])
  text(x=min(sort(var1)), y=max(sort(var3)), paste("rho = ", round(cor.test(var1, var3, method="spearman")$estimate, 2)), adj=c(0, 1))
  plot(gen.diff(var1, time), gen.diff(var3, time), cex=0.5, xlab=var1.name, ylab=var3.name)
  abline(a=lsfit(gen.diff(var1, time), gen.diff(var3, time))$coefficients[1], b=lsfit(gen.diff(var1, time), gen.diff(var3, time))$coefficients[2])
  lines(x=c(0, 0), y=c(-10, 10), col="grey", lty=2)
  lines(x=c(-10, 10), y=c(0, 0), col="grey", lty=2)
  text(x=min(sort(gen.diff(var1, time))), y=max(sort(gen.diff(var3, time))), paste("rho = ", round(cor.test(gen.diff(var1, time), gen.diff(var3, time), method="spearman")$estimate, 2)), adj=c(0, 1))
}

ts.diffs.plot <- function(var1.holder, var2.holder, var3.holder, log.var1=FALSE, log.var2=FALSE, log.var3=FALSE, diff.type="fd", mav, time)
{
	var1 <- var2 <- var3 <- vector(mode="numeric")
	for (i in 1:length(var1.holder[1, ])) var1[i] <- mean(sort(var1.holder[, i]))
	for (i in 1:length(var2.holder[1, ])) var2[i] <- mean(sort(var2.holder[, i]))
	for (i in 1:length(var3.holder[1, ])) var3[i] <- mean(sort(var3.holder[, i]))
	if(log.var1 == FALSE) var1 <- as.numeric(gsub(NaN, NA, var1))
	if(log.var1 == TRUE) var1 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, var1)))))
	if(log.var2 == FALSE) var2 <- as.numeric(gsub(NaN, NA, var2))
	if(log.var2 == TRUE) var2 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, var2)))))
	if(log.var3 == FALSE) var3 <- as.numeric(gsub(NaN, NA, var3))
	if(log.var3 == TRUE) var3 <- as.numeric(gsub(-Inf, NA, log(as.numeric(gsub(NaN, NA, var3)))))
	if(diff.type == "fd") {
		var1 <- diff(var1)
		var2 <- diff(var2)
		var3 <- diff(var3)
	}
	if(diff.type == "gd") {
		var1 <- gen.diff(var1, time)
		var2 <- gen.diff(var2, time)
		var3 <- gen.diff(var3, time)
	}
	if(diff.type == "mav") {
		trimtb <- floor(mav/2)
		mav.var1 <- vector(mode="numeric", length=length(var1)-(2*trimtb))
		mav.var2 <- vector(mode="numeric", length=length(var2)-(2*trimtb))
		mav.var3 <- vector(mode="numeric", length=length(var3)-(2*trimtb))
		for (i in (trimtb+1):(length(var1)-trimtb)) mav.var1[(i-trimtb)] <- mean(sort(var1[(i-trimtb):(i+trimtb)]))
		for (i in (trimtb+1):(length(var2)-trimtb)) mav.var2[(i-trimtb)] <- mean(sort(var2[(i-trimtb):(i+trimtb)]))
		for (i in (trimtb+1):(length(var3)-trimtb)) mav.var3[(i-trimtb)] <- mean(sort(var3[(i-trimtb):(i+trimtb)]))
		mav.var1 <- as.numeric(gsub(NaN, NA, mav.var1))
		mav.var2 <- as.numeric(gsub(NaN, NA, mav.var2))
		mav.var3 <- as.numeric(gsub(NaN, NA, mav.var3))
		mav.var1 <- c(rep(NA, trimtb), mav.var1, rep(NA, trimtb))
		mav.var2 <- c(rep(NA, trimtb), mav.var2, rep(NA, trimtb))
		mav.var3 <- c(rep(NA, trimtb), mav.var3, rep(NA, trimtb))
		var1 <- var1-mav.var1
		var2 <- var2-mav.var2
		var3 <- var3-mav.var3
	}
	par(mfrow=c(2, 1))
	if(diff.type == "fd") {
		time <- time[1:(length(time)-1)]+diff(time)/2
		plot(time, var1, xlim=c(max(time), 0), ylim=c(-max(sqrt(sort(c(var1, var2))^2)), max(sqrt(sort(c(var1, var2))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="First diff.")
		points(time, var2, type="l")
	}
	if(diff.type == "gd") {
		time <- time[1:(length(time)-1)]+diff(time)/2
		plot(time, var1, xlim=c(max(time), 0), ylim=c(-max(sqrt(sort(c(var1, var2))^2)), max(sqrt(sort(c(var1, var2))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Gen. diff.")
		points(time, var2, type="l")
	}
	if(diff.type == "mav") {
		plot(time, var1, xlim=c(max(time), 0), ylim=c(-max(sqrt(sort(c(var1, var2))^2)), max(sqrt(sort(c(var1, var2))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Moving average detrended")
		points(time, var2, type="l")
	}
	lines(x=c(max(time)+20, -20), y=c(0, 0), col="grey", lty=2)
	if(diff.type == "fd") plot(time, var1, xlim=c(max(time), 0), ylim=c(-max(sqrt(sort(c(var1, var3))^2)), max(sqrt(sort(c(var1, var3))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="First diff.")
	if(diff.type == "gd") plot(time, var1, xlim=c(max(time), 0), ylim=c(-max(sqrt(sort(c(var1, var3))^2)), max(sqrt(sort(c(var1, var3))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Gen. diff.")
	if(diff.type == "mav") plot(time, var1, xlim=c(max(time), 0), ylim=c(-max(sqrt(sort(c(var1, var3))^2)), max(sqrt(sort(c(var1, var3))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Moving average detrended")
	points(time, var3, type="l")
	lines(x=c(max(time)+20, -20), y=c(0, 0), col="grey", lty=2)
}

parse.nexustotnt <- function(file)
{
	X <- scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
	# Remove trees:
	deleteline <- grep("BEGIN TREES", X)-1
	X <- X[1:deleteline]
	# Replace #NEXUS with xread:
	nexusline <- grep("NEXUS", X)
	X[nexusline] <- "xread"
	# Replace [!...] with '...':
	textline <- grep("[", X, fixed=TRUE)
	X[textline] <- gsub("[!", "'", X[textline], fixed=TRUE)
	X[textline] <- gsub("]", "'", X[textline], fixed=TRUE)
	# Replace DATA block with char, ntax:
	ntaxline <- grep("NTAX", X)
	X[ntaxline] <- gsub("\tDIMENSIONS  NTAX=", "", X[ntaxline])
	X[ntaxline] <- gsub("NCHAR=", "", X[ntaxline])
	X[ntaxline] <- gsub(";", "", X[ntaxline])
	X[ntaxline] <- paste(strsplit(X[ntaxline], " ")[[1]][2], strsplit(X[ntaxline], " ")[[1]][1], sep=" ")
	# Delete now redundant starting lines:
	X <- X[-(3:(ntaxline-1))]
	formatline <- grep("FORMAT", X)
	matrixline <- grep("MATRIX", X)
	semicolonline <- grep(";", X, fixed=TRUE)[grep(TRUE, grep(";", X, fixed=TRUE) > matrixline)][1]
	if(length(grep("(", X[(matrixline+1):(semicolonline-1)], fixed=TRUE)) > 0) X[(matrixline+1):(semicolonline-1)] <- gsub("(", "[", X[(matrixline+1):(semicolonline-1)], fixed=TRUE) # Replace polymorphic ( with [
	if(length(grep(")", X[(matrixline+1):(semicolonline-1)], fixed=TRUE)) > 0) X[(matrixline+1):(semicolonline-1)] <- gsub(")", "]", X[(matrixline+1):(semicolonline-1)], fixed=TRUE) # Replace polymorphic ) with ]
	missingchar <- strsplit(strsplit(X[formatline], "MISSING=")[[1]][2], "")[[1]][1] # Find character used for missing states
	gapchar <- strsplit(strsplit(X[formatline], "GAP=")[[1]][2], "")[[1]][1] # Find character used for gaps
	if (missingchar != "?") X[(matrixline+1):(semicolonline-1)] <- gsub(missingchar, "?", X[(matrixline+1):(semicolonline-1)], fixed=TRUE) # Replace with ? if not ?
	if (gapchar != "-") X[(matrixline+1):(semicolonline-1)] <- gsub(gapchar, "-", X[(matrixline+1):(semicolonline-1)], fixed=TRUE) # Replace with - if not -
	X <- X[-matrixline]
	X <- X[-formatline]
	# Re-find end of matrix:
	semicolonline <- grep(";", X, fixed=TRUE)[1]
	# Search for ordered character list:
	if (length(grep("TYPESET", X)) > 0) {
		orderingline <- grep("TYPESET", X)
		X[orderingline] <- gsub("\tTYPESET * UNTITLED  = ", "", X[orderingline], fixed=TRUE) # Clean up
		X[orderingline] <- gsub(";", "", X[orderingline], fixed=TRUE) # Clean up
		unordered <- strsplit(X[orderingline], ", ")[[1]][1] # Get unordered list
		ordered <- strsplit(X[orderingline], ", ")[[1]][2] # Get ordered list
		unordered <- gsub("unord: ", "", unordered, fixed=TRUE) # Clean up
		ordered <- gsub("ord: ", "", ordered, fixed=TRUE) # Clean up
		if(length(grep("-", unordered)) > 0) {
			while(length(grep("-", unordered)) > 0) {
				unordered <- strsplit(unordered, " ")[[1]]
				unordered[grep("-", unordered)[1]] <- paste(c(strsplit(unordered[grep("-", unordered)[1]], "-")[[1]][1]:strsplit(unordered[grep("-", unordered)[1]], "-")[[1]][2]), collapse=" ")
				unordered <- paste(unordered, collapse=" ")
			}
		}
		if(length(grep("-", ordered)) > 0) {
			while(length(grep("-", ordered)) > 0) {
				ordered <- strsplit(ordered, " ")[[1]]
				ordered[grep("-", ordered)[1]] <- paste(c(strsplit(ordered[grep("-", ordered)[1]], "-")[[1]][1]:strsplit(ordered[grep("-", ordered)[1]], "-")[[1]][2]), collapse=" ")
				ordered <- paste(ordered, collapse=" ")
			}
		}
		orderline <- paste("ccode +", ordered, " -", unordered, ";", sep="", collapse="")
		endline <- c(orderline, "proc/;")
	} else {
		orderline <- paste("ccode -", paste(1:strsplit(X[3], " ")[[1]][1], collapse=" "), ";", sep="", collapse="")
		endline <- c(orderline, "proc/;")
	}
	X <- X[-((semicolonline+1):length(X))]
	X <- c(X, endline)
	return(X)
}

path.lengths <- function(phy)
{
	ntips <- length(phy$tip.label)
	pathlengths <- vector(mode="numeric")
	for (i in 1:ntips) {
		taxon <- i
		pathedges <- vector(mode="numeric")
		pathedges[1] <- grep(TRUE, phy$edge[, 2] == i)
		while (phy$edge[pathedges[length(pathedges)], 1] != (ntips+1)) {
			i <- grep(TRUE, phy$edge[, 2] == phy$edge[pathedges[length(pathedges)], 1])
			pathedges <- c(pathedges, i)
		}
		pathlengths[taxon] <- sum(phy$edge.length[pathedges])
	}
	names(pathlengths) <- phy$tip.label
	return(pathlengths)
}

pat.dist.phylo <- function(tree, comp, nchar)
{
	require(ape)
	nt <- Ntip(tree)
	cchar <- nchar-comp[tree$tip.label, 1]
	for (i in 1:length(tree$edge[, 1])) {
		if (tree$edge[i, 2] <= max(nt))  tree$edge.length[i] <- tree$edge.length[i]/cchar[tree$edge[i, 2]]
		if (tree$edge[i, 2] > max(nt))  tree$edge.length[i] <- tree$edge.length[i]/nchar
	}
	return(tree)
}

SCM <- function(names, ages, tbins)
{
	out <- matrix(NA, nrow=length(unique(names)), ncol=length(tbins)-1)
	rownames(out) <- sort(unique(names))
	for (i in 1:length(out[, 1])) {
		dates <- sort(ages[grep(TRUE, names == rownames(out)[i])])
		for(j in 1:length(dates)) {
			out[i, intersect(grep(TRUE, dates[j] <= tbins[1:(length(tbins)-1)]), grep(TRUE, dates[j] > tbins[2:length(tbins)]))] <- 1
		}
		start <- min(grep(TRUE, out[i, ] == 1))
		stop <- max(grep(TRUE, out[i, ] == 1))
		out[i, (start+grep(TRUE, is.na(out[i, start:stop]))-1)] <- 0
	}
	out.adj <- out
	for (i in length(out.adj[, 1]):1) {
		out.adj[i, min(grep(TRUE, out.adj[i, ] == 1))] <- NA
		out.adj[i, max(grep(TRUE, out[i, ] == 1))] <- NA
		if(length(sort(out.adj[i, ])) == 0) out.adj <- out.adj[-i, ]
	}
	SCM.taxa.adj <- SCM.tbins.adj <- SCM.taxa <- SCM.tbins <- vector(mode="numeric")
	for(i in 1:length(out[1, ])) {
		SCM.tbins[i] <- mean(sort(out[, i]))
	}
	for(i in 1:length(out[, 1])) {
		SCM.taxa[i] <- mean(sort(out[i, ]))
	}
	for(i in 1:length(out.adj[1, ])) {
		SCM.tbins.adj[i] <- mean(sort(out.adj[, i]))
	}
	for(i in 1:length(out.adj[, 1])) {
		SCM.taxa.adj[i] <- mean(sort(out.adj[i, ]))
	}
	names(SCM.taxa) <- rownames(out)
	names(SCM.taxa.adj) <- rownames(out.adj)
	tbin.midpoints <- ((tbins[1:(length(tbins)-1)]-tbins[2:length(tbins)])/2)+tbins[2:length(tbins)]
	result <- list(out, SCM.tbins, SCM.taxa, SCM.tbins.adj, SCM.taxa.adj, tbin.midpoints)
	names(result) <- c("out", "SCM.tbins", "SCM.taxa", "SCM.tbins.adj", "SCM.taxa.adj", "tbin.midpoints")
	return(result)
}

make.phan <- function(x, y, y_label, quartz=TRUE, main="")
{
	y.max <- max(sort(y))
	y.min <- 0
	young.bin <- -8 # Offset for youngest bin label...allows for age < 0 to center label
	plot.min <- y.min-.03*(y.max-y.min)
	seg.min <- y.min-.07*(y.max-y.min)
	text.min <- y.min-.036*(y.max-y.min)
	time.boundaries <- c(23.03, 65.5, 145.5, 199.6, 251, 299, 359.2, 416, 443.7)
	interval.names <- c("Ng", "Pg", "K", "J", "Tr", "P", "C", "D", "S")
	interval.midpoint <- time.boundaries-diff(c(0, time.boundaries))/2
	plot(1, 1, xlim=c(max(time.boundaries), 0), ylim=c(plot.min, y.max), type="n", xlab="Geologic time (Ma)", ylab=y_label, main=main)
	abline(h=y.min)
	segments(c(time.boundaries, 0), y.min, c(time.boundaries, 0), seg.min)
	text(interval.midpoint, text.min, labels=interval.names)
	lines(x, y)
}

randomisation.phylo <- function(tree, permutations)
{
	nchang <- sum(tree$edge.length)
	permat <- matrix(0, nrow=length(tree$edge.length), ncol=permutations)
	brk <- vector(mode="numeric", length=length(tree$edge.length))
	brk[1] <- 1/length(tree$edge.length)
	for (i in 2:length(brk)) brk[i] <- brk[i-1]+1/length(tree$edge.length)
	for (i in 1:permutations) {
		randno <- runif(nchang, min=0, max=1)
		randno <- sort(randno)
		for (j in 1:length(tree$edge.length)) {
			permat[j, i] <- length(grep(TRUE, randno <= brk[j]))
			randno <- randno[-grep(TRUE, randno <= brk[j])]
		}
	}
	probs <- vector(mode="numeric")
	for (i in 1:length(tree$edge.length)) {
		probs[i] <- wilcox.test(tree$edge.length[i], permat[i, ], alternative="greater")$p.value
	}
	probs
}

randombranch.phylo <- function(tree, ttree, comp, nchar, permutations)
{
	charcorr <- c((comp[ttree$tip.label, ]/nchar), rep(1, Nnode(ttree)))
	ttree$edge.length <- charcorr[tree$edge[, 2]]*ttree$edge.length
	nchang <- sum(tree$edge.length)
	permat <- matrix(0, nrow=length(tree$edge.length), ncol=permutations)
	brk <- vector(mode="numeric", length=length(ttree$edge.length))
	brk[1] <- ttree$edge.length[1]/length(ttree$edge.length)
	for (i in 2:length(brk)) brk[i] <- brk[i-1]+ttree$edge.length[i]/length(ttree$edge.length)
	brk <- brk/(sum(ttree$edge.length)/length(ttree$edge.length))
	for (i in 1:permutations) {
		randno <- runif(nchang, min=0, max=1)
		randno <- sort(randno)
		for (j in 1:length(tree$edge.length)) {
			breaker <- length(grep(TRUE, randno <= brk[j]))
			permat[j, i] <- breaker
			if (breaker > 0) randno <- randno[-(1:breaker)]
		}
	}
	sigs <- vector(mode="numeric")
	for (i in 1:length(tree$edge.length)) {
		ifelse(tree$edge.length[i] > sort(permat[i, ])[ceiling(0.95*length(permat[i, ]))], sigs[i] <- 1, sigs[i] <- 0)
	}
	sigs
}

get.branches <- function(node, froms, tos)
{
	branches <- vector(mode="numeric")
	for (i in 1:length(froms)) {
		for (j in 1:length(node)) {
			branches[length(branches)+1:length(branches)+2] <- grep(node[j], froms)
			branches <- sort(unique(branches))
		}
		node <- unique(c(node, froms[sort(match(tos[branches], froms))]))
	}
	branches <- vector(mode="numeric")
	for (i in 1:length(node)) {
		branches[length(branches)+1:length(branches)+2] <- grep(node[i], froms)
		branches <- sort(unique(branches))
	}
	branches
}

broken.stick <- function(timeseries, time)
{
	require(stats)
	require(nlme)
	require(paleoTS)
	if(length(grep("1 1", paste(diff(grep(TRUE, diff(timeseries) == 0)), collapse=" "))) > 0) { # Case if 4 values in a row are the same
		print("Some bins deleted due to consecutive equal values (confounding model comparison)")
		deletes <- vector(mode="numeric")
		for(i in 2:length(timeseries)) {
			if(timeseries[i-1] == timeseries[i]) deletes[length(deletes)+1] <- i
		}
		timeseries <- timeseries[-deletes]
		time <- time[-deletes]
	}
	tslength <- length(timeseries)
	onestick <- lm(timeseries~time) # First fit single stick model
	onestickAICc <- AICc(onestick, tslength)
	if(length(timeseries) < 15 && length(timeseries) >= 10) print("Time series too short (too few bins) for three stick model") # Case if broken stick pointless
	if(length(timeseries) < 10) print("Time series too short (too few bins) for two or three stick model") # Case if broken stick pointless
	if(length(timeseries) >= 10) { # Case if two stick doable
		splits <- vector(length=tslength-9, mode="numeric")
		for (i in 1:(tslength-9)) splits[i] <- (4+i) # Find all splits of at least five bins in length
		AICcsplits <- splits # Set up vector to store AICc of two stick models
		for(i in 1:length(splits)) { # Fill the above
			firstsplit <- 1:splits[i]
			secondsplit <- (splits[i]+1):tslength
			AICcsplits[i] <- AICc(lm(timeseries[firstsplit]~time[firstsplit]), length(firstsplit))+AICc(lm(timeseries[secondsplit]~time[secondsplit]), length(secondsplit))
		}
		AICcsplits <- as.numeric(AICcsplits)
		twostickAICc <- min(AICcsplits)
		besttwostick <- splits[grep(TRUE, AICcsplits == min(AICcsplits))] # Establish split with lowest combined AICc
		twostickone <- lm(timeseries[1:besttwostick]~time[1:besttwostick]) # Fit first stick of two stick model
		twosticktwo <- lm(timeseries[(besttwostick+1):tslength]~time[(besttwostick+1):tslength]) # Fit second stick of two stick model
	}
	if(length(timeseries) >= 15) { # Case if three stick doable
		splits <- vector(mode="character")
		for(i in 5:(tslength-10)) {
			for(j in (i+5):(tslength-5)) {
				splits[(length(splits)+1)] <- paste(i, ":", j, sep="", collapse="")
			}
		}
		AICcsplits <- splits
		for(i in 1:length(splits)) {
			AICcsplits[i] <- AICc(lm(timeseries[1:strsplit(splits[i], ":")[[1]][1]]~time[1:strsplit(splits[i], ":")[[1]][1]]), length(1:strsplit(splits[i], ":")[[1]][1]))+AICc(lm(timeseries[(as.numeric(strsplit(splits[i], ":")[[1]][1])+1):strsplit(splits[i], ":")[[1]][2]]~time[(as.numeric(strsplit(splits[i], ":")[[1]][1])+1):strsplit(splits[i], ":")[[1]][2]]), length((as.numeric(strsplit(splits[i], ":")[[1]][1])+1):strsplit(splits[i], ":")[[1]][2]))+AICc(lm(timeseries[(as.numeric(strsplit(splits[i], ":")[[1]][2])+1):tslength]~time[(as.numeric(strsplit(splits[i], ":")[[1]][2])+1):tslength]), length((as.numeric(strsplit(splits[i], ":")[[1]][2])+1):tslength))
		}
		AICcsplits <- as.numeric(AICcsplits)
		threestickAICc <- min(AICcsplits)
		bestthreestick <- splits[grep(TRUE, AICcsplits == min(AICcsplits))] # Establish split with lowest combined AICc
		threestickone <- lm(timeseries[1:as.numeric(strsplit(bestthreestick, ":")[[1]][1])]~time[1:as.numeric(strsplit(bestthreestick, ":")[[1]][1])]) # Fit first stick of three stick model
		threesticktwo <- lm(timeseries[(as.numeric(strsplit(bestthreestick, ":")[[1]][1])+1):as.numeric(strsplit(bestthreestick, ":")[[1]][2])]~time[(as.numeric(strsplit(bestthreestick, ":")[[1]][1])+1):as.numeric(strsplit(bestthreestick, ":")[[1]][2])]) # Fit second stick of three stick model
		threestickthree <- lm(timeseries[(as.numeric(strsplit(bestthreestick, ":")[[1]][2])+1):tslength]~time[(as.numeric(strsplit(bestthreestick, ":")[[1]][2])+1):tslength]) # Fit second stick of three stick model
	}
	modelcomp <- cbind(c(onestickAICc, twostickAICc, threestickAICc), round(akaike.wts(c(onestickAICc, twostickAICc, threestickAICc)), 3))
	rownames(modelcomp) <- c("One stick", "Two stick", "Three stick")
	colnames(modelcomp) <- c("AICc", "Akaike wt")
	result <- list(timeseries, time, onestick, twostickone, twosticktwo, threestickone, threesticktwo, threestickthree, besttwostick, bestthreestick, modelcomp)
	names(result) <- c("y.values", "x.values", "onestick.model", "twostick.model1", "twostick.model2", "threestick.model1", "threestick.model2", "threestick.model3", "twostick.split", "threestick.split", "model.comparison")
	return(result)
}

mars <- function(x, y)
{
	require(earth)
	AICs <- vector(mode="numeric")
	for(i in 3:21) { # Find optimal number of knots using AIC
		marsout <- earth(x, y, nk=i)
		cuts <- sort(x[match(sort(unique(marsout$cuts)), x)])
		rss <- marsout$rss
		AICs[(i-2)] <- (2*(length(cuts)+1))+(length(x)*log(rss))
	}
	nk <- max(grep(TRUE, min(AICs) == AICs))+2
	marsout <- earth(x, y, nk=nk)
	prediction <- y-marsout$residuals
	cuts <- sort(x[match(sort(unique(marsout$cuts)), x)])
	residuals <- marsout$residuals
	rss <- marsout$rss
	slopes <- vector(mode="numeric", length=length(cuts)+1)
	breaks <- sort(unique(c(1, match(cuts, x), length(x))))
	for(i in 2:length(breaks)) slopes[(i-1)] <- prediction[breaks[i]]-prediction[breaks[(i-1)]]
	slopes <- round(slopes, 4)
	AIC <- (2*(length(cuts)+1))+(length(x)*log(rss))
	linearRSS <- sum((y-((lsfit(x, y)$coefficients[2]*x)+lsfit(x, y)$coefficients[1]))^2)
	linearAIC <- 4+(length(x)*log(linearRSS))
	ifelse(linearAIC <= AIC, best.model <- "Linear", best.model <- paste("MARS with ", length(cuts), " hinges", sep="", collapse=""))
	result <- list(prediction, cuts, residuals, rss, slopes, AIC, linearRSS, linearAIC, best.model)
	names(result) <- c("prediction", "cuts", "residuals", "rss", "slopes", "AIC", "linearRSS", "linearAIC", "best.model")
	return(result)
}




# sqs version 2.0 by John Alroy
# performs shareholder quorum subsampling on an array of specimen counts
# can be used to perform classical rarefaction instead of SQS
# written 29 July 2010; version 2.0 completed 14 February 2011
# changes in version 2.0: improved subsampling algorithm; including the dominant 
#  taxon is now the default; improved reporting of errors and basic statistics
# warning: do not use this program with taxonomic occurrence data drawn from
#  multiple published references because it is not designed to count
#  single-reference taxa or adjust for long taxonomic lists
# warning: version 1.0 yields estimates that are downwards-biased when q < 0.6
#  and abundance distributions are highly uneven
#
# Modified by GTL 22/03/11 to include single publication occurrence correction
sqs <- function(ab, q, trials, method, dominant, p.1)	{
	
	params <- array(data=NA, dim=0, dimnames=c("raw richness"))
	if (missing(trials))	{
		trials <- 100
	}
	if (missing(method))	{
		method <- ""
	} else if (method != "" && method != "rarefaction" && method != "CR")	{
		return(print('If the method is rarefaction enter method="rarefaction" or "CR"', quote=F))
	}
	if ((q <= 0 || q >= 1) && method != "rarefaction" && method != "CR")	{
		return(print("If the method is SQS the quota must be greater than zero and less than one", quote=F))
	} else if (q < 1 && (method == "rarefaction" || method == "CR"))	{
		return(print("If the method is rarefaction the quota must be an integer", quote=F))
	}
	if (missing(dominant))	{
		dominant <- 0
	} else if (dominant != "" && dominant != "exclude" && dominant != "no")	{
		return(print('To exclude the dominant taxon, enter dominant="exclude" or "no"', quote=F))
	}
	
	# compute basic statistics
	specimens <- sum(ab)
	singletons <- 0
	doubletons <- 0
	highest <- 0
	for (i in 1:length(ab))	{
		if (ab[i] == 1)	{
			singletons <- singletons + 1
		} else if (ab[i] == 2)	{
			doubletons <- doubletons + 1
		}
		if (ab[i] > highest)	{
			highest <- ab[i]
			mostfrequent <- i
		}
	}
	
	u <- 1 - singletons / specimens
	
	# GTL modification starts
	if (missing(p.1) && dominant == "exclude") {
		u <- 1 - singletons/(specimens - highest)
	}
	if (missing(p.1) && dominant == "no") {
		u <- 1 - singletons/(specimens - highest)
	}
	if (!missing(p.1)) {
		u <- (sum(ab) - p.1)/sum(ab)
	}
	# GTL modification ends

	if (u == 0)	{
		return(print("Coverage is zero because all taxa are singletons", quote=F))
	}
	
	# compute raw taxon frequencies (temporarily)
	freq <- ab / specimens
	
	# standard recursive equation for Fishers alpha
	alpha <- 10
	oldalpha <- 0
	while (abs(alpha - oldalpha) > 0.0000001)	{
		oldalpha <- alpha
		alpha <- length(ab) / log(1 + specimens/alpha)
	}

	params["raw richness"] <- length(ab)
	params["Good's u"] <- u
	params["subsampled richness"] <- NA
	params["subsampled u"] <- NA
	params["Chao 1"] <- length(ab) + singletons**2/(2* doubletons)
	params["subsampled Chao 1"] <- NA
	# governing parameter of the geometric series distribution
	params["k"] <- abs(lm(log(sort(freq)) ~ c(1:length(freq)))$coefficients[2])
	params["Fisher's alpha"] <- alpha
	params["Shannon's H"] <- -1 * sum(freq * log(freq))
	params["Hurlbert's PIE"] <- (1 - sum(freq**2)) * length(ab) / (length(ab) - 1)
	params["dominance"] <- highest / specimens
	params["specimens"] <- specimens
	params["singletons"] <- singletons
	params["doubletons"] <- doubletons
	params["specimens drawn"] <- 0

	if (dominant != "exclude" && dominant != "no")	{
		highest <- 0
		mostfrequent <- 0
	}

	# return if the quorum target is higher than overall coverage
	if ((q > u && method != "rarefaction" && method != "CR") || (q >= sum(ab)))	{
		return(params)
	}
	# return if the rarefaction quota is equal to or higher than the
	#  specimen count
	if (method == "rarefaction" && q >= specimens - highest)	{
		return(params)
	}

	# compute adjusted taxon frequencies
	freq <- ab * u / (specimens - highest)

	# create an array in which each cell corresponds to one specimen
	ids <- array()
	n <- 0
	for (i in 1:length(ab))	{
		for (j in 1:ab[i])	{
			n <- n + 1
			ids[n] <- i
		}
	}

	# subsampling trial loop
	# s will be the subsampled taxon count
	s <- array(rep(0, trials))
	subsingle <- array(rep(0, trials))
	subdouble <- array(rep(0, trials))
	subchao <- array(rep(0, trials))
	mostfrequentdrawn <- 0
	for (trial in 1:trials)	{
		pool <- ids
		left <- length(pool)
		seen <-  array(data=rep(0, length(ab)))
		subfreq <- array(rep(0, length(ab)))
		if (method != "rarefaction" && method != "CR")	{
			udrawn <- 0
			while (udrawn < q)	{
				# draw a specimen
				x <- floor(runif(1, min=1, max=left+1))
				# add to frequency and taxon sums if species has
				#  not been drawn previously
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if (seen[pool[x]] == 0)	{
					if (pool[x] != mostfrequent)	{
						udrawn <- udrawn + freq[pool[x]]
					}
					seen[pool[x]] <- 1
					# randomly throw back some draws that put the sum over q
					#  (improved algorithm added in version 2.0)
					plus <- 1
					if (udrawn > q)	{
						plus <-  1 - q / udrawn
					}
					if (runif(1) <= plus)	{
						s[trial] <- s[trial] + 1
					} else	{
						subfreq[pool[x]] <- subfreq[pool[x]] - 1
					}
				}
				# decrease pool of specimens not yet drawn
				pool[x] <- pool[left]
				left <- left - 1
			}
		} else	{
			i <- 0
			draws <- 0
			while (i < q)	{
				draws <- draws + 1
				x <- floor(runif(1, min=1, max=length(ids)-draws+2))
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if (pool[x] != mostfrequent)	{
					i <- i + 1
				}
				if (seen[pool[x]] == 0)	{
					seen[pool[x]] <- 1
					s[trial] <- s[trial] + 1
				}
				pool[x] <- pool[length(ids)-draws+1]
			}
		}
		for (i in 1:length(ab))	{
			if (subfreq[i] == 1 && i != mostfrequent)	{
				subsingle[trial] <- subsingle[trial] + 1
			} else if (subfreq[i] == 2 && i != mostfrequent)	{
				subdouble[trial] <- subdouble[trial] + 1
			}
		}
		if (subsingle[trial] > 0 && subdouble[trial] > 0)	{
			subchao[trial] <- s[trial] + subsingle[trial]**2/(2*subdouble[trial])
		} else	{
			subchao[trial] <- s[trial]
		}
		params["specimens drawn"] <- params["specimens drawn"] + sum(subfreq)
		if (mostfrequent != 0)	{
			mostfrequentdrawn <- mostfrequentdrawn + subfreq[mostfrequent]
		}
	}
	params["specimens drawn"] <- params["specimens drawn"] / trials
	# compute vector of non-zero counts
	options(warn=-1)
	s2 <- sort(sqrt(s-1))^2+1
	options(warn=0)
	# compute geometric mean
	params["subsampled richness"] <- exp(mean(log(s2))) * length(s2)/length(s)
	# use of arithmetic means to compute Goods u is adequate
	mostfrequentdrawn <- mostfrequentdrawn / trials
	params["subsampled u"] <- 1 - mean(subsingle) / (params["specimens drawn"] - mostfrequentdrawn)
	params["subsampled Chao 1"] <- exp(mean(log(subchao)))
	return(params)
	
}

# sqs version 1.0 by John Alroy
# performs shareholder quorum subsampling on an array of specimen counts
# set method="rarefaction" or "CR" to perform classical rarefaction instead of SQS
# written 29 July 2010
# Modified by GTL 14/12/10 to include single publication occurrence correction
#sqs <- function(ab, q, trials, method, dominant, p.1)
#{
#	params <- array(data=NA, dim=0, dimnames=c("raw richness"))
#	if (missing(trials))  {
#		trials <- 100
#	}
#	if (missing(method))  {
#		method <- ""
#	}
#	if ((q <= 0 || q >= 1) && method != "rarefaction" && method != "CR")  {
#		print("If the method is SQS the quota must be greater than zero and less than one")
#		return(params)
#	} else if (q < 1 && (method == "rarefaction" || method == "CR"))  {
#		print("If the method is rarefaction the quota must be an integer")
#		return(params)
#	}
#	if (missing(dominant))  {
#		dominant <- 0
#	}
#	
#	# compute basic statistics
#	specimens <- sum(ab)
#	singletons <- 0
#	doubletons <- 0
#	highest <- 0
#	for (i in 1:length(ab))  {
#		if (ab[i] == 1)  {
#			singletons <- singletons + 1
#		} else if (ab[i] == 2)  {
#			doubletons <- doubletons + 1
#		}
#		if (ab[i] > highest)  {
#			highest <- ab[i]
#			mostfrequent <- i
#		}
#	}
#	# exclude dominant taxon unless told to include it
#	if (dominant == "include")  {
#		highest <- 0
#	}
#	
#	if (missing(p.1)) {
#		u <- 1 - singletons/(specimens - highest)
#	}
#	if (!missing(p.1)) {
#		u <- (sum(ab) - p.1)/sum(ab)
#	}
#	if (u == 0)  {
#		print("Coverage is zero because all taxa are singletons")
#		return(params)
#	}
#	
#	# compute raw taxon frequencies (temporarily)
#	freq <- ab / specimens
#	
#	# standard recursive equation for Fishers alpha
#	alpha <- 10
#	oldalpha <- 0
#	while (abs(alpha - oldalpha) > 0.0000001)  {
#		oldalpha <- alpha
#		alpha <- length(ab) / log(1 + specimens/alpha)
#	}
#	
#	params["raw richness"] <- length(ab)
#	params["Good's u"] <- u
#	params["subsampled richness"] <- NA
#	params["Chao 1"] <- length(ab) + singletons**2/(2* doubletons)
#	params["subsampled Chao 1"] <- NA
#	# governing parameter of the geometric series distribution
#	params["k"] <- abs(lm(log(sort(freq)) ~ c(1:length(freq)))$coefficients[2])
#	params["Fisher's alpha"] <- alpha
#	params["Shannon's H"] <- -1 * sum(freq * log(freq))
#	params["Hurlbert's PIE"] <- (1 - sum(freq**2)) * length(ab) / (length(ab) - 1)
#	params["dominance"] <- highest / specimens
#	params["specimens"] <- specimens
#	params["singletons"] <- singletons
#	params["doubletons"] <- doubletons
#	params["specimens drawn"] <- 0
#	
#	# return if the quorum target is higher than overall coverage
#	if ((q > u && method != "rarefaction" && method != "CR") || (q >= sum(ab)))  {
#		return(params)
#	}
#	
#	# compute adjusted taxon frequencies
#	freq <- ab * u / (specimens - highest)
#	
#	# create an array in which each cell corresponds to one specimen
#	ids <- array()
#	n <- 0
#	for (i in 1:length(ab))  {
#		for (j in 1:ab[i])  {
#			n <- n + 1
#			ids[n] <- i
#		}
#	}
#	
#	# subsampling trial loop
#	# s will be the subsampled taxon count
#	s <- array(rep(0, trials))
#	subchao <- array(rep(0, trials))
#	for (trial in 1:trials)  {
#		pool <- ids
#		left <- length(pool)
#		seen <- array(data=rep(0, length(ab)))
#		subfreq <- array(rep(0, length(ab)))
#		
#		if (method != "rarefaction" && method != "CR")  {
#			sumfreq <- 0
#			under <- q
#			while (sumfreq < q)  {
#				# draw a specimen
#				x <- floor(runif(1, min=1, max=left+1))
#				subfreq[pool[x]] <- subfreq[pool[x]] + 1
#				# add to frequency and taxon sums if species has
#				#  not been drawn previously
#				if (seen[pool[x]] == 0)  {
#					if (pool[x] != mostfrequent || dominant == "include")  {
#						sumfreq <- sumfreq + freq[pool[x]]
#					}
#					seen[pool[x]] <- 1
#					# count the taxon putting the sum over the target
#					#  only if the resulting overshoot would be less
#					#  than the undershoot created by not counting it
#					if ((sumfreq >= q && sumfreq - q < under) || sumfreq < q)  {
#						s[trial] <- s[trial] + 1
#					} else  {
#						subfreq[pool[x]] <- subfreq[pool[x]] - 1
#					}
#					under <- q - sumfreq
#				}
#				# decrease pool of specimens not yet drawn
#				pool[x] <- pool[left]
#				left <- left - 1
#			}
#		} else  {
#			for (i in 1:q)  {
#				x <- floor(runif(1, min=1, max=length(ids)-i+2))
#				subfreq[pool[x]] <- subfreq[pool[x]] + 1
#				if (seen[pool[x]] == 0)  {
#					seen[pool[x]] <- 1
#					s[trial] <- s[trial] + 1
#				}
#				pool[x] <- pool[length(ids)-i+1]
#			}
#		}
#		subsingle <- 0
#		subdouble <- 0
#		for (i in 1:length(ab))  {
#			if (subfreq[i] == 1)  {
#				subsingle <- subsingle + 1
#			} else if (subfreq[i] == 2)  {
#				subdouble <- subdouble + 1
#			}
#		}
#		if (subsingle > 0 && subdouble > 0)  {
#			subchao[trial] <- s[trial] + subsingle**2/(2*subdouble)
#		} else  {
#			subchao[trial] <- s[trial]
#		}
#		params["specimens drawn"] <- params["specimens drawn"] + sum(subfreq)
#	}
#	
#	params["subsampled richness"] <- exp(mean(log(s)))
#	params["subsampled Chao 1"] <- exp(mean(log(subchao)))
#	params["specimens drawn"] <- params["specimens drawn"] / trials
#	return(params)
#}

gen.diff <- function(x, time)
{
	#if(cor.test(time, x)$p.value > 0.05) print("Warning: variables not significantly correlated, generalised differencing not recommended")
	dt <- x-((lsfit(time, x)$coefficients[2]*time)+lsfit(time, x)$coefficients[1])
	m <- lsfit(dt[1:(length(dt)-1)], dt[2:length(dt)])$coefficients[2]
	gendiffs <- dt[1:(length(dt)-1)]-(dt[2:length(dt)]*m)
	gendiffs
}

getthedamndatain <- function(file, sep="\t", header=F)
{
	X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
	data <- matrix(nrow=length(X), ncol=length(strsplit(X[1], sep)[[1]]))
	for(i in 1:length(X)) {
		data[i, 1:length(strsplit(X[1], sep)[[1]])] <- strsplit(X[i], sep)[[1]]
	}
	if(header == TRUE) {
		colnames(data) <- data[1, ]
		data <- data[-1, ]
	}
	return(data)
}

my.ccf <- function(x, y, xlab, ylab) {
	completeforboth <- intersect(grep(TRUE, !is.na(x)), grep(TRUE, !is.na(y)))
	differences <- diff(completeforboth)
	if(length(grep(TRUE, differences != 1)) > 0) { # Does gap exist?
		for(i in 1:length(grep(TRUE, differences != 1))) {
			if(i == 1) {
				split <- completeforboth[1:grep(TRUE, differences != 1)[i]]
				if(length(grep(TRUE, differences != 1)) == 1) {
					split.a <- completeforboth[(grep(TRUE, differences != 1)[i]+1):length(completeforboth)]
					if(length(split.a) > length(split)) split <- split.a
				}
			}
			if(i > 1) {
				split.a <- completeforboth[(grep(TRUE, differences != 1)[(i-1)]+1):grep(TRUE, differences != 1)[i]]
				if(length(split.a) > length(split)) split <- split.a
			}
			if(i == length(grep(TRUE, differences != 1)) && i > 1) {
				split.a <- completeforboth[(grep(TRUE, differences != 1)[i]+1):length(completeforboth)]
				if(length(split.a) > length(split)) split <- split.a
			}
		}
	}
	if(length(grep(TRUE, differences != 1)) == 0) split <- completeforboth
	x <- x[split]
	y <- y[split]
	plot.title <- paste(xlab, "vs.", ylab, "\n(Continuous from bin", split[1], "to", split[length(split)], ")")
	ccf(y, x, main=plot.title)
}

get.anc.states <- function(nexus.matrix, tree) {
    require(ape)
    anc.lik.matrix <- matrix(nrow=Nnode(tree), ncol=length(nexus.matrix$matrix[1, ])) # Create matrix to record ancestral state estimates
    rownames(anc.lik.matrix) <- c((Ntip(tree)+1):(Ntip(tree)+Nnode(tree))) # Label matrix to record ancestral state estimates
    for(i in 1:length(nexus.matrix$matrix[1, ])) { # Cycle through characters:
        tipstogo <- rownames(nexus.matrix$matrix)[grep(TRUE, is.na(nexus.matrix$matrix[, i]))] # Find taxa with missing data:
        if(length(tipstogo) > 0) chartree <- drop.tip(tree, tipstogo) # Remove tips with missing data
        if(length(tipstogo) == 0) chartree <- tree # If no missing data use whole tree    
        tipvals <- nexus.matrix$matrix[chartree$tip.label, i]
        # Create state change matrix (important for ordered/unordered characters:
        minval <- as.numeric(nexus.matrix$min.vals[i])
        maxval <- as.numeric(nexus.matrix$max.vals[i])
        if(maxval-minval == 1) mymodel <- matrix(c(minval, maxval, maxval, minval), 2) # Case for binary character (ordering irrelevant)
        if(maxval-minval > 1 && nexus.matrix$ordering[i] == "unord") { # Case for unordered multistate character
            mymodel <- matrix(1, ncol=(maxval-minval)+1, nrow=(maxval-minval)+1)
            for(j in 1:length(mymodel[1, ])) mymodel[j, j] <- 0
        }
        if(maxval-minval > 1 && nexus.matrix$ordering[i] == "ord") { # Case for ordered multistate character
            mymodel <- matrix(0, ncol=(maxval-minval)+1, nrow=(maxval-minval)+1)
            for(j in 1:length(mymodel[1, ])) {
                for(k in 1:length(mymodel[1, ])) {
                    mymodel[j, k] <- sqrt((diff(c(j, k)))^2)
                }
            }
        }
        # Ascertain if polymorphisms are present and if so compute all possible prmutations
        if(length(grep("&", tipvals)) > 0) {
            permutation.counts <- vector(mode="numeric") # Vector to store permutation numbers for each taxon (i.e. number of possible states)
            for(j in 1:length(tipvals)) permutation.counts[j] <- length(strsplit(as.character(tipvals[j]), "&")[[1]]) # Get number of permutations for each taxon
            permutation.count <- prod(permutation.counts) # Get product of permutation counts (i.e. total number of permutations)
            permutations <- matrix(nrow=length(tipvals), ncol=permutation.count) # Make permutations matrix
            permutations[grep(TRUE, permutation.counts == 1), 1:permutation.count] <- as.numeric(tipvals[grep(TRUE, permutation.counts == 1)]) # Fill in data for non-polymorphic characters
            phase.count <- 1 # Permutation phase
            for(j in grep("&", tipvals)) {
                perm <- as.numeric(strsplit(as.character(tipvals[j]), "&")[[1]])
                if(phase.count == 1) permutations[j, ] <- perm
                if(phase.count > 1) {
                    perm <- sort(rep(perm, phase.count))
                    permutations[j, ] <- perm
                }
                phase.count <- length(perm)
            }
            rownames(permutations) <- names(tipvals) # Make sure tip value names are carried over
            # Ancestor state estimate for first permutation:
            if(length(unique(sort(permutations[, 1]))) > 1) { # If there is more than one state
                if(length(mymodel[, 1]) > length(sort(unique(as.numeric(permutations[, 1]))))) { # Case if model has larger spread than actual tipvalues (e.g. spans 0-5 when 2 is not recorded)
                    anc.lik <- ace(permutations[, 1], chartree, type="discrete", model=mymodel[sort(unique(as.numeric(permutations[, 1])))+1, sort(unique(as.numeric(permutations[, 1])))+1])$lik.anc
                }
                if(length(mymodel[, 1]) <= length(sort(unique(as.numeric(permutations[, 1]))))) { # Case if model has same spread as tip values (theoretically this should always be true, but in reality not so much)
                    anc.lik <- ace(permutations[, 1], chartree, type="discrete", model=mymodel[])$lik.anc
                }
            }
            if(length(unique(sort(permutations[, 1]))) == 1) {
                anc.lik <- matrix(0, nrow=length(permutations[, 1])-1, ncol=max(c((unique(sort(as.numeric(permutations[, 1])))+1), 2)))
                colnames(anc.lik) <- c(0:max(c(unique(sort(as.numeric(permutations[, 1]))), 1)))
                anc.lik[, as.character(unique(sort(as.numeric(permutations[, 1]))))] <- rep(1, length(permutations[, 1])-1)
            }
            for(j in 2:length(permutations[1, ])) { # Cycle through remaining permutations
                if(length(unique(sort(permutations[, j]))) > 1) {
                    if(length(mymodel[, 1]) > length(sort(unique(as.numeric(permutations[, j]))))) { # Case if model has larger spread than actual tipvalues (e.g. spans 0-5 when 2 is not recorded)
                        anc.lik <- ((anc.lik*(j-1))+ace(permutations[, j], chartree, type="discrete", model=mymodel[sort(unique(as.numeric(permutations[, j])))+1, sort(unique(as.numeric(permutations[, j])))+1])$lik.anc)/j # Get mean ancestral estimations
                    }
                    if(length(mymodel[, 1]) <= length(sort(unique(as.numeric(permutations[, j]))))) { # Case if model has same spread as tip values (theoretically this should always be true, but in reality not so much)
                        anc.lik <- ((anc.lik*(j-1))+ace(permutations[, j], chartree, type="discrete", model=mymodel)$lik.anc)/j # Get mean ancestral estimations
                    }
                }
                if(length(unique(sort(permutations[, j]))) == 1) {
                    anc.lik.perm <- matrix(0, nrow=length(permutations[, j])-1, ncol=max(c((unique(sort(as.numeric(permutations[, j])))+1), 2)))
                    colnames(anc.lik.perm) <- c(0:max(c(unique(sort(as.numeric(permutations[, j]))), 1)))
                    anc.lik.perm[, as.character(unique(sort(as.numeric(permutations[, j]))))] <- rep(1, length(permutations[, j])-1)
                    anc.lik <- ((anc.lik*(j-1))+anc.lik.perm)/j # Get mean ancestral estimations
                }
            }
        }
        # If no polymorphisms present just get ancestral likelihoods:
        if(length(grep("&", tipvals)) == 0 && length(unique(sort(as.numeric(tipvals)))) > 1) {
            if(length(mymodel[, 1]) > length(sort(unique(as.numeric(tipvals))))) { # Case if model has larger spread than actual tipvalues (e.g. spans 0-5 when 2 is not recorded)
                anc.lik <- ace(as.numeric(tipvals), chartree, type="discrete", model=mymodel[sort(unique(as.numeric(tipvals)))+1, sort(unique(as.numeric(tipvals)))+1])$lik.anc # Get ancestral state likelihoods for internal nodes
            }
            if(length(mymodel[, 1]) <= length(sort(unique(as.numeric(tipvals))))) { # Case if model has same spread as tip values (theoretically this should always be true, but in reality not so much)
                anc.lik <- ace(as.numeric(tipvals), chartree, type="discrete", model=mymodel)$lik.anc # Get ancestral state likelihoods for internal nodes
            }
        }
        if(length(grep("&", tipvals)) == 0 && length(unique(sort(as.numeric(tipvals)))) == 1) { # Case if all states are the same
            anc.lik <- matrix(0, nrow=length(tipvals)-1, ncol=max(c((unique(sort(as.numeric(tipvals)))+1), 2)))
            colnames(anc.lik) <- c(0:max(c(unique(sort(as.numeric(tipvals))), 1)))
            anc.lik[, as.character(unique(sort(as.numeric(tipvals))))] <- rep(1, length(tipvals)-1)
        }
        # Convert to just most likely states (i.e. 0 OR 1 for binary):
        anc.lik.vector <- vector(mode="numeric")
        max.lik <- apply(anc.lik, 1, max)
        for(j in 1:length(anc.lik[, 1])) {
            anc.lik.vector[j] <- as.numeric(colnames(anc.lik)[match(max.lik[j], anc.lik[j, ])]) # Case if single most likely state
            if(length(grep(max.lik[j], anc.lik[j, ])) > 1) anc.lik.vector[j] <- NA # Case if multiple equally most likely states (i.e. replace with NA)
        }
        chartreenodes <- (Ntip(chartree)+1):(Ntip(chartree)+Nnode(chartree)) # Nodes for character tree
        names(anc.lik.vector) <- chartreenodes
        # Copy information across to ancestral state matrix for whole tree:
        for(j in chartreenodes) {
            descs <- sort(chartree$tip.label[FindDescendants(j, chartree)]) # Find character tree descendants
            anc.lik.matrix[as.character(FindAncestor(descs, tree)), i] <- anc.lik.vector[as.character(j)]
        }
    }
    return(anc.lik.matrix)
}

get.bl <- function(state.matrix, tree, weights) {
    require(ape)
    tree.pd <- tree.nc <- tree.cc <- tree # Trees for patristic distance, number of changes and comparable characters
    for(i in 1:length(tree$edge[, 1])) {
        node.1 <- tree$edge[i, 1] # Get nodes for branch
        node.2 <- tree$edge[i, 2] # Get nodes for branch
        compchar <- sort(intersect(grep(TRUE, !is.na(state.matrix[node.1, ])), grep(TRUE, !is.na(state.matrix[node.2, ])))) # Characters that can be compared
        if(length(compchar) > 0) {
            polychar <- sort(unique(c(grep("&", state.matrix[node.1, compchar]), grep("&", state.matrix[node.2, compchar])))) # Characters with polymorphisms
            if(length(polychar) == 0) {
                tree.nc$edge.length[i] <- sum(sqrt(((as.numeric(state.matrix[node.1, compchar])-as.numeric(state.matrix[node.2, compchar]))*weights[compchar])^2))
                tree.pd$edge.length[i] <- sum(sqrt(((as.numeric(state.matrix[node.1, compchar])-as.numeric(state.matrix[node.2, compchar]))*weights[compchar])^2))/(length(compchar)*weights[compchar])
            }
            if(length(polychar) > 0) {
                diffs <- sqrt(( as.numeric(state.matrix[node.1, compchar])-as.numeric(state.matrix[node.2, compchar])) ^2)
                for (j in 1:length(polychar)) {
                    # Get minimum difference between polymorphic characters:
                    mindiff <- min(sqrt((as.numeric(strsplit(as.character(state.matrix[node.1, compchar[polychar[j]]]), "&")[[1]])-as.numeric(strsplit(as.character(state.matrix[node.2, compchar[polychar[j]]]), "&")[[1]]))^2))
                    diffs[polychar[j]] <- mindiff # Add to differences vector
                }
                tree.nc$edge.length[i] <- sum(diffs*weights[compchar])
                tree.pd$edge.length[i] <- sum(diffs*weights[compchar])/(length(compchar)*weights[compchar])
            }
            tree.cc$edge.length[i] <- length(compchar)*weights[compchar]
        }
        if(length(compchar) == 0) {
            tree.nc$edge.length[i] <- 0
            tree.pd$edge.length[i] <- 0
            tree.cc$edge.length[i] <- 0
        }
    }
    result <- list(tree.cc, tree.nc, tree.pd)
	names(result) <- c("tree.cc", "tree.nc", "tree.pd")
	return(result)
}

rbranch.phylo <- function(state.matrix, tree, ttree, cctree, permutations)
{
    charcorr <- (length(state.matrix[1, ])-apply(is.na(state.matrix), 1, sum))/length(state.matrix[1, ]) # Get proportion of completeness of nodes/tips
    ttree$edge.length <- (cctree$edge.length/length(state.matrix[1, ]))*ttree$edge.length # Normalise time tree branch lengths to reflect comparable characters number
    nchang <- sum(tree$edge.length) # Total number of character changes on tree
    permat <- matrix(0, nrow=length(tree$edge.length), ncol=permutations) # Matrix to store premutation results
    # Apportion branches to uniform distribution (between 0 and 1):
    brk <- vector(mode="numeric", length=length(ttree$edge.length))
	brk[1] <- ttree$edge.length[1]/length(ttree$edge.length)
    for (i in 2:length(brk)) brk[i] <- brk[i-1]+ttree$edge.length[i]/length(ttree$edge.length)
	brk <- brk/(sum(ttree$edge.length)/length(ttree$edge.length))
    # Assign changes to branches using a uniform distribution:
    for (i in 1:permutations) {
		randno <- runif(nchang, min=0, max=1)
		randno <- sort(randno)
		for (j in 1:length(tree$edge.length)) {
			breaker <- length(grep(TRUE, randno <= brk[j]))
			permat[j, i] <- breaker
			if (breaker > 0) randno <- randno[-(1:breaker)]
		}
	}
    sigs <- vector(mode="numeric") # Vector to store significant results
    for (i in 1:length(tree$edge.length)) ifelse(tree$edge.length[i] > sort(permat[i, ])[ceiling(0.95*length(permat[i, ]))], sigs[i] <- 1, sigs[i] <- 0) # Significance test for high rates
	return(sigs)
}

read.tnt <- function(file)
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
    if(length(grep("mxram", X)) > 0) { # If there is an nstates line
        X <- X[-grep("mxram", X)] # Remove nstates line
    }
    if(length(grep("nstates", X)) > 0) { # If there is an nstates line
        X <- X[-grep("nstates", X)] # Remove nstates line
    }
    if(length(grep("taxname", X)) > 0) { # If there is a taxname line
        X <- X[-grep("taxname", X)] # Remove taxname line
    }
    X <- X[-grep(TRUE, X == "xread")] # Remove xread
    if(length(grep(", ", X)) > 0) {
        header <- X[grep(", ", X)] # Find header text
        header <- gsub("'", "", header) # Remove quotes
        X <- X[-grep(", ", X)] # Remove header line(s)
    } else {
        header <- ""
    }
    nchar <- as.numeric(strsplit(X[1], " ")[[1]][1]) # Get number of characters
    ntax <- as.numeric(strsplit(X[1], " ")[[1]][2]) # Get number of taxa
    X <- X[-1] # Delete top line
    matrixblock <- X[1:ntax] # Get matrix block
    X <- X[-(1:ntax)] # Remove matrix block
    MATRIX <- matrix(nrow=ntax, ncol=nchar) # Make matrix
    rownames(MATRIX) <- c(1:ntax) # Dummy rownames
    for(i in 1:ntax) { # For each taxon
        rownames(MATRIX)[i] <- strsplit(matrixblock[i], " ")[[1]][1] # Extract taxon name
        matrixblock[i] <- strsplit(matrixblock[i], " ")[[1]][length(strsplit(matrixblock[i], " ")[[1]])] # And character block
        chars <- strsplit(matrixblock[i], "")[[1]] # Get character vector
        j <- 1 # Start value for j
        while(j <= length(chars)) { # Go through character vector until all characters have been looked at
            if(chars[j] != "[") { # If the character is not a polymorphism
                MATRIX[i, grep(TRUE, is.na(MATRIX[i, ]))[1]] <- chars[j] # Just store it in the matrix
                j <- j+1 # And move to the next one
			} else { # If it is a polymorphism
                MATRIX[i, grep(TRUE, is.na(MATRIX[i, ]))[1]] <- paste(chars[(j+1):(j+grep(TRUE, chars[j:length(chars)] == "]")[1]-2)], collapse="&") # Paste all states into the matrix
                j <- j+grep(TRUE, chars[j:length(chars)] == "]")[1] # And move to the next character
			}
		}
	}
	MATRIX <- gsub("\\?", NA, MATRIX) # Replace question marks with NAs
	unq.states <- sort(unique(as.vector(MATRIX))) # Find unique states list
	if(length(grep("&", unq.states)) > 0) { # If polymorphisms are present
		polys <- unq.states[grep("&", unq.states)] # Isolate them
		unq.states <- unq.states[-grep("&", unq.states)] # Remove them
		for(i in 1:length(polys)) { # Go through each polymorphism
			unq.states <- c(unq.states, strsplit(polys[i], "&")[[1]]) # Split adn add to unqiue states list
		}
		unq.states <- sort(unique(unq.states)) # Get unique states again
	}
	if(length(sort(match(LETTERS, unq.states))) > 0) { # If there are letters (greater than ten states)
		lets.found <- unq.states[sort(match(LETTERS, unq.states))] # List leters found
		for(i in 1:length(lets.found)) { # For each leter found
			MATRIX <- gsub(lets.found[i], match(lets.found[i], LETTERS)+10, MATRIX) # Replace with appropriate number
		}
	}
	if(length(grep("ccode", X)) > 0) { # If weights and/or orderings are present
		ordering <- weights <- vector(mode="character", length=nchar)
		wt.ord <- gsub(";", "", gsub("ccode ", "", X[grep("ccode", X)])) # Get character weights/orderings
		wt.ord <- strsplit(wt.ord, " ")[[1]]
		char.nos <- as.numeric(wt.ord[c(1:length(wt.ord))[grep(TRUE, c(1:length(wt.ord)) %% 2 == 0)]])+1 # Get character numbers (even cells)
		char.vals <- wt.ord[c(1:length(wt.ord))[grep(TRUE, c(1:length(wt.ord)) %% 2 != 0)]] # Get character values (odd cells)
		ordering[grep("-", char.vals)] <- "unord" # Get unordered characters
		ordering[grep("\\+", char.vals)] <- "ord" # Get ordered characters
		for(i in 1:nchar) { # For each character
			if(length(grep("\\[", char.vals[i])) > 0) { # If character is active
                weights[i] <- as.numeric(strsplit(char.vals[i], "/")[[1]][2]) # Set weight
			} else { # If character is inactive
                weights[i] <- 0 # Set weight to zero
			}
		}
	} else { # If they are not
		weights <- rep(1, nchar) # Make all characters weighted one
		ordering <- rep("unord", nchar) # Make all characters unordered
	}
	max.vals <- min.vals <- vector(mode="numeric", length=nchar) # Min and max value vectors
	for(i in 1:nchar) { # For each character
		unq.states <- sort(unique(as.vector(MATRIX[, i])))
		if(length(grep("&", unq.states)) > 0) { # If polymorphisms are present
			polys <- unq.states[grep("&", unq.states)] # Isolate them
			unq.states <- unq.states[-grep("&", unq.states)] # Remove them
			for(j in 1:length(polys)) { # Go through each polymorphism
				unq.states <- c(unq.states, strsplit(polys[j], "&")[[1]]) # Split and add to unqiue states list
			}
			unq.states <- sort(unique(unq.states)) # Get unique states again
		}
		min.vals[i] <- min(as.numeric(unq.states))
		max.vals[i] <- max(as.numeric(unq.states))
	}
	result <- list(header, MATRIX, ordering, weights, max.vals, min.vals)
	names(result) <- c("header", "matrix", "ordering", "weights", "max.vals", "min.vals")
	return(result)
}

read.paleodb.occs <- function(file, sep=", ", header=TRUE)
{
	# Not sure why I have to do this, but it seems to break otherwise:
	Sys.setlocale('LC_ALL', 'C')
	
	# Read in raw data:
	X <- scan(file=file, what="", sep="\n", quiet=TRUE)
	
	# Whilst there are gaps (empty cells) in the data:
	while(length(grep(paste(sep, sep, sep=""), X)) > 0) {
	
		# Insert NAs into gaps (excluding first and last columns):
		X <- gsub(paste(sep, sep, sep=""), paste(sep, "NA", sep, sep=""), X)
	
	}
	
	# To check ends of lines for gaps go line by line:
	for(i in 1:length(X)) {
		
		# If first value is empty replace with NA:
		if(strsplit(X[i], "")[[1]][1] == sep) X[i] <- paste(NA, X[i], sep="")
		
		# If last value is empty replace with NA:
		if(strsplit(X[i], "")[[1]][length(strsplit(X[i], "")[[1]])] == sep) X[i] <- paste(X[i], NA, sep="")

	}
	
	# Catch separators in text (i.e. when followed by a space):
	X <- gsub(paste(sep, " ", sep=""), paste("%sEpArAtOr%", " ", sep=""), X)

	# For each line of the text file:
	for(i in 1:length(X)) {
		
		# Break text at quotation marks:
		new.block <- strsplit(X[i], "\"")[[1]]
		
		# Find within quotes parts (i.e. even numbers):
		evens <- grep(TRUE, c(1:length(new.block)) %% 2 == 0)
		
		# Replace seperation values present within quotes:
		new.block[evens] <- gsub(sep, "%sEpArAtOr%", new.block[evens])
		
		# Reinsert quotes and store in X:
		X[i] <- paste(new.block, collapse="\"")
		
	}
	
	# Create vector to record the number of separator values on each line:
	seps <- vector(mode="numeric")
	
	# For each line:
	for(i in 1:length(X)) {
		
		# Count the number of separator values:
		seps[i] <- length(grep(TRUE, strsplit(X[i], "")[[1]] == sep))
		
	}
	
	# If there are still separator values inside fields:
	if(length(unique(seps)) > 0) {
		
		# Identify problem rows:
		problem.rows <- grep(TRUE, seps > min(seps))
		
		# For each problem row:
		for(i in problem.rows) {
			
			# Isolate values by separating them:
			row.values <- strsplit(X[i], sep)[[1]]
			
			# Set vectors for storing values where first or last characters are quotes:
			first.character.a.quote <- last.character.a.quote <- vector(mode="numeric")
			
			# For each value:
			for(j in 1:length(row.values)) {
				
				# If starts with a quote then store:
				if(strsplit(row.values[j], "")[[1]][1] == "\"") first.character.a.quote <- c(first.character.a.quote, j)
				
				# If ends with a quote then store:
				if(strsplit(row.values[j], "")[[1]][length(strsplit(row.values[j], "")[[1]])] == "\"") last.character.a.quote <- c(last.character.a.quote, j)
				
			}
			
			# Fields that begin and end with quotes are text and have already been dealt with above:
			text.fields <- intersect(first.character.a.quote, last.character.a.quote)
			
			# Other quotes are stray and may indicate problematic fields:
			stray.quotes <- sort(c(setdiff(first.character.a.quote, last.character.a.quote), setdiff(last.character.a.quote, first.character.a.quote)))
			
			# Identify fields that contain text (but may not be bookended by quotes):
			contains.text <- grep("[A-Z:a-z]", row.values)
			
			# Identify initial list of fields that are OK:
			OK.fields <- sort(c(text.fields, grep(TRUE, row.values == "NA")))
			
			# Identify initial list of fields that are potentially not OK:
			notOK.fields <- setdiff(c(1:length(row.values)), OK.fields)
			
			# Update Ok fields to exclude number fields:
			OK.fields <- sort(unique(c(OK.fields, setdiff(notOK.fields, contains.text))))

			# Update not OK fields:
			notOK.fields <- setdiff(c(1:length(row.values)), OK.fields)

			# If first value is not OK, but second value is OK:
			if(length(grep(TRUE, notOK.fields == 1)) > 0 && length(grep(TRUE, OK.fields == 2))) {
				
				# First value is now OK:
				OK.fields <- c(1, OK.fields) 
				
				# First value is no longer not OK:
				notOK.fields <- notOK.fields[-1]
				
			}
			
			# If first value is not OK, but second value is OK:
			if(length(grep(TRUE, notOK.fields == length(row.values))) > 0 && length(grep(TRUE, OK.fields == (length(row.values) - 1)))) {
				
				# First value is now OK:
				OK.fields <- c(OK.fields, length(row.values))
				
				# First value is no longer not OK:
				notOK.fields <- notOK.fields[-grep(TRUE, notOK.fields == length(row.values))]
				
			}
			
			# For each remaining not Ok field:
			for(j in notOK.fields) {
				
				# Check that it is not an end value:
				if(j != 1 && j != length(row.values)) {
					
					# If values either side are OK:
					if(length(grep(TRUE, OK.fields == (j + 1))) && length(grep(TRUE, OK.fields == (j - 1)))) {
						
						# First value is now OK:
						OK.fields <- sort(c(OK.fields, j))
						
						# First value is no longer not OK:
						notOK.fields <- notOK.fields[-grep(TRUE, notOK.fields == j)]
						
					}
					
				}
				
			}
			
			
			while(length(notOK.fields) > 1) {
				
				
				
			}
			
			
			
		
		}
		
	}
	
	
	
	
	
	if(header == TRUE) {
		
		out <- matrix(ncol=min(seps) + 1, nrow=length(X) - 1)
		
		colnames(out) <- strsplit(X[1], sep)[[1]]
		
		X <- X[-1]
		
	}
	
	
	if(header == FALSE) out <- matrix(ncol=min(seps) + 1, nrow=length(X))
	
	
	
	for(i in 1:length(X)) out[i, ] <- strsplit(X[i], sep)[[1]]
	
	
	
	
	# Reinsert non-separating separator characters:
	out <- gsub("%sEpArAtOr%", sep, out)
	
	# Remove now redundant quotes:
	out <- gsub("\"", "", out)
	
	# Return output as table:
	return(out)
	
}

# Function to find nodes that can be realistically dated by the Hedman technique:
find.dateable.nodes <- function(tree, tip.ages) {

	# Requires ape library for reading tree structure:
	require(ape)
	
	# List all internal nodes except root:
	list.nodes <- (Ntip(tree) + 2):(Ntip(tree) + Nnode(tree))
	
	# Create vector to store nodes that are dateable:
	defo.nodes <- vector(mode="numeric")
	
	# For each node in the list:
	for(i in length(list.nodes):1) {
		
		# Find descendant nodes:
		descs <- tree$edge[grep(TRUE, tree$edge[, 1] == list.nodes[i]), 2]
		
		# If all immediate descendants are internal nodes than it cannot be dated and is removed:
		if(length(which(descs > Ntip(tree))) == length(descs)) list.nodes <- list.nodes[-i]
		
		# If all immediate descendants are terminal nodes than it can definitely be dated...:
		if(length(which(descs <= Ntip(tree))) == length(descs)) {
			
			# ...and is retained...:
			defo.nodes <- c(defo.nodes, list.nodes[i])
			
			# ...but can be removed from list.nodes:
			list.nodes <- list.nodes[-i]
		
		}
	
	}
	
	# For each internal node with both a terminal and internal node descendant:
	for(i in length(list.nodes):1) {
		
		# Find node age:
		node.age <- max(tip.ages[FindDescendants(list.nodes[i], tree)])
		
		# List its descendants:
		descs.ages <- descs <- tree$edge[grep(TRUE, tree$edge[, 1] == list.nodes[i]), 2]
		
		# For each descendant:
		for(j in length(descs.ages):1) {
			
			# If an internal node date as oldest descendant:
			if(descs[j] > Ntip(tree)) descs.ages[j] <- max(tip.ages[FindDescendants(descs[j], tree)])
			
			# Remove if terminal node:
			if(descs[j] <= Ntip(tree)) descs.ages <- descs.ages[-j]
			
		}
		
		# If node age is older than any descendant internal nodes:
		if(node.age > max(descs.ages)) {
			
			# Then can be dated so add to list:
			defo.nodes <- c(defo.nodes, list.nodes[i])
			
		}
		
	}
	
	# Output all dateable nodes:
	return(defo.nodes)

}

# Function that retrieves dates for use as tnodes variable in Hedman function
get.tnodes <- function(node, tree, tip.ages) {
	
	# Require ape library:
	require(ape)
	
	# Number of root node (point at which searching for ancestral nodes stops):
	root.node <- Ntip(tree) + 1
	
	# If the input node is the root node:
	if(node == root.node) {
		
		# Error and warning:
		stop("ERROR: Input node is root node - tnodes has length one!")
		
	# If input node is not root node:
	} else {
		
		# Find first ancestor node:
		internodes <- tree$edge[match(node, tree$edge[, 2]), 1]
		
		# As long as we have not yet reached the root:
		while(internodes[length(internodes)] != root.node) {
			
			# Add next internode to set:
			internodes <- c(internodes, tree$edge[match(internodes[length(internodes)], tree$edge[, 2]), 1])
		}
		
		# Add in original node as that is counted too!:
		internodes <- c(node, internodes)
		
		# Create tnodes vector for storage:
		tnodes <- vector(mode="numeric")
		
		# Date first (original) node:
		tnodes[1] <- max(tip.ages[FindDescendants(internodes[1], tree)])
		
		# For each internode:
		for(i in 2:length(internodes)) {
			
			# Find descendants:
			descs <- tree$edge[grep(TRUE, tree$edge[, 1] == internodes[i]), 2]
			
			# Remove previously dated node:
			descs <- descs[-match(internodes[(i - 1)], descs)]
			
			# For each other descendant:
			for(j in 1:length(descs)) {
				
				# If descendant node is internal:
				if(descs[j] > Ntip(tree)) {
					
					# Get maximum possible age:
					descs[j] <- max(tip.ages[FindDescendants(descs[j], tree)])
					
				# If descendant node is terminal:
				} else {
					
					# Get maximum possible age:
					descs[j] <- tip.ages[descs[j]]
				
				}
			
			}
			
			# Store in tnodes:
			tnodes[i] <- max(descs)
		
		}
		
		# Vector for storing node ages to be deleted:
		deletes <- vector(mode="numeric")
		
		# For each potential tnode:
		for(i in 2:length(tnodes)) {
			
			# Find nodes too young for dating and store in deletes:
			if(max(tnodes[1:(i-1)]) > tnodes[i]) deletes <- c(deletes, i)
			
		}
		
		# Remove these nodes (if there are any) from tnodes:
		if(length(deletes) > 0) tnodes <- tnodes[-deletes]
		
		# Give in order from oldest to youngest:
		tnodes <- rev(tnodes)
		
		# Return tnodes:
		return(tnodes)
		
	}
	
}

# Hedman (2010) method for estimating confidence given ages of ougroups (tnodes is the sequence of ages, from oldest to youngest; t0 is the arbitrary lower stratigraphic bound; resolution is the number of steps to take between the FAD and the lower stratigraphic bound):
Hedman.2010 <- function (tnodes, t0, resolution) { # Outputs estimate with two-tailed 95% CIs
    
    # Check t0 is older than any other node age:
    if(!all(t0 > tnodes)) stop("t0 must be older than any tnode value.")
    
	# Store requested resolution:
	requested.resolution <- resolution
	
	# Function returns p.d.f. for node ages	(t0 is an arbitrary oldest age to consider, tnodes are the ages of oldest known fossil stemming from each node, and tsteps is a vector of arbitrary time steps on which the p.d.f. is calculated):
	nodeage <- function(tsteps, tnodes, t0) {

		# Get number of outgroups (tnodes):
		nn  <- length(tnodes)
		
		# Get number of time steps (resolution):
		nt  <- length(tsteps)
		
		# Initialize array:
		pnodes  <- matrix(0, nn, nt)
		
		# First get pdf for node 1, at discrete values:
		ii  <- which(tsteps > t0 & tsteps < tnodes[1])
		
		# Assume uniform distribution for oldest node:
		pnodes[1, ii] <- 1.0 / length(ii)
		
		# Cycle through remaining nodes:
		for (i in 2:nn) {
			
			# Cycle through series of time steps:
			for (j in 1:nt) {
				
				# Initialize vector:
				p21  <- rep(0, nt)
				
				# Get p.d.f. for ith node ay discrete values:
				ii <- which(tsteps >= tsteps[j] & tsteps < tnodes[i])
				
				
				p21[ii]  <- 1.0 / length(ii) * pnodes[(i - 1), j]
				
				# Conditional probability of this age, given previous node age, times probability of previous node age, added to cumulative sum:
				pnodes[i, ii]  <- pnodes[i, ii] + p21[ii]
				
			}
			
		}
		
		# Get just the p.d.f. vector for the youngest node:
		out <- pnodes[nn, ]
		
		# Return output:
		return(out)
		
	}
	
	# Calculate c.d.f. from p.d.f., find median and 95% credibility limits:
	HedmanAges <- function(tnodes, t0, resolution) {
	
		# Store input resolution for later reuse:
		old.resolution <- resolution
	
		# Get uniformly spaced p-values (including 95% limits adn median) for calculating CIs:
		CIs <- sort(unique(c(0.025, 0.975, 0.5, c(1:old.resolution) * (1 / old.resolution))))
	
		# Get enw resolution (may be longer by adding limits and/or median):
		resolution <- length(CIs)
	
		# Get oldest possible age to consider:
		first <- t0
	
		# Get youngest possible age to consider:
		last <- min(tnodes)
	
		# Make tnodes negative for Hedman function:
		tnodes  <- -tnodes
	
		# Get uniformly spaced time steps for Hedman function:
		tsteps <- seq(-first, -last, length=resolution)

		# Make t0 negative for Hedman function:
		t0 <- -t0
	
		# Run Hedman function to get p.d.f.:
		vector <- nodeage(tsteps, tnodes, t0)
	
		# Convert from p-value scale to age scale:
		integral.vector <- vector * ((abs(first - last)) / resolution)
	
		# Get sum of probabilities in order to re-scale as p.d.f.:
		probability.sum <- sum(integral.vector)
	
		# Get re-scaled p-values for CIs:
		ps.CIs <- probability.sum * CIs
	
		# Cretae empty vector to store dates:
		date.distribution <- vector(mode="numeric")
	
		# Set initial value:
		value <- 0
	
		# For each re-scaled p-value:
		for(i in length(ps.CIs):1) {
		
			# Update value:
			value <- value + integral.vector[i]
		
			# If in age window:
			if(length(which(value >= ps.CIs)) > 0) {

				# Store dates:
				date.distribution <- c(date.distribution, rep(tsteps[i], length(which(value >= ps.CIs))))
			
				# Update re-scaled p-values:
				ps.CIs <- ps.CIs[-grep(TRUE, value >= ps.CIs)]
			
			}
		
		}
	
		# Add t0 at end if length is short:
		while(length(date.distribution) < resolution) date.distribution <- c(date.distribution, t0)

		# Find median of distribution:
		Best.guess <- date.distribution[CIs == 0.5]
			
		# Get upper 95% CI:
		Best.guess.lower <- date.distribution[CIs == 0.975]
		
		# Get lower 95% CI:
		Best.guess.upper <- date.distribution[CIs == 0.025]
		
		# Combine results and reverse sign to give palaeo ages:
		results <- list(-Best.guess, -Best.guess.lower, -Best.guess.upper, -date.distribution[match(c(c(1:old.resolution) * (1 / old.resolution)), CIs)])
		
		# Name subvariables:
		names(results) <- c("Best.guess", "Best.guess.lower", "Best.guess.upper", "Age.distribution")
		
		# Return result
		return(results)
		
	}
	
	# Get Hedman ages:
	out <- HedmanAges(tnodes, t0, resolution)
	
	# If Hedman ages represent less than five unique values (leading to downstream problem of flat distributions):
	while(length(unique(out$Age.distribution[round(seq(1, length(out$Age.distribution), length.out=requested.resolution))])) < 5) {
		
		# Double the resolution size:
		resolution <- resolution * 2
		
		# Get Hedman ages:
		out <- HedmanAges(tnodes, t0, resolution)
		
	}
	
	# Update output:
	out$Age.distribution <- out$Age.distribution[round(seq(1, length(out$Age.distribution), length.out=requested.resolution))]
	
	# Return output:
	return(out)
	
}

# Over-arching Hedman tree-dating function:
Hedman.tree.dates <- function(tree, tip.ages, outgroup.ages, t0, resolution = 1000, conservative = TRUE) {
	
# MORE CONDITIONALS NEEDED FOR WHEN RESOLUTION IS LOW

	# Load libraries:
	require(ape)
	require(strap)
    
    # Check tree is fully bifurcating:
    if(!is.binary.tree(tree)) stop("Tree must be fully bifurcating.")
    
	# Ensure tip ages are in tree tip order:
	tip.ages <- tip.ages[tree$tip.label]

	# Find root node:
	root.node <- Ntip(tree) + 1
	
	# Create variables to store age estimates:
	age.estimates <- matrix(nrow = Nnode(tree), ncol = 3)
	
	# Create variables to store age distributions:
	age.distributions <- matrix(nrow = Nnode(tree), ncol = resolution)
	
	# Set column headings:
	colnames(age.estimates) <- c("Best.guess", "Best.guess.lower", "Best.guess.upper")

	# Set row names:
	rownames(age.estimates) <- rownames(age.distributions) <- c((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
	
	# Vector to store outgroup sequences:
	og.seq <- vector(mode = "numeric")

	# Report progress:
	cat("Identifying nodes that are date-able using the Hedman technique_")
	
	# Find dateable conservative nodes:
	if(conservative) dateable.nodes <- find.dateable.nodes(tree, tip.ages)
	
	# Find dateable non-conservative nodes:
	if(!conservative) dateable.nodes <- c((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))

	# Report progress:
	cat("Done\nDating nodes using Hedman technique_")

	# If not using the conservative approach:
	if(!conservative) {

		# Create ages matrix:
		ages <- cbind(tip.ages, tip.ages)

		# Add rownames:
		rownames(ages) <- names(tip.ages)

		# Add column names:
		colnames(ages) <- c("FAD", "LAD")

		# Sort by taxon order in tree:
		ages <- ages[tree$tip.label, ]

		# Date tree using basic algorithm:
		stree <- DatePhylo(tree, ages, 0, "basic", FALSE)

		# Get node ages:
		sages <- GetNodeAges(stree)

		# Set first outgroup sequence (i.e. root) as this is special case:
		og.seq[1] <- paste(c(outgroup.ages, sages[root.node]), collapse = "%%")

		# For each non-root node:
		for(i in 2:length(dateable.nodes)) {
			
			# Identify node:
			node <- dateable.nodes[i]
			
			# Find preceding node:
			internodes <- tree$edge[match(node, tree$edge[, 2]), 1] 
			
			# Keep going until we reach the root and add next internode to set:
			while(internodes[length(internodes)] != root.node) internodes <- c(internodes, tree$edge[match(internodes[length(internodes)], tree$edge[, 2]), 1])
			
			# Collate nodes and put in order:
			internodes <- sort(c(node, internodes))
			
			# Find outgroup age sequence and collaspe to single string and store:
			og.seq[i] <- paste(c(outgroup.ages, sages[internodes]), collapse="%%")

		}

	# If using the conservative approach:
	} else {

		# Find age estimates for each dateable node:
		for(i in 1:length(dateable.nodes)) {
		
			# Find tnodes for a specific node:
			tnodes <- get.tnodes(dateable.nodes[i], tree, tip.ages)
		
			# Add additional outgroup tnodes:
			tnodes <- c(outgroup.ages, tnodes)

			# Turn outgroup sequence into single string and store:
			og.seq[i] <- paste(tnodes, collapse="%%")
			
		}

	}
	
	# Get unique outgroup sequences:
	unq.og.seq <- unique(og.seq)

	# For each unique outgroup sequence:
	for(i in 1:length(unq.og.seq)) {
		
		# List nodes with this outgroup sequence:
		nodes <- dateable.nodes[which(og.seq == unq.og.seq[i])]
		
		# Define tnodes:
		tnodes <- as.numeric(strsplit(unq.og.seq[i], "%%")[[1]])

		# Get Hedman dates:
		hedman.out <- Hedman.2010(tnodes = tnodes, t0 = t0, resolution = resolution)
		
		# Update age distributions:
        for(j in 1:length(nodes)) age.distributions[as.character(nodes[j]), ] <- hedman.out$Age.distribution
		
	}
	
	# Separate check to see if root is dateable (first find root descendants that are tips, if any):
	root.desc.tips <- sort(tree$tip.label[tree$edge[which(root.node == tree$edge[, 1]), 2]])
	
	# Now check that there are descendants that are tips:
	if(length(root.desc.tips) > 0) {
		
		# Get maximum root descendant age:
		root.tip.age <- max(tip.ages[root.desc.tips])

		# Now check that descendants are older than any other tips and get Hedman dates for root and update age distributions:
		if(root.tip.age > max(tip.ages[setdiff(tree$tip.label, root.desc.tips)])) age.distributions[as.character(root.node), ] <- Hedman.2010(tnodes = c(outgroup.ages, root.tip.age), t0 = t0, resolution = resolution)$Age.distribution
		
	}

	# If not using the conservative approach:
	if(conservative == FALSE) {
	
		# Report progress:
		cat("Done\nIdentifying remaining undated nodes_")
		
		# Report progress:
		cat("Done\nDating remaining undated nodes using randomisation technique_")

	# If using the conservative approach:
	} else {

		# Report progress:
		cat("Done\nIdentifying remaining undated nodes_")
		
		# Find undated nodes by first establishing actually dated nodes:
		dated.nodes <- sort(dateable.nodes)
		
		# List all internal nodes:
		all.nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
		
		# Get undated nodes:
		undated.nodes <- setdiff(all.nodes, dated.nodes)
		
		# Find branches that define undated nodes:
		internal.branches <- tree$edge[which(tree$edge[, 2] > Ntip(tree)), ]
		
		# Is root dated? (default to "yes" before actually checking):
		root.dated <- "Yes"
		
		# If the root is undated (then it needs to be):
		if(is.na(match(root.node, dated.nodes))) {
			
			# Update root.dated:
			root.dated <- "no"
			
			# Create tnodes from outgroup ages ONLY (this will later be used to constrain the root age through the saem randomisation process as the other undated nodes):
			tnodes <- outgroup.ages
			
			# Date root node:
			hedman.out <- Hedman.2010(tnodes = tnodes, t0 = t0, resolution = resolution)
			
			# Add to age estimates:
			age.estimates <- rbind(age.estimates, c(hedman.out$Best.guess, hedman.out$Best.guess.lower, hedman.out$Best.guess.upper))
			
			# Add to age distributions:
			age.distributions <- rbind(age.distributions, hedman.out$Age.distribution)
			
			# Add node name 0:
			rownames(age.estimates)[length(rownames(age.estimates))] <- rownames(age.distributions)[length(rownames(age.distributions))] <- "0"
			
		}
		
		# Vector for storing strings of sequences of undated node(s) bound by dated nodes:
		anc.strings.trimmed <- anc.strings <- vector(mode = "character")
		
		# Work through dated nodes to find sequences of undated nodes bracketed by them:
		for(i in 1:length(dated.nodes)) {
			
			# As long as the dated node is not the root (which has no ancestors):
			if(dated.nodes[i] != root.node) {
				
				# Get ancestors:
				ancestors <- internal.branches[which(internal.branches[, 2] == dated.nodes[i]), 1]
				
				# As long as the ancestor is neither the root or a dated node:
				if(ancestors != root.node && length(sort(match(ancestors, dated.nodes))) == 0) {
					
					# And as long as it remains so find next ancestral node:
					while(ancestors[length(ancestors)] != root.node && length(sort(match(ancestors[length(ancestors)], dated.nodes))) == 0) ancestors <- c(ancestors, internal.branches[match(ancestors[length(ancestors)], internal.branches[, 2]), 1])
					
				}
				
				# Add initial dated node:
				ancestors <- c(dated.nodes[i], ancestors)
				
				# Only need retain those with at least one undated node between dated nodes:
				if(length(ancestors) > 2) {
					
					# Add to ancestor strings:
					anc.strings <- c(anc.strings, paste(ancestors, collapse = " "))
					
					# Add to trimmed ancestor strings:
					anc.strings.trimmed <- c(anc.strings.trimmed, paste(ancestors[2:(length(ancestors) - 1)], collapse = " "))
					
				}
				
			}
			
		}
		
		# If root is undated then anc.strings with the root present must be modified to include node 0 at the end:
		if(root.dated == "no") {
			
			# Go through each anc.string:
			for(i in 1:length(anc.strings)) {
				
				# Get ancestors vector:
				ancestors <- as.numeric(strsplit(anc.strings[i], " ")[[1]])
				
				# If last node is the root:
				if(ancestors[length(ancestors)] == root.node) {
					
					# Add 0 node to end:
					ancestors <- c(ancestors, 0)
					
					# Re-store in anc.strings:
					anc.strings[i] <- paste(ancestors, collapse = " ")
					
				}
				
			}
			
		}
		
		# Find oldest node from each ancestral sequence (first set up empty vector):
		oldest.nodes <- vector(mode="numeric")
		
		# Go through each anc.string:
		for(i in 1:length(anc.strings)) {
			
			# Get ancestors vector:
			ancestors <- as.numeric(strsplit(anc.strings[i], " ")[[1]])
			
			# Store oldest node:
			oldest.nodes <- c(oldest.nodes, ancestors[length(ancestors)])
			
		}
		
		# Collapse to unique oldest nodes only:
		oldest.nodes <- sort(unique(oldest.nodes))
		
		# Clump anc.strings by shared oldest node (i.e. blocks that will be dated at the same time) - first create empty vector:
		clumped.anc.strings <- vector(mode="character")
		
		# For each oldest node find anc.strings with oldest node and collapse with %%:
		for(i in 1:length(oldest.nodes)) clumped.anc.strings <- c(clumped.anc.strings, paste(anc.strings[grep(paste(" ", oldest.nodes[i], sep = ""), anc.strings)], collapse = "%%"))

		# Report progress:
		cat("Done\nDating remaining undated nodes using randomisation technique_")
		
		# Can now date undated nodes using random draws from bounding Hedman dated nodes:
		for(i in 1:length(clumped.anc.strings)) {
			
			# First retrieve anc.strings with shared oldest node:
			anc.strings <- strsplit(clumped.anc.strings[i], "%%")[[1]]
			
			# Get young dated nodes (for drawing from for upper bounds):
			young.dated.nodes <- vector(mode="numeric")
			
			# For each ancestor string find youngest dated node, i.e., lower bounds::
			for(j in 1:length(anc.strings)) young.dated.nodes <- c(young.dated.nodes, as.numeric(strsplit(anc.strings[j], " ")[[1]][1]))
			
			# Get age distributions for young dated nodes:
			young.distributions <- age.distributions[as.character(young.dated.nodes), ]
			
			# Force into a single row matrix if only a vector:
			if(!is.matrix(young.distributions)) young.distributions <- t(as.matrix(young.distributions))
			
			# Set young distribution half up from young distributions:
			young.distributions.old.half <- young.distributions
			
			# Get old half as age distributions for young dated nodes:
			for(j in 1:length(apply(young.distributions, 1, median))) {
				
				# As long as the median is not the maximum:
				if(max(young.distributions[j, ]) > median(young.distributions[j, ])) {
					
					# Split distribution into older half using median:
					temp.distribution <- young.distributions[j, which(young.distributions[j, ] > median(young.distributions[j, ]))]
					
				# If the median is the maximum:
				} else {
					
					# Just use the maximum:
					temp.distribution <- max(young.distributions[j, ])
					
				}
				
				# ???:
				young.distributions.old.half[j, ] <- c(temp.distribution, rep(NA, length(young.distributions[j, ]) - length(temp.distribution)))
			
			}
			
			# Get old dated node (for drawing from for lower bounds):
			old.dated.node <- as.numeric(strsplit(anc.strings[1], " ")[[1]][length(strsplit(anc.strings[1], " ")[[1]])])
			
			# Get age distributions for old dated node:
			old.distribution <- age.distributions[as.character(old.dated.node), ]
			
			# As long as the median is not the minimum:
			if(min(old.distribution) < median(old.distribution)) {
				
				# Split distribution using median:
				old.distribution.young.half <- old.distribution[which(old.distribution < median(old.distribution))]
				
			# If the median is the minimum:
			} else {
				
				# Just use the minimum:
				old.distribution.young.half <- min(old.distribution)
			
			}
			
			# Get undated nodes (vector for storing output):
			undated.nodes <- vector(mode="numeric")
			
			# For each anc.string:
			for(j in 1:length(anc.strings)) {
				
				# Get ancestors:
				ancestors <- strsplit(anc.strings[j], " ")[[1]]
				
				# Store undated nodes:
				undated.nodes <- c(undated.nodes, as.numeric(ancestors[2:(length(ancestors) - 1)]))
			
			}
			
			# Collapse to unique nodes only:
			undated.nodes <- sort(unique(undated.nodes))
			
			# Find young constraints for each undated node (vector for storing results):
			young.constraints <- vector(mode = "numeric")
			
			# For each undated node add young constraints to list:
			for(j in 1:length(undated.nodes)) young.constraints <- c(young.constraints, paste(young.dated.nodes[grep(paste(" ", undated.nodes[j], sep = ""), anc.strings)], collapse = " "))
			
			# Find old constraints for each undated node (i.e. all nodes that are lower in tree as these will potentially be dated first and replace the current oldest age; vector for storing results):
			old.constraints <- vector(mode="numeric")
			
			# For each undated node:
			for(j in 1:length(undated.nodes)) {
				
				# Get just anc.strings where undated node is present:
				node.strings <- anc.strings[grep(paste(" ", undated.nodes[j], sep=""), anc.strings)]
				
				# Vector for storing older nodes:
				older.nodes <- vector(mode="numeric")
				
				# For each anc.string where the undated node exists ?????:
				for(k in 1:length(node.strings)) older.nodes <- c(older.nodes, as.numeric(strsplit(strsplit(node.strings[k], undated.nodes[j])[[1]][2], " ")[[1]]))
				
				# ????:
				old.constraints <- c(old.constraints, paste(unique(sort(older.nodes)), collapse=" "))
			
			}
			
			# Find undated nodes constrained by dated nodes (vector for storage):
			young.dated.nodes.constrain <- vector(mode = "character")
			
			# For each young dated node constraint find undated nodes which it contrains:
			for(j in 1:length(young.dated.nodes)) young.dated.nodes.constrain[j] <- paste(as.numeric(strsplit(anc.strings[j], " ")[[1]][2:(length(strsplit(anc.strings[j], " ")[[1]]) - 1)]), collapse = " ")
			
			# Main dating loop (repeats random draw N times, where N is defined by the variable resolution):
			for(j in 1:resolution) {
				
				# Modify age distributions used if medians of undated nodes exceed bounds of medians of dated nodes (needs to have been at least two entries already):
				if(j >= 3 && j >= floor(resolution / 2)) {
					
					# If there is more than one undated node establish present medians for undated nodes:
					if(length(undated.nodes) > 1) undated.medians <- apply(age.distributions[as.character(undated.nodes), 1:(j - 1)], 1, median)
					
					# If there is only one undated node
					if(length(undated.nodes) == 1) {
						
						# Establish present median of undated node:
						undated.medians <- median(age.distributions[as.character(undated.nodes), 1:(j - 1)])
						
						# Add name for reference later:
						names(undated.medians) <- as.character(undated.nodes)
					
					}
					
					# Case if oldest median of the undated nodes exceeds the median of the bounding old dated node:
					if(max(undated.medians) >= median(sort(old.distribution))) {
						
						# Update active distribution with older half only (to force undated nodes towards median ages consistent with dated nodes):
						active.old.distribution <- old.distribution.young.half
						
					# If oldest median ages do not conflict:
					} else {
						
						# Draw from full distribution:
						active.old.distribution <- old.distribution
						
					}
					
					# Set default active young distributions:
					active.young.distributions <- young.distributions
					
					# For each young (top bounding) dated node:
					for(k in 1:length(young.dated.nodes)) {
						
						# Case if youngest median of the undated nodes exceeds the median of the bounding young dated node - update active distribution with younger half only (to force undated nodes towards median ages consistent with dated nodes):
						if(min(undated.medians[strsplit(young.dated.nodes.constrain[k], " ")[[1]]]) <= median(sort(age.distributions[as.character(young.dated.nodes[k]), ]))) active.young.distributions[k, ] <- young.distributions.old.half[k, ]
					
					}
					
				# Case if still in first half of resolution:
				} else {
					
					# Use uncorrected distribution for upper bound:
					active.old.distribution <- old.distribution
					
					# Use uncorrected distributions for lower bound:
					active.young.distributions <- young.distributions
				
				}
				
				# Special case of last iteration where we want to ensure we draw a date between the constraining medians:
				if(j == resolution) {
					
					# Set old distribution to young half only:
					active.old.distribution <- old.distribution.young.half
					
					# Set default active young distributions:
					active.young.distributions <- young.distributions
					
					# For each young (top bounding) dated node set all young distribution to old half only:
					for(k in 1:length(young.dated.nodes)) active.young.distributions[k, ] <- young.distributions.old.half[k, ]
					
				}
				
				# Modify young and old distributions to remove tails which will always violate node order (young node older than old node and vice versa) to speed up the random draw step below; if part of young distributions are older than the oldest part of the old distribution:
				if(max(sort(as.vector(active.young.distributions))) >= max(sort(active.old.distribution))) {
				
					# Find oldest part of old distribution
					upper.limit <- max(sort(active.old.distribution))
					
					# ????:
					active.young.distributions[which(active.young.distributions >= upper.limit)] <- NA
				
				}
				
				# Vector to store young distribution minima:
				active.young.distribution.mins <- vector(mode="numeric")
				
				# For each young node store minima:
				for(k in 1:length(active.young.distributions[, 1])) active.young.distribution.mins <- c(active.young.distribution.mins, min(sort(active.young.distributions[k, ])))
				
				# If part of old distribution is younger than the youngest part of the youngest young distribution:
				if(min(sort(active.old.distribution)) <= min(active.young.distribution.mins)) {
					
					# Find lower limit (minimum of minima):
					lower.limit <- min(active.young.distribution.mins)
					
					# Remove values less than or equal to the lower limit from the old distribution:
					active.old.distribution <- active.old.distribution[grep(FALSE, active.old.distribution <= lower.limit)]
				
				}
				
				# Draw upper and lower bounds at random from old and young nodes:
				young.random.ages <- vector(mode="numeric")
				
				# For each young node draw a random age from its distributon:
				for(k in 1:length(young.dated.nodes)) young.random.ages <- c(young.random.ages, sort(active.young.distributions[k, ])[ceiling(runif(1, 0, length(sort(active.young.distributions[k, ]))))])
				
				# Repeat for old age:
				old.random.age <- sort(active.old.distribution)[ceiling(runif(1, 0, length(sort(active.old.distribution))))]

				# Ensure old node is older than all young nodes; while old node is younger or equal in age to oldest young nodes:
				while(old.random.age <= max(young.random.ages)) {
					
					# Re-draw ages as above:
					young.random.ages <- vector(mode="numeric")
					
					# For each young node draw a random age from its distributon:
					for(k in 1:length(young.dated.nodes)) young.random.ages <- c(young.random.ages, sort(active.young.distributions[k, ])[ceiling(runif(1, 0, length(sort(active.young.distributions[k, ]))))])
					
					# Repeat for old age:
					old.random.age <- sort(active.old.distribution)[ceiling(runif(1, 0, length(sort(active.old.distribution))))]
				
				}
				
				# Need to add names so can extract data later:
				names(young.random.ages) <- young.dated.nodes
				
				# Find oldest (i.e. constraining) young node for each undated node (vector for storing output):
				young.constraints.node <- young.constraints.age <- vector(mode="numeric")
				
				# For each undated node:
				for(k in 1:length(undated.nodes)) {
					
					# Store young age:
					young.constraints.age <- c(young.constraints.age, max(young.random.ages[strsplit(young.constraints[k], " ")[[1]]]))
					
					# Store young node:
					young.constraints.node <- c(young.constraints.node, as.numeric(names(young.random.ages[strsplit(young.constraints[k], " ")[[1]]])[which(young.random.ages[strsplit(young.constraints[k], " ")[[1]]] == max(young.random.ages[strsplit(young.constraints[k], " ")[[1]]]))[1]]))
				
				}
				
				# Vector for storing node dates:
				node.dates <- rep(NA, length(undated.nodes))
				
				# Date undated nodes using randomisation process:
				for(k in 1:length(unique(young.constraints.node)) ) {
					
					# Find top.node:
					top.node <- unique(young.constraints.node)[k]
					
					# Find target (.e. undated) nodes:
					target.nodes <- undated.nodes[which(young.constraints.node == top.node)]
					
					# Find bounding dates (potential old constraints):
					old.constraining.nodes <- as.numeric(sort(unique(strsplit(paste(old.constraints[match(target.nodes, undated.nodes)], collapse=" "), " ")[[1]])))
					
					# Remove any target nodes present
					if(length(sort(match(target.nodes, old.constraining.nodes))) > 0) old.constraining.nodes <- old.constraining.nodes[-sort(match(target.nodes, old.constraining.nodes))]
					
					# Find youngest old node as actual constraint:
					old.constraining.node <- max(old.constraining.nodes)
					
					# If lower bound is the previously dated old node:
					if(old.constraining.node == old.dated.node) {
						
						# Use that as maximum limit
						max <- old.random.age
					
					# If lower bound is a previously undated node:
					} else {
						
						# Use that as maximum limit
						max <- node.dates[match(old.constraining.node, undated.nodes)]
						
					}
					
					# Set minimum age:
					min <- young.random.ages[as.character(top.node)]
					
					# Get node dates:
					node.dates[match(target.nodes, undated.nodes)] <- sort(runif(length(target.nodes), min = min, max = max), decreasing = TRUE) # Draw node ages from bounding minima and maxima
				
				}
				
				# Store results:
				age.distributions[as.character(undated.nodes), j] <- node.dates
				
			}
			
		}
		
		# Can now remove node "0" (i.e., if root is undated by Hedman method) if present:
		if(root.dated == "no") {
			
			# Remove zero node from distributions:
			age.distributions <- age.distributions[-which(rownames(age.distributions) == "0"), ]
		
			# Remove zero node from estimates:
			age.estimates <- age.estimates[-which(rownames(age.estimates) == 0), ]
	
		}

	}
	
	# Report progress:
	cat("Done\nTidying up and returning results_")

	# Create vector of nodes estimated using Hedman method:
	Hedman.estimated <- rep(1, length(rownames(age.estimates)))
	
	# Update those not estimated using Hedman method:
	if(conservative == TRUE) Hedman.estimated[match(setdiff(rownames(age.estimates), dateable.nodes), rownames(age.estimates))] <- 0
	
	# Add to age estimates:
	age.estimates <- cbind(age.estimates, Hedman.estimated)
	
	# Sort age distributions in advance of picking confidence intervals:
	age.distributions <- t(apply(age.distributions, 1, sort))
	
	# Add median values to age estimates:
	age.estimates[, "Best.guess"] <- apply(age.distributions, 1, median)
	
	# Add upper CI to age estimates:
	age.estimates[, "Best.guess.upper"] <- age.distributions[, ceiling(0.975 * resolution)]
	
	# Add lower CI to age estimates:
	age.estimates[, "Best.guess.lower"] <- age.distributions[, max(c(1, floor(0.025 * resolution)))]

	# Create vectors to store node ages for each branch:
	to.ages <- from.ages <- age.estimates[match(tree$edge[, 1], rownames(age.estimates)), "Best.guess"]
	
	# Get to ages for terminal branches:
	to.ages[which(tree$edge[, 2] <= Ntip(tree))] <- tip.ages[tree$edge[which(tree$edge[, 2] <= Ntip(tree)), 2]]
	
	# Get to ages for internal branches:
	to.ages[which(tree$edge[, 2] > Ntip(tree))] <- age.estimates[match(tree$edge[which(tree$edge[, 2] > Ntip(tree)), 2], rownames(age.estimates)), "Best.guess"]
	
	# Update tree with branch lengths scaled to time using median ages:
	tree$edge.length <- from.ages-to.ages
	
	# Update tree to include root time:
	tree$root.time <- age.estimates[as.character(root.node), "Best.guess"]
	
	# Collate results:
	results <- list(age.estimates, age.distributions, tree)
	
	# Add names:
	names(results) <- c("age.estimates", "age.distributions", "tree")
	
	# Repprt progress
	cat("Done")
	
	# Return output:
	return(results)
	
}

# Alroys metric for the "wobbliness" of diversity curves:
wobble.index <- function(diversity.vector) {
	
	# Calculate wobble index:
	out <- median(abs(log(diversity.vector[2:(length(diversity.vector) - 1)] ** 2 / (diversity.vector[1:(length(diversity.vector) - 2)] * diversity.vector[3:length(diversity.vector)]))))

	# Return result:
	return(out)

}

# Matts intervalFitting function:
intervalFitting <- function(phy, data, divisions, rootAge, interval.names=NULL, data.names=NULL, bounds=NULL, meserr=NULL) {
    
	require(mvtnorm)
	
	# sort is T because sub-functions assume data are in
	# this particular order
    
	name.check  <- 
	function(phy, data, data.names=NULL)
	{
		if(is.null(data.names)) 
		{
			if(is.vector(data))
			data.names=names(data)
			else
			data.names <- rownames(data)
		}
		t <- phy$tip.label
		r1 <- t[is.na(match(t, data.names))]
		r2 <- data.names[is.na(match(data.names, t))]
		
		r <- list(sort(r1), sort(r2))
		
		names(r) <- cbind("Tree.not.data", "Data.not.tree")
		if(length(r1)==0 && length(r2)==0) return("OK")
		else return(r)
	}
	
	treedata <- function(phy, data, data.names=NULL, sort=F, warnings=T)
	{
		
		if(is.vector(data)) data <- as.matrix(data)
		if(is.factor(data)) data <- as.matrix(data)
		if(is.array(data) & length(dim(data))==1) data <- as.matrix(data)
		
		if(is.null(data.names)) {
			if(is.null(rownames(data))) {
				data.names <- phy$tip.label[1:dim(data)[1]]
				if(warnings)
				cat("Warning: no tip labels, order assumed to be the same as in the tree\n")
			} else
			data.names <- rownames(data)
		}
		nc <- name.check(phy, data, data.names)
		if(is.na(nc[[1]][1]) | nc[[1]][1]!="OK") {
			if(length(nc[[1]]!=0)) {
				phy=drop.tip(phy, as.character(nc[[1]]))
				if(warnings) {
					cat("Dropped tips from the tree because there were no matching names in the data:\n")
					print(nc[[1]])
					cat("\n")
				}
			}
			
			if(length(nc[[2]]!=0)) {
				m <- match(data.names, nc[[2]])
				data=as.matrix(data[is.na(m), ])
				data.names <- data.names[is.na(m)]
				if(warnings) {
					cat("Dropped rows from the data because there were no matching tips in the tree:\n")
					print(nc[[2]])
					cat("\n")
				}
			}
		}
		order <- match(data.names, phy$tip.label)	
		
		rownames(data) <- phy$tip.label[order]
		
		if(sort) {
			
			index <- match(phy$tip.label, rownames(data))
			data <- as.matrix(data[index, ])
		}
		
		phy$node.label=NULL
		
		return(list(phy=phy, data=data))
	}
	
	td <- treedata(phy, data, data.names, sort=T)
    
	ntax=length(td$phy$tip.label)
    
	if(is.null(meserr)) {
		me=td$data
		me[]=0
		meserr=me
	} else if(length(meserr)==1) {
		me=td$data
		me[]=meserr
		meserr=me
	} else if(is.vector(meserr)) {
		if(!is.null(names(meserr))) {
			o <- match(rownames(td$data), names(meserr))
			if(length(o)!=ntax) stop("meserr is missing some taxa from the tree")
			meserr <- as.matrix(meserr[o, ])
		} else {
			if(length(meserr)!=ntax) stop("No taxon names in meserr, and the number of taxa does not match the tree")
			me <- td$data
			me[]=meserr
			meserr=me
		}
	} else {
		if(!is.null(rownames(meserr))) {
			o <- match(rownames(td$data), rownames(meserr))
			meserr=meserr[o, ]
		} else {
			if(sum(dim(meserr)!=dim(td$data))!=0)
            stop("No taxon names in meserr, and the number of taxa does not match the tree")
			print("No names in meserr; assuming that taxa are in the same order as tree")
		}
	}
	
	if (is.null(interval.names)) {
        interval.names <- as.character (c(1:(length(divisions)-1)))
	} else if (is.vector(interval.names)) {
        if (length(interval.names) != (length(divisions)-1)) {
            stop ("Number of interval names does not match number of intervals specified by 'divisions'")
            print ("Number of interval names does not match numbr of intervals specified by 'divisions'")
        } else {
            interval.names <- interval.names
        }
    } else {
        stop ("List of interval names is not a vector")
        print ("List of interval names is not a vector")
    }
    
	
	
    
    #--------------------------------
    #---    PREPARE DATA LIST     ---
    #--------------------------------
    ds			 <-  list()
    ds$tree 		 <-  td$phy          # TIP data
    ds$sliceEdge <- timesliceEdge (ds$tree, rootAge, divisions)  #time slice tree once!
    ds$interval.names <- interval.names
    #--------------------------------
    #--- IDENTIFY ANCHORING BIN   ---
    #--------------------------------
	
	ds$anchorBin <- which.max (colSums (ds$sliceEdge))
    
    #--------------------------------
    #--- SET MODEL SPECIFICATIONS ---
    #--------------------------------
    cat("Fitting interval model")
    #-----------------------------
    #---  SET PARAMETER BOUNDS ---
    #-----------------------------
    #---- DEFAULT BOUNDS
    bounds.default			 <- matrix(c(0.00000001, 20, rep (c(0.000001, 10000), (length(interval.names)-1))), nrow=(length(interval.names)), ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("beta", interval.names[c(2:length(interval.names))]);
    colnames(bounds.default) <- c("min", "max")
    
 	#---- USER DEFINED PARAMETER BOUNDS
 	if (is.null(bounds)) {
 		bounds <- bounds.default       # USE DEFAULTS
 	}else{
 		if (class(bounds)!="list"){
 			stop("Please specify user defined parameter bounds as a list()")
 		}else{
 			specified   <- !c(is.null(bounds$beta)
            )
 			bounds.user <- matrix(c(bounds$beta), 
            nrow=sum(specified), ncol=2, byrow=TRUE
            )
 			rownames(bounds.user) <- c("beta")[specified]
   	 		colnames(bounds.user) <- c("min", "max")
            
   	 		#----  SET FINAL SEARCH BOUNDS
 			bounds <- bounds.default
 			bounds[specified, ] <- bounds.user     # Final Bounds
   		} # END if list
   	}  # END user bound if loop
   	#--------------------------------
    #---   APPEND MODEL SETTINGS  ---
    #--------------------------------
  	ds$bounds <- data.frame(t(bounds))
    
  	#--------------------------------
    #---        FIT MODEL         ---
    #--------------------------------
    result <- list()
    for(i in 1:ncol(td$data)) {
    	ds$data=td$data[, i]
    	ds$meserr=meserr[, i]
  		result[[i]] <- fitContinuousModel(ds, print=print)
  		if(!is.null(colnames(td$data))) names(result)[i] <- colnames(td$data)[i] else names(result)[i] <- paste("Trait", i, sep="")
        
  	}
  	result
}

# Matts insert function:
#----------------------------------
#-----   INSERTION FUNCTION   -----
#----------------------------------
insert <- function(v, e, pos){
    return(c(v[1:(pos-1)], e, v[(pos):length(v)]))
}

# Matts fitContinuousModel function:
fitContinuousModel <- function(ds, print=TRUE)
{
	bounds 	 <-  ds$bounds
	n 		 <-  length(ds$data)
    
	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
    
	beta.start <- var(ds$data)/max(branching.times(ds$tree))
    
    
	out         <- NULL
    
	y			 <-  ds$data				# TIP data
	tree		 <-  ds$tree			# Tree
	meserr		 <-  ds$meserr
	n			 <-  length(y)
	nFactors <- ncol (ds$sliceEdge)-1 #specify the number of factors
	sliceEdge <- ds$sliceEdge #sliced branches
  	sliceTree <- ds$tree
	anchorBin <- ds$anchorBin
    
    
    
	#----------------------------------
	#-----       DEFAULT FIT      -----
	#----------------------------------
    
    k <- 2+nFactors
    
    
    start=log(c(beta.start, rep (1, nFactors)))
    lower=log(bounds["min", ])
    upper=log(bounds["max", ])
    
    foo <- function(x) {  #x is a vector of parameters
        
        
        if (anchorBin==1) {
            
            factors <- exp(c(log(1), x[-1])) #place anchor at start
            
        }
        
        else if (anchorBin==length(x)) {
            
            factors <- exp(c(x[-1], log(1))) #place anchor at end
            
        }
        
        else {
            
            x.x <- x [-1]
            
            x.x <- insert (x.x, log(1), anchorBin)
            
            factors <- exp (x.x)
            
        }
        
		
        
        newBranchLengths <- rowSums (aperm(factors*aperm(sliceEdge)))
        sliceTree$edge.length <- newBranchLengths
		
        vcv <- vcv.phylo(sliceTree)
        vv <- exp(x[1])*vcv
        diag(vv) <- diag(vv)+meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        dmvnorm(y, mu, vv, log=T)
    }
    
    control <- list()
    
    control$fnscale <- -1
    
    o <- optim(foo, p=start, lower=lower, upper=upper, method="L", hessian = T, control = control)
    
    #----------------------------------
    #---Caclulate CIs (approximate)----
    #----------------------------------
    
    fisherInfo <- solve(-o$hessian)
    
    propSigma <- sqrt(diag(fisherInfo))
    
    propSigma2 <- propSigma [-1] #removes information for variance parameter
    
    
    
    if (anchorBin==1) {
        
        factors <- c(1, exp(o$par[-1])) #anchor at start
        
        errors <- c(NA, propSigma2)
        
        
        
    }
    
    else if (anchorBin==(nFactors+1)) {
        
        factors <- c(exp(o$par[-1]), 1) #anchor at end
        
        errors <- c(propSigma2, NA)
        
        
    }
    
    else {
        
        factors <- exp(o$par[-1])
        
        factors <- insert (factors, 1, anchorBin)
        
        errors <- insert (propSigma2, NA, anchorBin)
        
    }
    
    
    names (factors) <- ds$interval.names
    
    upperSE <- factors + errors
	
    lowerSE <- factors - errors
    
    upper95 <- factors + 1.96*errors
	
    lower95 <- factors- 1.96*errors
    
    
    
    
    results <- list(lnl=o$value, beta= exp(o$par[1]), factors = factors)
    
    
    
    
	#----------------------------------
	#-----    Collect results     -----
	#----------------------------------
    
    
	results$aic <- 2*k-2*results$lnl
	results$aicc <- 2*k*(n-1)/(n-k-2)-2*results$lnl
	results$k <- k
	results$hessian <- o$hessian
	results$upperSE <- upperSE
	results$lowerSE <- lowerSE
	results$upper95 <- upper95
	results$lower95 <- lower95
	return(results)
    
}

# Matts phylogMean function:
phylogMean <- function(phyvcv, data)
{
	o <- rep(1, length(data))
	ci <- solve(phyvcv)
    
	m1 <- solve(t(o) %*% ci %*% o)
	m2 <- t(o) %*% ci %*% data
    
	return(m1 %*% m2)
}

# Matts timesliceedge function:
timesliceEdge <- function (tree, rootAge, divisions)
{
    require(ape)
    require(phytools)
    
    #Create a new matrix i x j matrix, where i is equal to the number of branches and j is equal to the number of intervals
    
    sliceEdge <- matrix (nrow = nrow(tree$edge), ncol = (length(divisions)-1))
    
    nodeAges <-  -(nodeHeights (tree) - rootAge)
    
    #Run some simple checks
    
    
    
    for (i in 1:nrow (sliceEdge)) {
        
        for (j in 1:ncol(sliceEdge)){
            
            #entirely below interval
            if (nodeAges[i, 2] > divisions[j]){
                
                sliceEdge[i, j] <- 0
                
            }
            
            #entirely above interval
            else if (nodeAges[i, 1] < divisions[j+1]){
                
                sliceEdge[i, j] <- 0
                
            }
            
            #entirely within interval
            else if (nodeAges[i, 1] < divisions[j] && nodeAges[i, 2] > divisions[j+1]){
                
                sliceEdge[i, j] <- nodeAges[i, 1] - nodeAges[i, 2]
                
            }
            
            #extends through entire interval
            else if (nodeAges[i, 1] >= divisions[j] && nodeAges[i, 2] <= divisions[j+1]){
                sliceEdge[i, j] <- divisions[j]-divisions[j+1]
            } else {
                #includes beginning of interval, but ends in interval
                if (nodeAges[i, 2] > divisions[j+1]){
                    sliceEdge[i, j] <- divisions[j]-nodeAges[i, 2]
                } else { #includes end of interval, but begins in interval
                    sliceEdge[i, j] <- nodeAges[i, 1]-divisions[j+1]
                }
            }
        }
    }
    return(sliceEdge)
}

collapse.clade <- function(n, tree)
{
    # If tree has no branch lengths make all branch lengths equal one:
	if(is.null(tree$edge.length)) tree$edge.length <- rep(1, length(tree$edge[, 1]))
    
    # If there are zero length branches make all branches equal to one:
    #if(min(tree$edge.length) == 0) tree$edge.length <- rep(1, length(tree$edge[, 1]))
    
    # Set the collapse cutoff branch length as less than the smallest branch length:
	tol <- min(tree$edge.length)/2
    
    # Set all branches ionside clade we want to collapse to zero:
	for(i in n) tree$edge.length[GetDescendantEdges(i, tree)] <- 0
    
    # Collapse those branches:
	tree <- di2multi(tree, tol)
    
    # Return the collapsed tree:
	return(tree)
}

# Function to calculate all possible combinations for variables given in var.names:
var.combns <- function(var.names) {
	require(utils)
	combos <- vector(mode="character")
	for(i in 1:(length(var.names)-1)) {
		mat <- apply(combn(var.names, i), 1, paste)
		if(is.matrix(mat) == FALSE) mat <- as.matrix(mat)
		for(j in 1:length(mat[, 1])) combos <- c(combos, paste(mat[j, ], collapse=" + "))
	}
	combos <- c(combos, paste(var.names, collapse=" + "))
	return(combos)
}

# Function to reinsert safely deleted taxa:
str.reinsert <- function(treefile.in, treefile.out, str.list, multi.placements="exclude")
{
	
	# Ensure str list is formatted as characters:
	str.list[, "Junior"] <- as.character(str.list[, "Junior"])
	str.list[, "Senior"] <- as.character(str.list[, "Senior"])
	str.list[, "Rule"] <- as.character(str.list[, "Rule"])
	
	# Read in tree file as text:
	text <- scan(treefile.in, what="\n", quiet=TRUE)
	
	# Resort str list by Junior taxon:
	str.list <- str.list[order(str.list[, "Junior"]), ]
	
	# Find number of seniors for each junior:
	names.and.numbers <- rle(as.character(str.list[, "Junior"]))
	
	# Make list of taxa that have a single senior:
	single.replacements <- names.and.numbers$values[grep(TRUE, names.and.numbers$lengths == 1)]
	
	# Collapse to just those not in a polytomy with other taxa in the str list:
	single.replacements <- single.replacements[is.na(match(single.replacements, str.list[, "Senior"]))]
	
	# For taxa that have a single replacement:
	if(length(single.replacements) > 0) {
		
		# For each single replacement taxon:
		for(i in 1:length(single.replacements)) {
			
			# Reinsert into tree next to its senior:
			text <- gsub(str.list[match(single.replacements[i], str.list[, "Junior"]), "Senior"], paste("(", paste(str.list[match(single.replacements[i], str.list[, "Junior"]), c("Junior", "Senior")], collapse=", "), ")", sep=""), text)
			
		}
		
		# Remove single replacement taxa from str list:
		str.list <- str.list[-match(single.replacements, str.list[, "Junior"]), ]
		
		# Update names and numbers now taxa have been removed:
		names.and.numbers <- rle(as.character(sort(str.list[, "Junior"])))
	}
	
	# Vector to store taxa that only occur in a single polytomy:
	polytomy.taxa <- vector(mode="character")
	
	# For each taxon:
	for(i in 1:length(names.and.numbers$values)) {
		
		# Get taxon name:
		taxon.name <- names.and.numbers$values[i]
		
		# Find its seniors:
		seniors <- str.list[grep(TRUE, str.list[, "Junior"] == taxon.name), "Senior"]
		
		# If its seniors, except for one, are all also juniors then record it:
		# (This finds taxa that only exist in a single polytomy in the original tree)
		if(length(grep(TRUE, is.na(match(seniors, str.list[, "Junior"])))) <= 1) polytomy.taxa <- c(polytomy.taxa, taxon.name)
	}
	
	# If there are taxa that only occur in a single polytomy:
	if(length(polytomy.taxa) > 0) {
		
		# Reorder from most seniors to least:
		taxa.to.delete <- polytomy.taxa <- polytomy.taxa[order(names.and.numbers$lengths[match(polytomy.taxa, names.and.numbers$values)], decreasing=TRUE)]
		
		# Whilst there are still polytomous taxa in the list:
		while(length(polytomy.taxa) > 0) {
			
			# Get taxon name with most juniors:
			taxon.name <- polytomy.taxa[1]
			
			# Find all of its seniors:
			seniors <- str.list[grep(TRUE, str.list[, "Junior"] == taxon.name), "Senior"]
			
			# Find senior taxon (already in tree):
			senior.taxon <- seniors[grep(TRUE, is.na(match(seniors, str.list[, "Junior"])))]
			
			# Replace senior with all juniors in polytomy:
			text <- gsub(senior.taxon, paste("(", paste(sort(c(taxon.name, seniors)), collapse=", "), ")", sep=""), text)
			
			# Remove taxa just dealt with from polytomy.taxa:
			polytomy.taxa <- polytomy.taxa[-sort(match(c(taxon.name, seniors), polytomy.taxa))]
		}
		
		# Trims str list down to remaining taxa:
		for(i in 1:length(taxa.to.delete)) str.list <- str.list[-grep(TRUE, str.list[, "Junior"] == taxa.to.delete[i]), ]
	}
	
	# Only keep going if there are taxa still to reinsert:
	if(length(str.list[, 1]) > 0) {
		
		# If the user wishes to reinsert remaining taxa at random:
		if(multi.placements == "random") {
			
			# List unique juniors:
			unique.juniors <- rle(str.list[, "Junior"])$values[order(rle(str.list[, "Junior"])$lengths)]
			
			# For each junior taxon remaining:
			for(i in 1:length(unique.juniors)) {
				
				# Isolate junior taxon name:
				junior <- unique.juniors[i]
				
				# Get a senior for each tree:
				seniors <- sample(str.list[grep(TRUE, str.list[, "Junior"] == junior), "Senior"], length(text), replace=TRUE)
				
				# Remove junior from str list:
				str.list <- str.list[-grep(TRUE, str.list[, "Junior"] == junior), ]
				
				# Make replacements:
				replacements <- paste("(", paste(junior, seniors, sep=", "), ")", sep="")
				
				# For each tree:
				for(j in 1:length(text)) {
					
					# Isolate tree:
					tree.text <- text[j]
					
					# Case if senior ends in a parenthesis:
					if(length(grep(paste(seniors[j], ")", sep=""), tree.text)) > 0) {
						
						# Replace senior with combination of junior and senior:
						tree.text <- gsub(paste(seniors[j], ")", sep=""), paste(replacements[j], ")", sep=""), tree.text)
					}
					
					# Case if senior ends in a comma:
					if(length(grep(paste(seniors[j], ", ", sep=""), tree.text)) > 0) {
						
						# Replace senior with combination of junior and senior:
						tree.text <- gsub(paste(seniors[j], ", ", sep=""), paste(replacements[j], ", ", sep=""), tree.text)
					}
					
					# Update tree text:
					text[j] <- tree.text
					
				}
				
			}
			
			# Set uninserted list to null (i.e. empty):
			uninserted <- NULL
		}
		
		# If the user wishes to exclude remaining taxa:
		if(multi.placements == "exclude") {
			
			# Set uninserted list to remaining juniors:
			uninserted <- unique(str.list[, "Junior"])
		}
		
	# If taxa all have a single reinsertion position:
	} else {
		
		# Set uninserted list to null (i.e. empty):
		uninserted <- NULL
	}
	
	# Write out tree file:
	write(text, treefile.out)
	
	# Compile output:
	output <- list(uninserted)
	
	# Name output:
	names(output) <- "unreinserted.taxa"
	
	# Return output:
	return(output)
}

GetNBifurcatingResolutions <- function(tree) {
	
	# Require the ape library:
	require(ape)
	
	# Get a list of the number of descendant edges for each node:
	Nfurcations <- rle(sort(tree$edge[, 1]))$lengths
	
	# Reduce this list to just the polytomies (nodes with three or more descendant edges):
	polytomies <- Nfurcations[grep(TRUE, Nfurcations > 2)]
	
	# If there is at least one polytomy:
	if(length(polytomies) > 0) {
		
		# Work through each polytomy:
		for(i in 1:length(polytomies)) {
			
			# Record the number of descendnat edges for the ith polytomy:
			n <- polytomies[i]
			
			# If tree is rooted:
			if(is.rooted(tree)) {
				
				# Calculate the number of possible rooted bifurcations for the ith polytomy (equation 1 in Felsenstein 1978):
				polytomies[i] <- factorial((2 * n) - 3) / ((2 ^ (n - 2)) * factorial(n - 2))
				
			# If tree is unrooted:
			} else {
				
				# Calculate the number of possible unrooted bifurcations for the ith polytomy (equation from Casey Dunns slides):
				polytomies[i] <- factorial((2 * n) - 5) / ((2 ^ (n - 3)) * factorial(n - 3))
				
			}
			
		}
		
		# Find product of all polytomies (total number of bifurcating resolutions):
		out <- prod(polytomies)
		
	# If there are no polytomies:
	} else {
		
		# Warn the user:
		print("Tree is already fully bifurcating")
		
		# Record the number of bifurcating resolutions (has to be one):
		out <- 1
		
	}
	
	# Return the number of bifurcating resolutions:
	return(out)
	
}

# Compile functions; load compiler library:
library(compiler) 

# For each object in memory:
for(i in 1:length(objects())) {
	
	# If the object is a function (i.e. do not try to compile non-function objects):
    if(is.function(get(objects()[i]))) {
		
		# Compile function and reassign it to itself:
        assign(objects()[i], cmpfun(get(objects()[i]), options=c(suppressAll=TRUE)))
		
    }
	
}