## ---- include = FALSE----------------------------------------------------
set.seed(1)
is_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")
if (!is_CRAN) {
   options(mc.cores = parallel::detectCores())
} else {
   q()  # takes too long, pretend everything is fine
}
library(knitr)
library(ggplot2)
library(reshape2)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE------------------------------------------------------
library(rstan)
library(pcFactorStan)

## ---- results='hide'-----------------------------------------------------
head(phyActFlowPropensity)

## ---- results='asis', echo=FALSE-----------------------------------------
kable(head(phyActFlowPropensity))

## ------------------------------------------------------------------------
dl <- prepData(phyActFlowPropensity[,c(paste0('pa',1:2), 'skill')])
dl$varCorrection <- 2.0

## ----pcStan, message=FALSE, results='hide', cache=TRUE-------------------
fit1 <- pcStan("unidim_adapt", data=dl)

## ------------------------------------------------------------------------
check_hmc_diagnostics(fit1)

## ------------------------------------------------------------------------
allPars <- summary(fit1, probs=c())$summary
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ----skill, cache=TRUE---------------------------------------------------
library(ggplot2)

theta <- summary(fit1, pars=c("theta"), probs=c())$summary[,'mean']
palist <- levels(filterGraph(phyActFlowPropensity)$pa1)

ggplot(data.frame(x=theta, activity=palist, y=0.47)) +
  geom_point(aes(x=x),y=0) +
  geom_text(aes(label=activity, x=x, y=y),
            angle=85, hjust=0, size=2,
            position = position_jitter(width = 0, height = 0.4)) + ylim(0,1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## ------------------------------------------------------------------------
s50 <- summary(fit1, pars=c("scale"), probs=c(.5))$summary[,'50%']
print(s50)

## ----calibrateItems, message=FALSE, results='hide', cache=TRUE-----------
result <- calibrateItems(phyActFlowPropensity, iter=1000L)


## ---- results='hide'-----------------------------------------------------
print(result)

## ---- results='asis', echo=FALSE-----------------------------------------
kable(result)

## ------------------------------------------------------------------------
excl <- match(c('goal1','feedback1'), result$item)
median(result[-excl,'scale'])

## ------------------------------------------------------------------------
pafp <- phyActFlowPropensity
excl <- match(c('goal1','feedback1'), colnames(pafp))
pafp <- pafp[,-excl]
dl <- prepData(pafp)
dl$scale <- median(result[-excl,'scale'])

## ----covariance, message=FALSE, results='hide', cache=TRUE---------------
fit2 <- pcStan("covariance", data=dl, include=FALSE, pars=c('rawTheta', 'rawThetaCorChol'))

## ------------------------------------------------------------------------
check_hmc_diagnostics(fit2)

allPars <- summary(fit2, probs=0.5)$summary
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ------------------------------------------------------------------------
head(allPars[order(allPars[,'sd']),])

## ------------------------------------------------------------------------
excl <- match(paste0('thetaCor[',1:dl$NITEMS,',', 1:dl$NITEMS,']'), rownames(allPars))
allPars <- allPars[-excl,]
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ------------------------------------------------------------------------
itemNames <- colnames(pafp)[-(1:2)]
tc <- summary(fit2, pars=c("thetaCor"), probs=c(.5))$summary[,'50%']
tcor <- matrix(tc, length(itemNames), length(itemNames))
dimnames(tcor) <- list(itemNames, itemNames)

library(qgraph)
qgraph(tcor, layout = "spring", graph = "cor", labels=colnames(tcor),
       legend.cex = 0.3,
       cut = 0.3, maximum = 1, minimum = 0, esize = 20,
       vsize = 7, repulsion = 0.8)

## ------------------------------------------------------------------------
itemNames <- setdiff(itemNames, c('spont','control','evaluated','waiting'))
pafp <- pafp[,c(paste0('pa',1:2), itemNames)]
dl <- prepData(pafp)
dl$scale <- median(result[-excl,'scale'])  # close enough

## ----factor, message=FALSE, results='hide', cache=TRUE-------------------
fit3 <- pcStan("factor", data=dl, include=FALSE,
               pars=c('rawUnique', 'rawUniqueTheta', 'rawFactor', 'rawLoadings'))

## ------------------------------------------------------------------------
check_hmc_diagnostics(fit3)

allPars <- summary(fit3, probs=0.5)$summary
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ----factor2, message=FALSE, results='hide', cache=TRUE------------------
fit3 <- pcStan("factor", data=dl, include=FALSE, iter=4000L,
               pars=c('rawUnique', 'rawUniqueTheta', 'rawFactor', 'rawLoadings'))

## ------------------------------------------------------------------------
check_hmc_diagnostics(fit3)

allPars <- summary(fit3, probs=0.5)$summary
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ----factorProp, cache=TRUE----------------------------------------------
library(reshape2)

prop <- summary(fit3, pars='factorProp', probs=c(0.1, 0.5, 0.9))$summary[,c('10%','50%','90%')]
colnames(prop) <- paste0('q',c(10,50,90))
prop <- as.data.frame(prop)
prop$item <- factor(itemNames)

rawProp <- extract(fit3, pars=c("factorProp"))[[1]]
colnames(rawProp) <- itemNames
rawProp <- rawProp[sample.int(nrow(rawProp), 500),]
rawPropTall <- melt(rawProp, variable.name='item')
colnames(rawPropTall)[2] <- c('item')

ggplot() +   geom_vline(xintercept=0, color="green") +
  geom_jitter(data=rawPropTall, aes(value, item), height = 0.35, alpha=.05) +
  geom_segment(data=prop, aes(y=item, yend=item, x=q10, xend=q90),
               color="yellow", alpha=.5) +
  geom_point(data=prop, aes(x=q50, y=item), color="red", size=1) +
  theme(axis.title.y=element_blank())

## ----responseCurve-------------------------------------------------------

thresholdVec <- summary(fit3, pars='threshold', probs=c())$summary[,'mean']
thresholds <- matrix(thresholdVec, nrow=2, dimnames=list(NULL, itemNames))

softmax <- function(y) exp(y) / sum(exp(y))
draw1 <- function(scale, th, item) {
  tdiff <- seq(-2.5/scale, 2.5/scale, .05/scale)
  gr <- expand.grid(tdiff=tdiff, category=c("much more","somewhat more", 'equal',
    "somewhat less", "much less"), p=NA, item=item)
  gg <- matrix(c(0, rev(cumsum(th)), -cumsum(th)), ncol=5, nrow=length(tdiff), byrow=TRUE)
  gg[,2:5] <-  (gg[,2:5] + scale * tdiff)
  gg <- t(apply(gg, 1, cumsum))
  gg <- t(apply(gg, 1, softmax))
  for (lev in 1:length(levels(gr$category))) {
    gr[gr$category == levels(gr$category)[lev],'p'] <- gg[,lev]
  }
  geom_line(data=gr, aes(x=tdiff,y=p,color=category,linetype=category), size=.1, alpha=.2)
}

rawThreshold <- extract(fit3, pars=c("threshold"))[[1]]
rawThreshold <- rawThreshold[sample.int(nrow(rawThreshold), 50),]

pl <- ggplot() + xlab("difference in latent score (logits)") + ylab("probability") +
  ylim(0,1) + facet_wrap(~item)
for (cx in 1:nrow(rawThreshold)) {
    for (ix in 1:length(itemNames)) {
       pl <- pl + draw1(dl$scale, rawThreshold[cx,c(ix*2-1,ix*2)], itemNames[ix])
    }
}
print(pl)

## ---- results='hide'-----------------------------------------------------
print(thresholds)

## ---- results='asis', echo=FALSE-----------------------------------------
kable(thresholds)

## ----activities, cache=TRUE----------------------------------------------
orig <- filterGraph(pafp)$pa1
pa11 <- filterGraph(pafp, minDifferent=11L)$pa1

fs <- summary(fit3, pars='factor', probs=c(0.1, 0.5, 0.9))$summary[,c('10%','50%','90%')]
colnames(fs) <- paste0('q',c(10,50,90))
fs <- fs[match(levels(pa11), levels(orig)),]
fs <- as.data.frame(fs)
fs$activity <- levels(pa11)
fs$activity <- factor(fs$activity, levels=levels(pa11)[order(fs$q50)])

rawFs <- extract(fit3, pars=c("factor"))[[1]]
rawFs <- rawFs[,match(levels(pa11), levels(orig))]
colnames(rawFs) <- levels(pa11)
rawFs <- rawFs[sample.int(nrow(rawFs), 500),]
rawFsTall <- melt(rawFs, variable.name='activity')
colnames(rawFsTall)[2] <- c('activity')
rawFsTall$activity <- factor(as.character(rawFsTall$activity),
                         levels=levels(pa11)[order(fs$q50)])

ggplot() +   geom_vline(xintercept=0, color="green") +
  geom_jitter(data=rawFsTall, aes(value, activity), height = 0.35, alpha=.05) +
  geom_segment(data=fs, aes(y=activity, yend=activity, x=q10, xend=q90),
               color="yellow", alpha=.5) +
  geom_point(data=fs, aes(x=q50, y=activity), color="red", size=1) +
  theme(axis.title.y=element_blank())


## ------------------------------------------------------------------------
sessionInfo()

