## ---- include = FALSE---------------------------------------------------------
library(knitr)
set.seed(1)
is_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")
if (!is_CRAN) {
   options(mc.cores = parallel::detectCores())
} else {
  knitr::opts_chunk$set(eval = FALSE)
  knitr::knit_hooks$set(evaluate.inline = function(x, envir) x)
}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache.lazy = FALSE  # https://github.com/yihui/knitr/issues/572
)
options(digits=4)
options(scipen=2)

## ---- message=FALSE-----------------------------------------------------------
library(rstan)
library(pcFactorStan)
library(loo)
library(qgraph)
library(ggplot2)
library(Matrix)

## ---- results='hide'----------------------------------------------------------
head(phyActFlowPropensity)

## ---- results='asis', echo=FALSE----------------------------------------------
kable(head(phyActFlowPropensity))

## -----------------------------------------------------------------------------
dl <- prepData(phyActFlowPropensity[,c(paste0('pa',1:2), 'skill')])
dl$varCorrection <- 5.0

## ----pcStan, message=FALSE, results='hide', cache=TRUE------------------------
fit1 <- pcStan("unidim_adapt", data=dl)

## ----pcStanDiag1, cache=TRUE--------------------------------------------------
check_hmc_diagnostics(fit1)

## ----pcStanDiag2, cache=TRUE--------------------------------------------------
allPars <- summary(fit1, probs=c())$summary
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ----skill, cache=TRUE--------------------------------------------------------
theta <- summary(fit1, pars=c("theta"), probs=c())$summary[,'mean']

ggplot(data.frame(x=theta, activity=dl$nameInfo$pa, y=0.47)) +
  geom_point(aes(x=x),y=0) +
  geom_text(aes(label=activity, x=x, y=y),
            angle=85, hjust=0, size=2,
            position = position_jitter(width = 0, height = 0.4)) + ylim(0,1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## -----------------------------------------------------------------------------
s50 <- summary(fit1, pars=c("scale"), probs=c(.5))$summary[,'50%']
print(s50)

## ---- results='hide', include = FALSE-----------------------------------------
rm(fit1)  # free up some memory

## ----calibrateItems, message=FALSE, results='hide', cache=TRUE----------------
result <- calibrateItems(phyActFlowPropensity, iter=1000L)

## ---- results='hide'----------------------------------------------------------
print(result)

## ---- results='asis', echo=FALSE----------------------------------------------
kable(result)

## ----covarianceData, cache=TRUE-----------------------------------------------
pafp <- phyActFlowPropensity
excl <- match(c('goal1','feedback1'), colnames(pafp))
pafp <- pafp[,-excl]
dl <- prepData(pafp)
dl$scale <- result[match(dl$nameInfo$item, result$item), 'scale']

## ----covariance, message=FALSE, results='hide', cache=TRUE--------------------
fit2 <- pcStan("correlation", data=dl, include=FALSE, pars=c('rawTheta', 'rawThetaCorChol'))

## ----covarianceDiag1, cache=TRUE----------------------------------------------
check_hmc_diagnostics(fit2)

allPars <- summary(fit2, probs=0.5)$summary
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ----covarianceDiag2, cache=TRUE----------------------------------------------
head(allPars[order(allPars[,'sd']),])

## ----covarianceDiag3, cache=TRUE----------------------------------------------
allPars <- allPars[allPars[,'sd'] > 1e-6,]
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ----corPlot, cache=TRUE------------------------------------------------------
corItemNames <- dl$nameInfo$item
tc <- summary(fit2, pars=c("thetaCor"), probs=c(.1,.5,.9))$summary
tcor <- matrix(tc, length(corItemNames), length(corItemNames))
tcor[sign(tc[,'10%']) != sign(tc[,'90%'])] <- 0  # delete faint edges
dimnames(tcor) <- list(corItemNames, corItemNames)
tcor <- nearPD(tcor, corr=TRUE)$mat

qgraph(tcor, layout = "spring", graph = "cor", labels=colnames(tcor),
       legend.cex = 0.3,
       cut = 0.3, maximum = 1, minimum = 0, esize = 20,
       vsize = 7, repulsion = 0.8, negDashed=TRUE, theme="colorblind")

## ----responseCurves, cache=TRUE-----------------------------------------------
df <- responseCurve(dl, fit2,
  item=setdiff(dl$nameInfo$item, c('spont','control','evaluated','waiting')),
  responseNames=c("much more","somewhat more", 'equal',
                  "somewhat less", "much less"))
ggplot(df) +
  geom_line(aes(x=worthDiff, y=prob, color=response,linetype=response,
                group=responseSample), size=.2, alpha=.2) +
  xlab("difference in latent worths") + ylab("probability") +
  ylim(0,1) + facet_wrap(~item, scales="free_x") +
  guides(color=guide_legend(override.aes=list(alpha = 1, size=1)))

## ----factorData1, cache=TRUE--------------------------------------------------
pafp <- pafp[,c(paste0('pa',1:2),
             setdiff(corItemNames, c('spont','control','evaluated','waiting')))]
pafp <- normalizeData(filterGraph(pafp))
dl <- prepCleanData(pafp)
dl <- prepSingleFactorModel(dl)
dl$scale <- result[match(dl$nameInfo$item, result$item), 'scale']

## ---- results='hide', include = FALSE-----------------------------------------
rm(fit2)  # free up some memory

## ----factor, message=FALSE, results='hide', cache=TRUE------------------------
fit3 <- pcStan("factor1_ll", data=dl, include=FALSE,
               pars=c('rawUnique', 'rawUniqueTheta', 'rawPerComponentVar',
	       'rawFactor', 'rawLoadings', 'rawFactorProp', 'rawThreshold',
         'rawPathProp', 'rawCumTh'))

## ----factorDiag1, cache=TRUE--------------------------------------------------
check_hmc_diagnostics(fit3)

interest <- c("threshold", "alpha", "pathProp", "factor", "residualItemCor", "lp__")

allPars <- summary(fit3, pars=interest)$summary
allPars <- allPars[allPars[,'sd'] > 1e-6,]
print(min(allPars[,'n_eff']))
print(max(allPars[,'Rhat']))

## ----factorLoo, cache=TRUE----------------------------------------------------
options(mc.cores=1)  # otherwise loo consumes too much RAM
kThreshold <- 0.1
l1 <- toLoo(fit3)
print(l1)

## -----------------------------------------------------------------------------
pa11 <- levels(filterGraph(pafp, minDifferent=11L)$pa1)
ot <- outlierTable(dl, l1, kThreshold)
ot <- subset(ot, pa1 %in% pa11 & pa2 %in% pa11)

## ---- results='hide'----------------------------------------------------------
print(ot[1:6,])

## ---- results='asis', echo=FALSE----------------------------------------------
kable(ot[1:6,], row.names=TRUE)

## ---- results='hide'----------------------------------------------------------
xx <- which(ot[,'pa1'] == 'mountain biking' & ot[,'pa2'] == 'climbing' & ot[,'item'] == 'predict' & ot[,'pick'] == -2)

## ---- results='asis', echo=FALSE----------------------------------------------
if (length(xx) == 0) {
  xx <- 1
  warning("Can't find outlier")
}
kable(ot[xx,,drop=FALSE], row.names=TRUE)

## -----------------------------------------------------------------------------
pafp[pafp$pa1 == ot[xx,'pa1'] & pafp$pa2 == ot[xx,'pa2'],
     c('pa1','pa2', as.character(ot[xx,'item']))]

## ----outlier, cache=TRUE------------------------------------------------------
loc <- sapply(ot[xx,c('pa1','pa2','item')], unfactor)
exam <- summary(fit3, pars=paste0("theta[",loc[paste0('pa',1:2)],
                          ",", loc['item'],"]"))$summary
rownames(exam) <- c(as.character(ot[xx,'pa1']), as.character(ot[xx,'pa2']))

## ---- results='asis', echo=FALSE----------------------------------------------
#exam <- data.frame(mean=c(0,0), '2.5%'=c(0,0), '97.5%'=c(0,0))
kable(exam)

## -----------------------------------------------------------------------------
sum(c(pafp$pa1 == ot[xx,'pa1'], pafp$pa2 == ot[xx,'pa1']))
sum(c(pafp$pa1 == ot[xx,'pa2'], pafp$pa2 == ot[xx,'pa2']))

## ----pathProp, cache=TRUE-----------------------------------------------------
pi <- parInterval(fit3, 'pathProp', dl$nameInfo$item, label='item')
pi <- pi[order(abs(pi$M)),]

ggplot(pi) +
  geom_vline(xintercept=0, color="green") +
  geom_jitter(data=parDistributionFor(fit3, pi),
              aes(value, item), height = 0.35, alpha=.05) +
  geom_segment(aes(y=item, yend=item, x=L, xend=U),
               color="yellow", alpha=.5) +
  geom_point(aes(x=M, y=item), color="red", size=1) +
  theme(axis.title.y=element_blank())

## ----activities, cache=TRUE---------------------------------------------------
pick <- paste0('factor[',match(pa11, dl$nameInfo$pa),',1]')
pi <- parInterval(fit3, pick, pa11, label='activity')
pi <- pi[order(pi$M),]

ggplot(pi) +
  geom_vline(xintercept=0, color="green") +
  geom_jitter(data=parDistributionFor(fit3, pi, samples=200),
              aes(value, activity), height = 0.35, alpha=.05) +
  geom_segment(aes(y=activity, yend=activity, x=L, xend=U),
               color="yellow", alpha=.5) +
  geom_point(aes(x=M, y=activity), color="red", size=1) +
  theme(axis.title.y=element_blank())

## ----residualItemCor, cache=TRUE----------------------------------------------
m <- matrix(apply(expand.grid(r=1:dl$NITEMS, c=1:dl$NITEMS), 1,
      function(x) paste0("residualItemCor[",x['r'],",",x['c'],"]")),
      dl$NITEMS, dl$NITEMS)
n <- matrix(apply(expand.grid(r=dl$nameInfo$item, c=dl$nameInfo$item), 1,
                  function(x) paste0(x['r'],":",x['c'])),
            dl$NITEMS, dl$NITEMS)
pi <- parInterval(fit3, m[lower.tri(m)], n[lower.tri(n)], label='cor')
pi <- pi[abs(pi$M) > .08,]
pi <- pi[order(-abs(pi$M)),]
ggplot(pi) +
  geom_vline(xintercept=0, color="green") +
  geom_jitter(data=parDistributionFor(fit3, pi, samples=800),
              aes(value, cor), height = 0.35, alpha=.05) +
  geom_segment(aes(y=cor, yend=cor, x=L, xend=U),
               color="yellow", alpha=.5) +
  geom_point(aes(x=M, y=cor), color="red", size=1) +
  theme(axis.title.y=element_blank())

## -----------------------------------------------------------------------------
sessionInfo()

