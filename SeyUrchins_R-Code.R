setwd("C:/Users/jan/OneDrive/Box/PhD Feedbacks/3 - Data Collection/Data 2018/R Repos 2019")
urchins2017 <- read.csv("Urchins2017.csv")
str(urchins2017)

#################### Load packages

library(dplyr)
library(lattice)
source("HighstatLibV10.R")
library(lme4)
library(pscl)
library(glmmTMB)
library(bbmle) ## for AICtab
library(ggplot2)
library(tibble)
library(gridExtra)
library(emmeans)

boxplot(urchins ~ habitat, data = urchins2017)
boxplot(urchins ~ location, data = urchins2017)

# Mean-center covariates

urchins2017$complexity.sc <- scale(urchins2017$complexity)
urchins2017$cca.sc <- scale(urchins2017$cca)
urchins2017$sand.sc <- scale(urchins2017$sand)
urchins2017$rubble.sc <- scale(urchins2017$rubble)
urchins2017$branchdead.sc <- scale(urchins2017$branchdead)
urchins2017$macroalgae.sc <- scale(urchins2017$macroalgae)
urchins2017$coral.sc <- scale(urchins2017$coral)
urchins2017$habitatgra.sc <- scale(urchins2017$habitatgra)
urchins2017$habitatcar.sc <- scale(urchins2017$habitatcar)
urchins2017$habitatpat.sc <- scale(urchins2017$habitatpat)
urchins2017$herbivores.sc <- scale(urchins2017$herbivores)
urchins2017$invertivores.sc <- scale(urchins2017$invertivores)

head(urchins2017)

#################### Data exploration

# Outliers

par(mfrow = c(1, 2))
boxplot(urchins2017$urchins, 
        main = "Total Count")
dotchart(urchins2017$urchins, 
         xlab = "Range of data", 
         ylab = "Order of the data")

# Collinearity X

MyVar <- c("complexity.sc", "cca.sc", "sand.sc", "rubble.sc", "branchdead.sc", "macroalgae.sc", 
           "coral.sc", "habitatpat.sc", "herbivores.sc", "invertivores.sc")
pairs(urchins2017[,MyVar], 
      lower.panel = panel.cor)
corvif(urchins2017[,MyVar])

# Relationship X & Y

Myxyplot(urchins2017,MyVar, "urchins")

# Interactions

Myxyplot(urchins2017,MyVar, "rubble.sc")
Myxyplot(urchins2017,MyVar, "sand.sc")
Myxyplot(urchins2017,MyVar, "macroalgae.sc")
Myxyplot(urchins2017,MyVar, "branchdead.sc")
Myxyplot(urchins2017,MyVar, "cca.sc")
Myxyplot(urchins2017,MyVar, "complexity.sc")
Myxyplot(urchins2017,MyVar, "coral.sc")
Myxyplot(urchins2017,MyVar, "habitatpat.sc")
Myxyplot(urchins2017,MyVar, "herbivores.sc")
Myxyplot(urchins2017,MyVar, "invertivores.sc")

# Zero inflation
sum(urchins2017$urchins == 0)
100 * sum(urchins2017$urchins == 0) / nrow(urchins2017)

################ Mixed Effects ZIP

formula1 <- urchins ~ macroalgae.sc + complexity.sc + coral.sc + habitatpat.sc + 
  herbivores.sc + invertivores.sc + 
  macroalgae.sc:complexity.sc + habitatpat.sc:complexity.sc + (1|location)

fit_zipoisson <- glmmTMB(formula1,
                         data=urchins2017,
                         ziformula=~1,
                         family=poisson)

summary(fit_zipoisson)

fit_zinbinom <- update(fit_zipoisson,family=nbinom2)

AICtab(fit_zipoisson,fit_zinbinom)

#model validation with Pearson residuals not possible with package glmmTMB

plot(resid(fit_zinbinom)~predict(fit_zinbinom))
plot(resid(fit_zinbinom)~urchins2017$urchins)

#worrying patterns in both residual plots

################ Fixed Effects ZIP/ZINB

formula3 <- urchins ~ macroalgae.sc + complexity.sc + coral.sc + habitatpat.sc + 
  herbivores.sc + invertivores.sc + 
  macroalgae.sc:complexity.sc + habitatpat.sc:complexity.sc

formula4 <- urchins ~ macroalgae.sc + complexity.sc + coral.sc + habitatpat.sc + 
  macroalgae.sc:complexity.sc

formula5 <- urchins ~ macroalgae.sc + complexity.sc + habitatpat.sc + macroalgae.sc:complexity.sc

modelZIP3 <- zeroinfl(formula = formula3, dist = "poisson", data = urchins2017)

summary(modelZIP3)

modelZINB3 <- zeroinfl(formula = formula1, dist = "negbin", link = "logit", data = urchins2017)

AIC(modelZIP3, modelZINB3)

step(modelZINB3)
step(modelZIP3)

modelZINB3 <- zeroinfl(formula = formula3, dist = "negbin", link = "logit", data = urchins2017)
modelZINB4 <- zeroinfl(formula = formula4, dist = "negbin", link = "logit", data = urchins2017)
modelZINB5 <- zeroinfl(formula = formula5, dist = "negbin", link = "logit", data = urchins2017)

step(modelZINB3)
step(modelZINB4)
step(modelZINB5)

AIC(fit_zipoisson,fit_zinbinom, modelZIP3, modelZINB3, modelZINB4, modelZINB5)

summary(modelZINB3)
summary(modelZINB4)
summary(modelZINB5)

################### Model validation

E1 <- resid(modelZINB3, type = "pearson")
N  <- nrow(urchins2017)
p  <- length(coef(modelZINB3))
sum(E1^2) / (N - p)

E1 <- resid(modelZINB4, type = "pearson")
N  <- nrow(urchins2017)
p  <- length(coef(modelZINB4))
sum(E1^2) / (N - p)

E1 <- resid(modelZINB5, type = "pearson")
N  <- nrow(urchins2017)
p  <- length(coef(modelZINB5))
sum(E1^2) / (N - p)

par(mfrow = c(1, 3))
plot(predict(modelZINB3), resid(modelZINB3, type = "pearson"))
plot(predict(modelZINB4), resid(modelZINB4, type = "pearson"))
plot(predict(modelZINB5), resid(modelZINB5, type = "pearson"))
plot(resid(modelZINB3, type = "pearson"))
plot(resid(modelZINB4, type = "pearson"))
plot(resid(modelZINB5, type = "pearson"))

#modelZINB5 best performing considering AIC, df & residual plots

#################### Model plotting

# ggplot theme

theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}

theme_set(theme_sleek())

library(wesanderson)
par(mfrow = c(1, 1))
pal <- wesanderson::wes_palette("Zissou1", 5, type = "continuous")
pal
cols<-c(pal[1], pal[2], pal[3], pal[4], pal[5])

#################### Effect sizes

# Extract coefficients and standard errors from ZERO model summary
coefs = as.data.frame(summary(modelZINB5)$coefficients$zero[,1:2])
names(coefs)[2] = "se" 
coefs$vars = rownames(coefs)
coefs$Estimate
coefs<-coefs[!coefs$vars %in% c('(Intercept)', 'Log(theta)'),]
coefs

# Coefficient plot
g1 <- ggplot(coefs, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour=cols[2]) +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se), 
                lwd=1, colour=cols[4], width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour=cols[1], width=0) +
  geom_point(size=4, pch=21, fill=cols[3]) +
  theme(axis.text.x = element_text(angle=1)) + 
  coord_flip() + 
  scale_y_reverse() + 
  labs(y = 'Coefficient estimate', x = '') +
  scale_x_discrete(limits= c('macroalgae.sc', 'complexity.sc', 'macroalgae.sc:complexity.sc', 'habitatpat.sc'), 
                 labels = c('Macroalgae', 'Complexity', 'Macroalgae * complexity', 'Patch reef type'))

# Extract coefficients and standard errors from COUNT model summary
coefs = as.data.frame(summary(modelZINB5)$coefficients$count[,1:2])
names(coefs)[2] = "se" 
coefs$vars = rownames(coefs)
coefs$Estimate


coefs<-coefs[!coefs$vars %in% c('(Intercept)', 'Log(theta)'),]

# Coefficient plot
g2 <- ggplot(coefs, aes(vars, Estimate)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour=cols[2]) +
  geom_errorbar(aes(ymin=Estimate - 1.96*se, ymax=Estimate + 1.96*se), 
                lwd=1, colour=cols[4], width=0) +
  geom_errorbar(aes(ymin=Estimate - se, ymax=Estimate + se), 
                lwd=2.5, colour=cols[1], width=0) +
  geom_point(size=4, pch=21, fill=cols[3]) +
  theme(axis.text.x = element_text(angle=0)) + coord_flip() +
  labs(y = 'Coefficient estimate', x = '') +
  scale_x_discrete(limits= c('macroalgae.sc:complexity.sc', 'complexity.sc', 'macroalgae.sc', 'habitatpat.sc'), 
                   labels = c('Macroalgae * complexity', 'Complexity', 'Macroalgae', 'Patch reef type'))

# Figure 1
cowplot::plot_grid(g1, g2, nrow=1, align = 'v', rel_heights=c(1, 0.5))

#################### Predictor trends

# Setup master prediction dataframe. All effects = 0
master<-data.frame(macroalgae.sc = rep(0,30),
                   complexity.sc = rep(0,30),
                   habitatpat.sc = rep(0,30),
                   rubble.sc = rep(0,30),
                   cca.sc = rep(0,30))

source("predict.zeroinfl.R")


# Macroalgae as focal covariate
pred1<-master
pred1$macroalgae.sc<-seq(min(urchins2017$macroalgae.sc),max(urchins2017$macroalgae.sc), length.out=30)
pred1$macroalgae<-seq(min(urchins2017$macroalgae),max(urchins2017$macroalgae), length.out=30)
pred1$estimate<-predict(modelZINB5, newdata=pred1, type='response')
cis<-predict.zeroinfl(modelZINB5, newdata = pred1, MC = 1000, se = TRUE, type='response')[[2]]
pred1$upper<-cis$upper
pred1$lower<-cis$lower

# Complexity 
pred2<-master
pred2$complexity.sc<-seq(min(urchins2017$complexity.sc),max(urchins2017$complexity.sc), length.out=30)
pred2$complexity<-seq(min(urchins2017$complexity),max(urchins2017$complexity), length.out=30)
pred2$estimate<-predict(modelZINB5, newdata=pred2, type='response')
cis<-predict.zeroinfl(modelZINB5, newdata = pred2, MC = 1000, se = TRUE, type='response')[[2]]
pred2$upper<-cis$upper
pred2$lower<-cis$lower

# Habitat Patch
pred3 <- master
pred3$habitatpat.sc<-seq(min(urchins2017$habitatpat.sc),max(urchins2017$habitatpat.sc), length.out=30)
pred3$estimate<-predict(modelZINB5, newdata=pred3, type='response')
cis<-predict.zeroinfl(modelZINB5, newdata = pred3, MC = 1000, se = TRUE, type='response')[[2]]
pred3$upper<-cis$upper
pred3$lower<-cis$lower


# Final figures of focal covariates with all other effects held to zero

p1 <- ggplot(pred1, aes(x = macroalgae, y = estimate)) +
  geom_line(colour=cols[1]) + labs(x = 'Macroalgae cover in %', y = expression(paste("Predicted sea urchin abundance")))
p1 <- p1 + geom_ribbon(data = pred1, aes(ymin = lower, ymax = upper), fill=cols[1],
                       alpha = 0.1)

p2 <- ggplot(pred2, aes(x = complexity, y = estimate)) +
  geom_line(colour=cols[1]) + labs(x = 'Structural complexity', y = '')
p2 <- p2 + geom_ribbon(data = pred2, aes(ymin = lower, ymax = upper), fill=cols[1],
                       alpha = 0.1)

p3 <- ggplot(pred3, aes(x = habitatpat.sc, y = estimate)) +
  geom_line(colour=cols[1]) + labs(x = expression(paste("Granitic & carbonate vs. patch habitat")), y = '')
p3 <- p3 + geom_ribbon(data = pred3, aes(ymin = lower, ymax = upper), fill=cols[1],
                       alpha = 0.1)

# Figure 2A

fig3 <- list(p1, p2, p3)
marrangeGrob(fig3, nrow=1, ncol=3, top='')

# Macroalgae * complexity interactions
pred6<-master
pred6$macroalgae.sc<-seq(min(urchins2017$macroalgae.sc),max(urchins2017$macroalgae.sc), length.out=30)
pred6$macroalgae<-seq(min(urchins2017$macroalgae),max(urchins2017$macroalgae), length.out=30)
pred6<-do.call('rbind', replicate(3, pred6, simplify=FALSE))
pred6$complexity.sc<-rep(c(0, 2, 3.5), each = 30)
pred6$complexity<-rep(c(0, 2, 3.5), each = 30)

pred6$estimate<-predict(modelZINB5, newdata=pred6, type='response')
cis<-predict.zeroinfl(modelZINB5, newdata = pred6, MC = 1000, se = TRUE, type='response')[[2]]
pred6$upper<-cis$upper	
pred6$lower<-cis$lower

# Figure 2B

p6 <- ggplot(pred6, aes(macroalgae, estimate)) + geom_line(colour=cols[1]) +
  geom_ribbon(aes(ymin = lower, ymax=upper), fill=cols[1], alpha=0.1) + 
  labs(x = 'Macroalgae cover in %', y = expression(paste("Predicted sea urchin abundance"))) + 
  facet_wrap(~complexity)
p6

#################### 2018 penning study

urchins2018 <- read.csv("Urchins2018.csv")

library(lattice)
library(lme4)
library(emmeans)

par(mfrow = c(1, 1))
xyplot(urchins2018$macroalgae ~ urchins2018$time | urchins2018$plot.id)
boxplot(macroalgae ~ time : plot.id, data = urchins2018)
dotchart(urchins2018$macroalgae)
hist(urchins2018$macroalgae)

modelurchins2018 <- glmer(macroalgae ~ time * plot.id + (1|name), data = urchins2018, family = Gamma)

plot(modelurchins2018)

summary(modelurchins2018)
emmeans(modelurchins2018, pairwise ~ plot.id * time)

# Figure 3

p7<-ggplot(urchins2018, aes(x=time, y=macroalgae)) +
  geom_boxplot(fill = cols[1]) + facet_grid(~plot.id) + 
  labs(x = '', y = 'Macroalgae cover (%)')
p7