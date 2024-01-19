library(tidyverse)
library(ggplot2)
library(LaplacesDemon)

n.occs <- 2 #number of occasions for occupancy data collection
n.states <- 3 #number of states

gamma.0s <- c(0.2, 0.5, 0.8) #intrinsic invasion probability
gamma.1s <- c(0.2, 0.5, 0.8) #effect of site characteristics
gamma.2s <- c(0.2, 0.5, 0.8) #effect of neighboring invasion state

eps.l0s <- c(0.9, 0.6, 0.3) #base eradication at low state
eps.l1s <- c(0.9, 0.6, 0.3)  #effect of eradication at low state
eps.h0s <- c(0.6, 0.4, 0.2)  #base eradication at high state
eps.h1s <- c(0.6, 0.4, 0.2) #effect of eradication at high state

phi0.lhs <- c(0.2, 0.6, 0.8) #transition from low to high
phi0.hhs <- c(0.2, 0.6, 0.8) #transition from high to high

phi1.lhs <- c(0.6, 0.4, 0.2) #transition from low to high
phi1.hhs <- c(0.6, 0.4, 0.2) #transition from high to high

TPM.48s <- array(NA, dim = c(n.states, n.states, 3))
TPM.48s[1,,1] <- c(1,0,0)
TPM.48s[2,,1] <- c(0.2,0.8,0)
TPM.48s[3,,1] <- c(0.1,0.2,0.7)

TPM.48s[1,,2] <- c(1,0,0)
TPM.48s[2,,2] <- c(0,1,0)
TPM.48s[3,,2] <- c(0,0,1)

TPM.48s[1,,3] <- c(0.8,0.2,0)
TPM.48s[2,,3] <- c(0,0.8,0.2)
TPM.48s[3,,3] <- c(0,0,1)

p.l0s <- c(0.2, 0.4, 0.6) #base detection for low state
p.l1s <- c(0.2, 0.4, 0.6) #effect of effort
p.h0s <- c(0.4, 0.6, 0.8) #base detection for high state
p.h1s <- c(0.4, 0.6, 0.8) #effect of effort
alpha.ls <- c(0.1, 0.2, 0.3) #difference in baseline detection between dat D and M
alpha.hs <- c(0.1, 0.2, 0.3)  #difference in baseline detection between dat D and M
deltas <- c(0.2, 0.6, 0.8)  # Probability of observing the high state given the species

# has been detected and the true state is high
search.hourss <- c(0.5, 1.1, 2) #effort is fixed
n.resource <- 40 #total resources we can use each week (hours)

#### Data ####
#---Habitat data---#
# effect of habitat quality on occupancy
n.sites <- 40
n.states <- 3

set.seed(03222021)
site.char <- runif(n.sites)

#---Initial state data---#
#Code that generated initial true state
State.init <- rep(NA, n.sites)
rate.init <- rep(NA, n.sites)
occ.init <- rep(NA, n.sites)
init.matrix <- array(NA, c(n.sites, n.states))
for(i in 1:n.sites){
  rate.init[i] <- mean(rbinom(100000,1,invlogit(gamma.0 + gamma.1*site.char[i]))) #invasion rate
  p.high <- 0.5 #say the probability of being in high state is 0.5

  occ.init[i] <- round(mean(rbern(100000,rate.init[i]))) #being invaded or not

  init.matrix[i,1] <- (1-rate.init[i])*occ.init[i] + (1-occ.init[i]) #empty
  init.matrix[i,2] <- (rate.init[i])*occ.init[i]*(1-p.high) #low state
  init.matrix[i,3] <- (rate.init[i])*occ.init[i]*(p.high) #high state

  State.init[i] <- rcat(1,init.matrix[i,1:3])
}

#### Empty jags arrays ####
S.init <- array(NA, c(n.sites,n.years, n.sims))
D.init <- array(NA, c(n.sites,n.years, n.sims))

#priors:
eps.l0.a <- array(NA, c(n.years, n.sims))
eps.l0.b <- array(NA, c(n.years, n.sims))
eps.l1.mean <- array(NA, c(n.years, n.sims))
eps.l1.sd <- array(NA, c(n.years, n.sims))

eps.h0.a <- array(NA, c(n.years, n.sims))
eps.h0.b <- array(NA, c(n.years, n.sims))
eps.h1.mean <- array(NA, c(n.years, n.sims))
eps.h1.sd <- array(NA, c(n.years, n.sims))

gamma.0.mean <- array(NA, c(n.years, n.sims))
gamma.0.sd <- array(NA, c(n.years, n.sims))
gamma.1.mean <- array(NA, c(n.years, n.sims))
gamma.1.sd <- array(NA, c(n.years, n.sims))
gamma.2.mean <- array(NA, c(n.years, n.sims))
gamma.2.sd <- array(NA, c(n.years, n.sims))

phi.lh.a <- array(NA, c(n.years, n.sims))
phi.lh.b <- array(NA, c(n.years, n.sims))
phi.lh1.mean <- array(NA, c(n.years, n.sims))
phi.lh1.sd <- array(NA, c(n.years, n.sims))

phi.hh.a <- array(NA, c(n.years, n.sims))
phi.hh.b <- array(NA, c(n.years, n.sims))
phi.hh1.mean <- array(NA, c(n.years, n.sims))
phi.hh1.sd <- array(NA, c(n.years, n.sims))

p.l0.a <- array(NA, c(n.years, n.sims))
p.l0.b <- array(NA, c(n.years, n.sims))
p.l1.mean <- array(NA, c(n.years, n.sims))
p.l1.sd <- array(NA, c(n.years, n.sims))
l.mean <- array(NA, c(n.years, n.sims))
l.sd <- array(NA, c(n.years, n.sims))
p.h0.a <- array(NA, c(n.years, n.sims))
p.h0.b <- array(NA, c(n.years, n.sims))
p.h1.mean <- array(NA, c(n.years, n.sims))
p.h1.sd <- array(NA, c(n.years, n.sims))
h.mean <- array(NA, c(n.years, n.sims))
h.sd <- array(NA, c(n.years, n.sims))

x <- list()
rhat_vals <- array(NA, c(n.years, n.sims))
sites <- list()
my.data <- list()
outs <- rep(NA,n.sims)
outputsfull <- rep(NA, n.sims)
outputs <- rep(NA, n.sims)
mcmcs <- rep(NA, n.sims)

alpha.eps.l0 <- rep(NA, n.sims)
beta.eps.l0 <- rep(NA, n.sims)
alpha.eps.h0 <- rep(NA, n.sims)
beta.eps.h0 <- rep(NA, n.sims)
alpha.phi.lh <- rep(NA, n.sims)
beta.phi.lh <- rep(NA, n.sims)
alpha.phi.hh <- rep(NA, n.sims)
beta.phi.hh <- rep(NA, n.sims)
alpha.p.l0 <- rep(NA, n.sims)
beta.p.l0 <- rep(NA, n.sims)
alpha.p.h0 <- rep(NA, n.sims)
beta.p.h0 <- rep(NA, n.sims)

State.est <- rep(NA, n.sims)
eps.l0.est <- rep(NA, n.sims)
eps.l1.est <- rep(NA, n.sims)
eps.h0.est <- rep(NA, n.sims)
eps.h1.est <- rep(NA, n.sims)
gamma.0.est <- rep(NA, n.sims)
gamma.1.est <- rep(NA, n.sims)
gamma.2.est <- rep(NA, n.sims)
phi0.lh.est <- rep(NA, n.sims)
phi1.lh.est <- rep(NA, n.sims)
phi0.hh.est <- rep(NA, n.sims)
phi1.hh.est <- rep(NA, n.sims)
p.l0.est <- rep(NA, n.sims)
p.l1.est <- rep(NA, n.sims)
alpha.l.est <- rep(NA, n.sims)
p.h0.est <- rep(NA, n.sims)
p.h1.est <- rep(NA, n.sims)
alpha.h.est <- rep(NA, n.sims)

all.State.est <- rep(NA, n.sims)
all.eps.l0.est <- rep(NA, n.sims)
all.eps.l1.est <- rep(NA, n.sims)
all.eps.h0.est <- rep(NA, n.sims)
all.eps.h1.est <- rep(NA, n.sims)
all.gamma.0.est <- rep(NA, n.sims)
all.gamma.1.est <- rep(NA, n.sims)
all.gamma.2.est <- rep(NA, n.sims)
all.phi0.lh.est <- rep(NA, n.sims)
all.phi1.lh.est <- rep(NA, n.sims)
all.phi0.hh.est <- rep(NA, n.sims)
all.phi1.hh.est <- rep(NA, n.sims)
all.p.l0.est <- rep(NA, n.sims)
all.p.l1.est <- rep(NA, n.sims)
all.alpha.l.est <- rep(NA, n.sims)
all.p.h0.est <- rep(NA, n.sims)
all.p.h1.est <- rep(NA, n.sims)
all.alpha.h.est <- rep(NA, n.sims)

initial.values <- list()

TPM.est <- array(NA, c(n.states, n.sites,n.sims,n.states))
D.est <- array(NA, c(n.sites,n.sims))
gamma.est <- array(NA, c(n.sites,n.sims))
eps.l.est <- array(NA, c(n.sites,n.sims))
eps.h.est <- array(NA, c(n.sites,n.sims))
phi.lh.est <- array(NA, c(n.sites,n.sims))
phi.hh.est <- array(NA, c(n.sites,n.sims))

prev.state <- array(NA, c(n.sites, n.sims))
States.mean.round <- array(NA, c(n.sites, n.years, n.sims))
States.mean <- array(NA, c(n.sites, n.years, n.sims))
S.end <- array(NA, c(n.sites, n.sims))

#### Save workspace ####
save.image(file = "parameters.RData")

#load("parameters.RData")

#### DETECTION CURVES ####
search.hoursvals <- c(search.hourss,5) 
#search.hoursvals <- c(0.25, 0.5, 1, 2, 4, 8, 16) 
n <- length(search.hoursvals)
logsearch.efforts <- log(search.hoursvals) #log search effort

pM.l.lowd <- invlogit(p.l0s[1] + p.l1s[1]*logsearch.efforts + alpha.ls[1])
pM.l.highd <- invlogit(p.l0s[3] + p.l1s[3]*logsearch.efforts + alpha.ls[3])

pM.h.lowd <- invlogit(p.h0s[1] + p.h1s[1]*logsearch.efforts + alpha.hs[1])
pM.h.highd <- invlogit(p.h0s[3] + p.h1s[3]*logsearch.efforts + alpha.hs[3])

det.df <- data.frame(detection = c(pM.l.lowd, pM.l.highd, pM.h.lowd, pM.h.highd),
           val = c(rep(1,n), rep(2,n), rep(1,n), rep(2,n)),
           state = c(rep(1,n), rep(1,n), rep(2,n), rep(2,n)),
           effort = rep(search.hoursvals,4),
           color = c(rep('low state, low detection',n),rep('low state, high detection',n),
                         rep('high state, low detection',n), rep('high state, high detection',n)),
           colorvals = c(rep('grey',n),rep('black',n),
                     rep('red1',n), rep('red4',n))
           )

colorvals <- c(rep('lightblue',n), rep('darkblue',n),
              rep('lightgreen',n),rep('darkgreen',n)
              )

ggplot(det.df)+
  geom_point(aes(x = effort, y = detection, group = val, color = color), color = colorvals)+
  geom_line(aes(x = effort, y = detection, group = interaction(val, state), color = color), color = colorvals)+
  xlab("Search effort (hours)") +
  ylab("Detection probability")+
  theme_bw() + theme(#panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.x = element_text(size = 16),
                     axis.text.x = element_text(size = 14),
                     axis.title.y = element_text(size = 16),
                     axis.line = element_line(colour = "black"))

