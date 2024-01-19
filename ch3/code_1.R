library(LaplacesDemon)
library(rjags)
library(R2jags)
library(MCMCvis)
library(tidyverse)
library(strex)
library(plyr)


#Remove at downstream 5 locations but collect monitoring data
#at next 5 downstream locations

start.time <- Sys.time()

#### JAGS model ####
sink("Flower_multistate_datM.txt")
cat("
model{

# -------------------------------------------------
# Parameters:
# gamma: invasion probability
# eps: erradication probability
# p: probability of observing the plant
# -------------------------------------------------
# States (S):
# 1 empty
# 2 low abundance
# 3 high abundance 
# 
# Observations Dat M:  
# 1 not detected
# 2 low abundance
# 3 high abundance
#
# -------------------------------------------------

#### PRIORS ####
  #Erradication:
  eps.l0 ~ dbeta(eps.l0.a,eps.l0.b)T(0.0001,0.9999) #eradication when at low state
  eps.l1 ~ dnorm(eps.l1.mean, eps.l1.tau) #effect of eradication 
  eps.l1.tau <- 1/(eps.l1.sd * eps.l1.sd) #precision parameter
  
  eps.h0 ~ dbeta(eps.h0.a,eps.h0.b)T(0.0001,0.9999) #eradication when at high state
  eps.h1 ~ dnorm(eps.h1.mean, eps.h1.tau) #effect of eradication 
  eps.h1.tau <- 1/(eps.h1.sd * eps.h1.sd) #precision parameter
  
  #Invasion:
  gamma.0 ~dnorm(gamma.0.mean,gamma.0.tau) #intrinsic invasion probability
  gamma.0.tau <- 1/(gamma.0.sd*gamma.0.sd) #precision parameter
  gamma.1 ~dnorm(gamma.1.mean, gamma.1.tau) #effect of site characteristics on invasion probability
  gamma.1.tau <- 1/(gamma.1.sd*gamma.1.sd) #precision parameter
  gamma.2 ~dnorm(gamma.2.mean, gamma.2.tau) #effect of Neighboring invasion state
  gamma.2.tau <- 1/(gamma.2.sd*gamma.2.sd) #precision parameter
  
  #State transition:
  phi0.lh ~ dbeta(phi.lh.a, phi.lh.b)T(0.0001,0.9999) #transition from low to high
  phi1.lh ~ dnorm(phi.lh1.mean, phi.lh1.tau) #effect of removal on transition
  phi.lh1.tau <- 1/(phi.lh1.sd*phi.lh1.sd) #precision parameter
  
  phi0.hh ~ dbeta(phi.hh.a, phi.hh.b)T(0.0001,0.9999) #transition from high to high
  phi1.hh ~ dnorm(phi.hh1.mean, phi.hh1.tau) #effect of removal on transition
  phi.hh1.tau <- 1/(phi.hh1.sd*phi.hh1.sd) #precision parameter
  
  #Detection low state:
  p.l0 ~ dbeta(p.l0.a, p.l0.b)T(0.0001,0.9999) #base detection for low state
  p.l1 ~ dnorm(p.l1.mean, p.l1.tau) #effect of effort 
  p.l1.tau <- 1/(p.l1.sd * p.l1.sd) #precision parameter
  alpha.l ~ dnorm(l.mean,l.tau) #difference in baseline detection between dat D and M
  l.tau <- 1/(l.sd * l.sd) #precision
  
  #Detection high state:
  p.h0 ~ dbeta(p.h0.a, p.h0.b)T(0.0001,0.9999) #base detection for high state
  p.h1 ~ dnorm(p.h1.mean, p.h1.tau) #effect of effort 
  p.h1.tau <- 1/(p.h1.sd * p.h1.sd) #precision parameter
  alpha.h ~ dnorm(h.mean,h.tau) #difference in baseline detection between dat D and M
  h.tau <- 1/(h.sd * h.sd) #precision

  logit(pM.l) <- p.l0 + p.l1*logeffort + alpha.l #detection low state
  logit(pM.h) <- p.h0 + p.h1*logeffort + alpha.h #detection high state
  
  #initial occupancy probabilities
  psi[1:3] ~ ddirch(alpha) #alpha = rep(1,3)
  
#--------------------------------------------------#
# STATE TRANSITION
for (i in 1:n.sites){  
  # State transition probabilities (TPM): probability of S(t+1) given S(t)
  for (t in 1:n.weeks){
  
    logit(gamma[i,t]) <-gamma.0 + gamma.1*site.char[i] + gamma.2*D[i,t] #invasion probability
    logit(eps.l[i,t]) <- eps.l0 + eps.l1*rem.vec[i,t]*removal.hours[2] #erradication low state
    logit(eps.h[i,t]) <- eps.h0 + eps.h1*rem.vec[i,t]*removal.hours[3] #erradication high state
                                        # rem.vec[i] = 0,1 if 0, then no removal and no erradiction
    
    logit(phi.lh[i,t]) <- phi0.lh - phi1.lh*rem.vec[i,t]*removal.hours[2]
    logit(phi.hh[i,t]) <- phi0.hh - phi1.hh*rem.vec[i,t]*removal.hours[3]
    
    #index = [current state, location, time, future state]
    #empty stay empty
    TPM[1,i,t,1] <- 1-gamma[i,t] #1-gamma = not invasion probability
    
    #empty to low abundance
    TPM[1,i,t,2] <- gamma[i,t] #invasion probability
    
    #empty to high abundance
    TPM[1,i,t,3] <- 0 #invasion probability

    #low abundance to empty
    TPM[2,i,t,1] <- eps.l[i,t] #erradication probability
                                      
    #low abundance to low abundance
    TPM[2,i,t,2] <- (1- eps.l[i,t])*(1-phi.lh[i,t]) #erradication failure probability
    
    #low abundance to high abundance
    TPM[2,i,t,3] <- (1- eps.l[i,t])*(phi.lh[i,t])
    
    #high abundance to empty
    TPM[3,i,t,1] <- eps.h[i,t] #erradication probability
    
    #high abundance to low abundance
    TPM[3,i,t,2] <- (1- eps.h[i,t])*(1-phi.hh[i,t]) #erradication failure probability
    
    #high abundance to high abundance
    TPM[3,i,t,3] <- (1- eps.h[i,t])*(phi.hh[i,t])
    
    #--------------------------------------------------#
    # OBSERVATION PROBABILITIES (for multi state detection/nondetection data)
    
    for(j in 1:n.occs){

      #Empty and not observed  
      P.datM[1,i,j,t,1] <- 1
      
      #Empty and observed low 
      P.datM[1,i,j,t,2] <- 0
      
      #Empty and observed high
      P.datM[1,i,j,t,3] <- 0
   
      #Low state and not observed
      P.datM[2,i,j,t,1] <- 1-pM.l #not detected probability low state
      
      #Low state and observed low
      P.datM[2,i,j,t,2] <- pM.l #detection probability low state
      
      #Low state and observed high
      P.datM[2,i,j,t,3] <- 0 #detection probability low state
      
      #High state and not observed
      P.datM[3,i,j,t,1] <- 1-pM.h #not detected probability high state
      
      #High state and observed low
      P.datM[3,i,j,t,2] <- 0
      
      #High state and observed high
      P.datM[3,i,j,t,3] <- pM.h
      
    } #j 
  } #t
} #i

  #### Likelihood ####
  for (i in 1:n.sites){
      
    #-- Initial State: --# 
    
      State[i,1] ~ dcat(psi) #psi is written above in the priors
    
      D[i,1] <- sum(State[neighbors[i,], 1])/n.neighbors[i] #state of neighbors 
    
    #-- State Model: --#
    for (t in 2:n.weeks){ 
      # State process: state given previous state and transition probability
      State[i,t] ~ dcat(TPM[State[i,t-1], i, t-1, ]) 
      
      #below is not correct because it does not lead to a true state: 
      D[i,t] <- sum(State[neighbors[i,], t])/n.neighbors[i] #state of neighbors 
       
     
    } #t loop

    #----- Observation Model -----#
    for(j in 1:n.occs){
      for(t in 1:n.weeks){
        # Observation process: draw observation given current state
        yM[i,j,t] ~ dcat(P.datM[State[i,t], i, j, t,]) 
        
      } #t
    } #j

    #Derived parameter: final estimated state
    State.fin[i] <- State[i,n.weeks] #state after 4 weeks
    
  } #i


} #end model
", fill = TRUE)
sink()


#### Data and parameters ####
load("parameters.RData")

n.sims <-  2 #number of simulations
n.sites <- 40 #number of sites
n.years <- 10 #number of years
n.weeks <- 4 #number of weeks
n.occs <- 2 #number of occasions for occupancy data collection
n.states <- 3 #number of states

##### STATE VALUES ####
gamma.0 <- gamma.0s[1] #intrinsic invasion probability
gamma.1 <- gamma.1s[1] #effect of site characteristics
gamma.2 <- gamma.2s[1] #effect of neighboring invasion state

eps.l0 <- eps.l0s[1] #base eradication at low state
eps.l1 <- eps.l1s[1] #effect of eradication at low state
eps.h0 <- eps.h0s[1] #base eradication at high state
eps.h1 <- eps.h1s[1] #effect of eradication at high state
phi0.lh <- phi0.lhs[1] #transition from low to high base
phi1.lh <- phi1.lhs[1] #effect of removal on transition to low and high
phi0.hh <- phi0.hhs[1] #transition from high to high base
phi1.hh <- phi1.hhs[1] #effect of removal on transition to high and high

##### OBSERVATION VALUES ####
p.l0 <- p.l0s[1] #base detection for low state
p.l1 <- p.l1s[1] #effect of effort
alpha.l <- alpha.ls[1] #difference in baseline detection between dat D and M

p.h0 <- p.h0s[1] #base detection for high state
p.h1 <- p.h1s[1] #effect of effort
alpha.h <- alpha.hs[1] #difference in baseline detection between dat D and M

TPM.48 <- TPM.48s[,,1] #TPM matrix for 48 week period

search.hours <- search.hourss[1] #search effort

removal.hours <- c(0, 2, 3) #it removal takes 2 hours if in low state and 3 hours if in high state
n.resource <- 40 #total hours per week

#---- arrays ----#
gamma <- array(NA, c(n.sites, n.weeks, n.years, n.sims))
eps.l <- array(NA, c(n.sites, n.weeks, n.years, n.sims))
eps.h <- array(NA, c(n.sites, n.weeks, n.years, n.sims))
phi.lh <- array(NA, c(n.sites, n.weeks, n.years, n.sims))
phi.hh<- array(NA, c(n.sites, n.weeks, n.years, n.sims))

TPM<- array(NA, c(n.states,n.sites,n.weeks, n.years + 1,n.sims, n.states)) 

#---Habitat data---#
# effect of habitat quality on occupancy
site.char <- site.char
State.init <- State.init
State <- array(NA,c(n.sites, n.weeks, n.years+1, n.sims)) #state array

#---Neighbor data---#
D <- array(NA, c(n.sites, n.weeks, n.years+1, n.sims)) #neighbors array
num.neighbors <- 2 #one upstream, one downstream
neighbors <- matrix(NA, nrow = n.sites, ncol = num.neighbors) #neighbors matrix, each row (site) identifies the neighbors for that site 
neighbors[1,1] <- 2 #site 1 only has site 2 as its neighbor 
neighbors[2:n.sites, 1] <- seq(1,n.sites-1) #filling in upstream neighbors
neighbors[n.sites,2] <- n.sites-1 #end site only has end site -1 as its neighbor
neighbors[1:(n.sites-1), 2] <- seq(2,n.sites) #filling in downstream neighbors
n.neighbors <- rep(2,n.sites)
n.neighbors[1] <- n.neighbors[n.sites] <- 1

#--- removal data and occupancy data ---#

sites.rem.M <- array(NA, c(n.sites, n.weeks, n.years, n.sims)) 

for(s in 1: n.sims){
  sites.rem.M[,1,1,s] <- sample(n.sites, n.sites, replace = F)
}


yM <- array(NA, c(n.sites, n.occs, n.weeks, n.years, n.sims)) 
resource.total <- array(0, c(n.weeks, n.years, n.sims)) 

logsearch.effort <- log(search.hours) #log search effort

pM.l <- invlogit(p.l0 + p.l1*logsearch.effort + alpha.l) #low state detection probability (base detection + effect of effort)

pM.h <- invlogit(p.h0 + p.h1*logsearch.effort + alpha.h) #high state detection probability (base detection + effect of effort)


P.datM <- array(NA, dim = c(n.states, n.states))
P.datM[1,] <- c(1,0,0)
P.datM[2,] <- c(1-pM.l, pM.l, 0)
P.datM[3,] <- c(1-pM.h, 0, pM.h)
  
rem.vec <- array(NA, c(n.sites, n.weeks, n.years, n.sims)) #removal sites array

#### JAGS arrays ####
#empty arrays are loaded with parameter data

####################################################################################
#### Run Adaptive Management ####

year <- 1
#n.years <- 2

for(year in 1:n.years){
  #--------------------------------------------------------------------------------#
  #### 1. Simulate the truth ####
  
  ### Steps: 
  #---1. Simulate the truth
  #---2. Simulate occupancy data collection (include removal data)
  #---3. Go back to step 1 and take into account removal that previously occurred into state process
  
  ##### Week 1 State model only #####
  week <- 1
  ###### Week 1 year 1 #####
  if(year == 1){

    State[,1,year,1:n.sims] <- State.init #first week state is from data
    
    for(s in 1:n.sims){
      for(i in 1:n.sites){
        D[i,1,1,s] <- sum(State[neighbors[i,], 1,1,s])/n.neighbors[i] #state of neighbors
      }
    }
    
    gamma[,1,1,] <-invlogit(gamma.0 + gamma.1*site.char + gamma.2*D[,1,1,]) #invasion (week 1 year 1)
    eps.l[,1,1,] <- invlogit(eps.l0) #eradication low (week 1 year 1)
    eps.h[,1,1,] <- invlogit(eps.h0) #eradication high (week 1 year 1)
    phi.lh[,1,1,] <- invlogit(phi0.lh) #transiotion low to high
    phi.hh[,1,1,] <- invlogit(phi1.lh) #transition high to high
    
    # TPM used for week 2
    TPM[1,1:n.sites,1,1,,1] <- 1-gamma[,1,1,] #empty to empty (week 1 year 1)
    TPM[1,1:n.sites,1,1,,2] <- gamma[,1,1,] #empty to low (week 1 year 1)
    TPM[1,1:n.sites,1,1,,3] <- 0 #empty to high (week 1 year 1)
    
    TPM[2,1:n.sites,1,1,,1] <- eps.l[,1,1,] #low to empty (week 1 year 1)
    TPM[2,1:n.sites,1,1,,2] <- (1- eps.l[,1,1,])*(1-phi.lh[,1,1,]) #low to low (week 1 year 1)
    TPM[2,1:n.sites,1,1,,3] <- (1- eps.l[,1,1,])*(phi.lh[,1,1,]) #low to high (week 1 year 1)
    
    TPM[3,1:n.sites,1,1,,1] <- eps.h[,1,1,] #high to empty (week 1 year 1)
    TPM[3,1:n.sites,1,1,,2] <- (1- eps.h[,1,1,])*(1-phi.hh[,1,1,]) #high to low (week 1 year 1)
    TPM[3,1:n.sites,1,1,,3] <- (1- eps.h[,1,1,])*(phi.hh[,1,1,]) #high to high (week 1 year 1)
    
  } #ends year = 1 loop
  
  ###### Week 1 year >1 #####
  #for all years > 1 we need to project 48 weeks forward
  
  if(year > 1){
    for(s in 1:n.sims){

      
      for(i in 1:n.sites){ #State process: state given previous state and transition probability
        State[i,week,year,s] <- rcat(1,TPM.48[State[i,4,(year-1),s], ]) 
      }
        
      for(i in 1:n.sites){
        D[i,week,year,s] <- sum(State[neighbors[i,], week,year,s])/n.neighbors[i] #state of neighbors
      }
      
    #--- Data for the TPM for the next week: week 2 ---#
    #prev.rem.vec = vector of 0 and 1s indicating where removal previously occurred
    prev.rem.vec <- replace(rem.vec[,4,(year-1),s], is.na(rem.vec[,4,(year-1),s]), 0) 
    
    #invasion probability =  base invasion + effect of site habitat + effect of neighbor being invaded
    gamma[,week,year,s] <-invlogit(gamma.0 + gamma.1*site.char + gamma.2*D[,week,year,s]) 
    
    # eradication probability = base + effect of previous removal (removal*removal hours)
    eps.l[,week,year,s] <- invlogit(eps.l0 + eps.l1*prev.rem.vec*removal.hours[2]) #low eradication 
    eps.h[,week,year,s] <- invlogit(eps.h0 + eps.h1*prev.rem.vec*removal.hours[3]) #high eradication
    
    #transition rates
    phi.lh[,week,year,s] <- invlogit(phi0.lh - phi1.lh*prev.rem.vec*removal.hours[2])
    phi.hh[,week,year,s] <- invlogit(phi0.hh - phi1.hh*prev.rem.vec*removal.hours[3])
    
    
    TPM[1,1:n.sites,week,year,s,1] <- 1-gamma[,week,year,s] #empty to empty
    TPM[1,1:n.sites,week,year,s,2] <- gamma[,week,year,s] #empty to low
    TPM[1,1:n.sites,week,year,s,3] <- 0 #empty to high
    
    TPM[2,1:n.sites,week,year,s,1] <- eps.l[,week,year,s] #low to empty (eradication)
    TPM[2,1:n.sites,week,year,s,2] <- (1- eps.l[,week,year,s])*(1-phi.lh[,week,year,s]) #low to low (eradication failure)
    TPM[2,1:n.sites,week,year,s,3] <- (1- eps.l[,week,year,s])*(phi.lh[,week,year,s]) #low to high 
    
    TPM[3,1:n.sites,week,year,s,1] <- eps.h[,week,year,s] #high to empty (eradication)
    TPM[3,1:n.sites,week,year,s,2] <- (1- eps.h[,week,year,s])*(1-phi.hh[,week,year,s]) #high to low 
    TPM[3,1:n.sites,week,year,s,3] <- (1- eps.h[,week,year,s])*(phi.hh[,week,year,s]) #high to high
    
    } #ends s loop
  } #ends year > 1 loop
  
  ##### Week 1+ State and Observation model #####
  for(s in 1:n.sims){
    for(week in 1:n.weeks){
      
      ###### State process ######
      
      if(week > 1){
        for(i in 1:n.sites){ #State process: state given previous state and transition probability
          State[i,week,year,s] <- rcat(1,TPM[State[i,week-1,year,s], i, week-1, year, s,]) 
        }
        
        for(i in 1:n.sites){ #state of neighbors
          D[i,week,year,s] <- sum(State[neighbors[i,], week,year,s])/n.neighbors[i] #state of neighbors
        }
        
        #--- Data for the TPM for the next week: week 2+ ---#
        #prev.rem.vec = vector of 0 and 1s indicating where removal previously occurred
        prev.rem.vec <- replace(rem.vec[,week-1,year,s], is.na(rem.vec[,week-1,year,s]), 0)
        
        #invasion probability =  base invasion + effect of site habitat + effect of neighbor being invaded
        gamma[,week,year,s] <-invlogit(gamma.0 + gamma.1*site.char + gamma.2*D[,week,year,s]) 
        
        # eradication probability = base + effect of previous removal (removal*removal hours)
        eps.l[,week,year,s] <- invlogit(eps.l0 + eps.l1*prev.rem.vec*removal.hours[2]) #low eradication
        eps.h[,week,year,s] <- invlogit(eps.h0 + eps.h1*prev.rem.vec*removal.hours[3]) #high eradication
        
        #transition rates
        phi.lh[,week,year,s] <- invlogit(phi0.lh - phi1.lh*prev.rem.vec*removal.hours[2])
        phi.hh[,week,year,s] <- invlogit(phi0.hh - phi1.hh*prev.rem.vec*removal.hours[3])
        
        TPM[1,1:n.sites,week,year,s,1] <- 1-gamma[,week,year,s] #empty to empty
        TPM[1,1:n.sites,week,year,s,2] <- gamma[,week,year,s] #empty to low
        TPM[1,1:n.sites,week,year,s,3] <- 0 #empty to high
        
        TPM[2,1:n.sites,week,year,s,1] <- eps.l[,week,year,s] #low to empty (eradication)
        TPM[2,1:n.sites,week,year,s,2] <- (1- eps.l[,week,year,s])*(1-phi.lh[,week,year,s]) #low to low (eradication failure)
        TPM[2,1:n.sites,week,year,s,3] <- (1- eps.l[,week,year,s])*(phi.lh[,week,year,s]) #low to high
        
        TPM[3,1:n.sites,week,year,s,1] <- eps.h[,week,year,s] #high to empty (eradication)
        TPM[3,1:n.sites,week,year,s,2] <- (1- eps.h[,week,year,s])*(1-phi.hh[,week,year,s]) #high to low
        TPM[3,1:n.sites,week,year,s,3] <- (1- eps.h[,week,year,s])*(phi.hh[,week,year,s]) #high to high
        
        #week: Identify the sites where removal will occur 
        n.pre.visit <- length(which(rem.vec[,week-1,year,s] >= 0)) #number of sites that were sampled last week
        #put last weeks sampling sites at the end of the sampling queue 
        sites.rem.M[,week,year,s] <- c(sites.rem.M[,(week-1),year,s][-1:-n.pre.visit],
                                     sites.rem.M[,(week-1),year,s][1:n.pre.visit])
      } #week > 1
      
      ##### Observation process #######
      # Observation process: draw observation given current state
      
      for(i in sites.rem.M[,week,year,s]){ #order of sites where removal occurs
        
        #A. while we still have resources to spend:
        if(resource.total[week,year,s] < n.resource){
          
          #1. first occasion occupancy data (1 = not detected, 2 = detected)
          yM[i,1,week, year, s] <- rcat(1, P.datM[State[i,week,year,s], ])
          
          #2. second occasion occupancy data
          #2a. if seen in first occasion, do not search again and remove the rush
          if(yM[i,1,week, year, s] > 1){ 
            yM[i,2, week,year, s] <- NA #no occupancy data because we did not need to search again
            rem.vec[i,week,year,s] <- 1 #notes that removal occurred that week at that site
            
            #Calculating resources used = resources already used + search hours + removal hours
            resource.total[week,year,s] <- resource.total[week,year,s] + search.hours + removal.hours[State[i,week,year,s]]
            
          }else{
            #2b. If not seen the first occasion, we need to search again:
            #Second occasion occupancy data
            yM[i,2, week, year, s] <- rcat(1, P.datM[State[i,week,year,s], ])
            
            #2bi. If seen at the second occasion:
            if(yM[i,2, week, year, s] > 1){ #if seen (state observed > 1) the second time
              rem.vec[i,week,year,s] <- 1 #notes that removal occurred that week at that site
              
              #Calculating resources used = resources already used + 2*search hours + removal hours
              resource.total[week,year,s] <- resource.total[week,year,s] + 2*search.hours + removal.hours[State[i,week,year,s]]
            } 
            
            #2bi. If we do not detect flowering rush during the second occasion:
            if(yM[i,2, week, year, s]==1){ #if not seen (state observed = 1)
              rem.vec[i,week,year,s] <- 0 #notes removal did not occur
              
              #Calculating resources used = resources already used + 2*search hours
              resource.total[week,year,s] <- resource.total[week,year,s] + 2*search.hours 
            } 
          }
        
        #B. if we do not have any more resources to spend:
        }else{
          yM[i,1:2, week, year, s] <- NA #no occupancy data
          rem.vec[i,week,year,s] <- NA #removal did not occur
        }
        
      } #ends site loop
    } #ends week loop
  } #ends sims loop  

  #--------------------------------------------------------------------------------#
  #### 2. Learning: ####
  
  ##### 2a. Update priors #####
  ###### year 1 priors #####
  #------------------------Year 1 Priors------------------------#
  if(year == 1){
    
    # --- eps.l ---  eradication low state -------------------- #
      #eps.l0 = base eradication at low state (beta distribution)
      eps.l0.a[year,] <- 1 #alpha shape
      eps.l0.b[year,] <- 1 #beta shape
      
      #eps.l1 = effect of eradication at low state (normal distribution)
      eps.l1.mean[year,] <- 0 #mean
      eps.l1.sd[year,] <-  10 #sd
    
    # --- eps.h ---  eradication high state ------------------- #
      #eps.h0 = base eradication at high state (beta distribution)
      eps.h0.a[year,] <- 1 #alpha shape
      eps.h0.b[year,] <- 1 #beta shape
    
      #eps.h1 = effect of eradication at high state (normal distribution)
      eps.h1.mean[year,] <- 0 #mean
      eps.h1.sd[year,] <- 10 #sd
    
    # --- gamma ---  invasion -------------------------------- #  
      #gamma.0 = intrinsic invasion (normal distribution)
      gamma.0.mean[year,] <- 0 #mean
      gamma.0.sd[year, ] <- 10 #sd
      
      #gamma.1 = effect of site characteristics (normal distribution)
      gamma.1.mean[year,] <- 0 #mean
      gamma.1.sd[year,] <- 10 #sd
    
      #gamma.2 = effect of neighboring state (normal distribution)
      gamma.2.mean[year,] <- 0 #mean
      gamma.2.sd[year,] <- 10 #sd
    
    
    # --- phi ---  transition rates -------------------------- #
      #phi.lh = base transition low to high (beta distribution)
      phi.lh.a[year,] <- 1 #alpha shape
      phi.lh.b[year,] <- 1 #beta shape
      
      #effect of removal on transition from low to high
      phi.lh1.mean[year,] <- 0 #mean
      phi.lh1.sd[year,] <- 10 #sd
      
      #phi.hh = transition high to high (beta distribution)
      phi.hh.a[year,] <- 1 #alpha shape
      phi.hh.b[year,] <- 1 #beta shape
      
      #effect of removal on transition from high to high
      phi.hh1.mean[year,] <- 0 #mean
      phi.hh1.sd[year,] <- 10 #sd
   
    # --- p.l ---  detection low state ----------------------- #
      #p.l.0 = base detection low state (beta distribution)
      p.l0.a[year,] <- 1 #alpha shape
      p.l0.b[year,] <- 1 #beta shape
      
      #p.l.1 = effect of effort (normal distribution)
      p.l1.mean[year,] <- 0 #mean
      p.l1.sd[year,] <- 10 #sd
    
      
    # --- alpha.l --- difference in baseline detection btwn dat D and M -- #  
      l.mean[year,] <- 0 #mean
      l.sd[year,] <- 1 #sd
      
    # --- p.h ---  detection high state ---------------------- #
      #p.h.0 = base detection high state (beta distribution)
      p.h0.a[year,] <- 1 #alpha shape
      p.h0.b[year,] <- 1 #beta shape
      
      #p.h.1 = effect of effort (normal distribution)
      p.h1.mean[year,] <- 0 #mean
      p.h1.sd[year,] <- 10 #sd
      
    # --- alpha.h --- difference in baseline detection btwn dat D and M -- #   
      h.mean[year,] <- 0 #mean
      h.sd[year,] <- 1 #sd
      
    ##### UNSURE ####
    # --- S.init and D.init ---  Initial states ------------ #
      alpha <- rep(1,n.states) #initial state probability vector
      
  } else{
    
    ###### year 1+ priors #####
    for(s in 1:n.sims){
    #------------------------Year 1+ Priors------------------------#
    # --- eps.l ---  eradication low state ----------------------- #
    #eps.l0 = base eradication at low state (beta distribution)
    alpha.eps.l0[s] <- paste("alpha.eps.l0", s, sep = "_")
    #assigning alpha values for beta: alpha = (1-mean)*(1+cv^2)/cv^2
    assign(alpha.eps.l0[s],
           (1 - get(eps.l0.est[s])$mean*(1 + get(eps.l0.est[s])$cv^2))/(get(eps.l0.est[s])$cv^2))
    
    #assigning beta values for beta: beta = (alpha)*(1-mean)/(mean)
    beta.eps.l0[s]<- paste("beta.eps.l0", s, sep = "_")
    
    assign(beta.eps.l0[s],
           get(alpha.eps.l0[s])*(1 - get(eps.l0.est[s])$mean)/get(eps.l0.est[s])$mean)
    
    eps.l0.a[year,s] <- get(alpha.eps.l0[s]) #alpha shape
    eps.l0.b[year,s] <- get(beta.eps.l0[s]) #beta shape
    
    #eps.l1 = effect of eradication at low state (normal distribution)
    eps.l1.mean[year,s] <- get(eps.l1.est[s])$mean #mean
    eps.l1.sd[year,s] <-  get(eps.l1.est[s])$sd #sd
    
    # --- eps.h ---  eradication high state ------------------------- #
    #eps.h0 = base eradication at high state (beta distribution)
    alpha.eps.h0[s] <- paste("alpha.eps.h0", s, sep = "_")
    #assigning alpha values for beta: alpha = (1-mean)*(1+cv^2)/cv^2
    assign(alpha.eps.h0[s],
           (1 - get(eps.h0.est[s])$mean*(1 + get(eps.h0.est[s])$cv^2))/(get(eps.h0.est[s])$cv^2))
    
    #assigning beta values for beta: beta = (alpha)*(1-mean)/(mean)
    beta.eps.h0[s]<- paste("beta.eps.h0", s, sep = "_")
    
    assign(beta.eps.h0[s],
           get(alpha.eps.h0[s])*(1 - get(eps.h0.est[s])$mean)/get(eps.h0.est[s])$mean)
    
    eps.h0.a[year,s] <- get(alpha.eps.h0[s]) #alpha shape
    eps.h0.b[year,s] <- get(beta.eps.h0[s]) #beta shape
    
    #eps.h1 = effect of eradication at high state (normal distribution)
    eps.h1.mean[year,s] <- get(eps.h1.est[s])$mean #mean
    eps.h1.sd[year,s] <-  get(eps.h1.est[s])$sd #sd
    
    # --- gamma ---  invasion -------------------------------------- #  
    #gamma.0 = intrinsic invasion (normal distribution)
    gamma.0.mean[year,s] <- get(gamma.0.est[s])$mean #mean
    gamma.0.sd[year,s ] <- get(gamma.0.est[s])$sd
    
    #gamma.1 = effect of site characteristics (normal distribution)
    gamma.1.mean[year,s] <- get(gamma.1.est[s])$mean #mean
    gamma.1.sd[year,s ] <- get(gamma.1.est[s])$sd
    
    #gamma.2 = effect of neighboring state (normal distribution)
    gamma.2.mean[year,s] <- get(gamma.2.est[s])$mean #mean
    gamma.2.sd[year,s ] <- get(gamma.2.est[s])$sd
    
    # --- phi ---  transition rates ------------------------------- #
    #phi.lh = transition low to high (beta distribution)
    alpha.phi.lh[s] <- paste("alpha.phi.lh", s, sep = "_")
    #assigning alpha values for beta: alpha = (1-mean)*(1+cv^2)/cv^2
    assign(alpha.phi.lh[s],
           (1 - get(phi0.lh.est[s])$mean*(1 + get(phi0.lh.est[s])$cv^2))/(get(phi0.lh.est[s])$cv^2))
    
    #assigning beta values for beta: beta = (alpha)*(1-mean)/(mean)
    beta.phi.lh[s]<- paste("beta.phi.lh", s, sep = "_")
    
    assign(beta.phi.lh[s],
           get(alpha.phi.lh[s])*(1 - get(phi0.lh.est[s])$mean)/get(phi0.lh.est[s])$mean)
    
    
    phi.lh.a[year,s] <- get(alpha.phi.lh[s]) #alpha shape
    phi.lh.b[year,s] <- get(beta.phi.lh[s]) #beta shape
    
    #effect of removal on transition
    phi.lh1.mean[year,s] <-  get(phi1.lh.est[s])$mean
    phi.lh1.sd[year,] <- get(phi1.lh.est[s])$sd
    
    #phi.hh = transition high to high (beta distribution)
    alpha.phi.hh[s] <- paste("alpha.phi.hh", s, sep = "_")
    #assigning alpha values for beta: alpha = (1-mean)*(1+cv^2)/cv^2
    assign(alpha.phi.hh[s],
           (1 - get(phi0.hh.est[s])$mean*(1 + get(phi0.hh.est[s])$cv^2))/(get(phi0.hh.est[s])$cv^2))
    
    #assigning beta values for beta: beta = (alpha)*(1-mean)/(mean)
    beta.phi.hh[s]<- paste("beta.phi.hh", s, sep = "_")
    
    assign(beta.phi.hh[s],
           get(alpha.phi.hh[s])*(1 - get(phi0.hh.est[s])$mean)/get(phi0.hh.est[s])$mean)
    
    
    phi.hh.a[year,s] <- get(alpha.phi.hh[s]) #alpha shape
    phi.hh.b[year,s] <- get(beta.phi.hh[s]) #beta shape
    
    #phi.hh1
    phi.hh1.mean[year,s] <-  get(phi1.hh.est[s])$mean
    phi.hh1.sd[year,] <- get(phi1.hh.est[s])$sd
    
    # --- p.l ---  detection low state ----------------------------- #
    #p.l.0 = base detection low state (beta distribution)
    alpha.p.l0[s] <- paste("alpha.p.l.0", s, sep = "_")
    #assigning alpha values for beta: alpha = (1-mean)*(1+cv^2)/cv^2
    assign(alpha.p.l0[s],
           (1 - get(p.l0.est[s])$mean*(1 + get(p.l0.est[s])$cv^2))/(get(p.l0.est[s])$cv^2))
    
    #assigning beta values for beta: beta = (alpha)*(1-mean)/(mean)
    beta.p.l0[s]<- paste("beta.p.l0", s, sep = "_")
    
    assign(beta.p.l0[s],
           get(alpha.p.l0[s])*(1 - get(p.l0.est[s])$mean)/get(p.l0.est[s])$mean)
    
    p.l0.a[year,s] <- get(alpha.p.l0[s]) #alpha shape
    p.l0.b[year,s] <- get(beta.p.l0[s]) #beta shape
    
    #p.l.1 = effect of effort (normal distribution)
    p.l1.mean[year,s] <- get(p.l1.est[s])$mean #mean
    p.l1.sd[year,s] <- get(p.l1.est[s])$sd #sd
    
    # --- alpha.l --- difference in baseline detection btwn dat D and M -- #   
    l.mean[year,] <- get(alpha.l.est[s])$mean  #mean
    l.sd[year,] <- get(alpha.l.est[s])$sd
    
    # --- p.h ---  detection high state -------------------------- #
    alpha.p.h0[s] <- paste("alpha.p.h.0", s, sep = "_")
    #assigning alpha values for beta: alpha = (1-mean)*(1+cv^2)/cv^2
    assign(alpha.p.h0[s],
           (1 - get(p.h0.est[s])$mean*(1 + get(p.h0.est[s])$cv^2))/(get(p.h0.est[s])$cv^2))
    
    #assigning beta values for beta: beta = (alpha)*(1-mean)/(mean)
    beta.p.h0[s]<- paste("beta.p.h0", s, sep = "_")
    
    assign(beta.p.h0[s],
           get(alpha.p.h0[s])*(1 - get(p.h0.est[s])$mean)/get(p.h0.est[s])$mean)
    
    p.h0.a[year,s] <- get(alpha.p.h0[s]) #alpha shape
    p.h0.b[year,s] <- get(beta.p.h0[s]) #beta shape
    
    #p.h.1 = effect of effort (normal distribution)
    p.h1.mean[year,s] <- get(p.h1.est[s])$mean #mean
    p.h1.sd[year,s] <- get(p.h1.est[s])$sd #sd
    
    # --- alpha.h --- difference in baseline detection btwn dat D and M -- #   
    h.mean[year,] <- get(alpha.h.est[s])$mean  #mean
    h.sd[year,] <- get(alpha.h.est[s])$sd
    
    # --- S.init and D.init ---  Initial states -------------------- #
    alpha <- rep(1,n.states) #initial state probability vector
    
    } #ends simulation loop 
    
  } #ends year 1+ priors
  
  #--------------------------------------------------------------------------------#
  ###### 2b. JAGS data ######

  #sites where removal occurred
  rem.vec.dat <- rem.vec[,,year,] 
  rem.vec.dat[is.na(rem.vec.dat)] <- 0 #replaces na with 0

  #Parameters monitored
  parameters.to.save <- c("eps.l0", "eps.l1", "eps.h0", "eps.h1", "gamma.0", "gamma.1",
                          "gamma.2", "phi0.lh", "phi1.lh", "phi0.hh", "phi1.hh", "phi.hh", 
                          "p.l0", "p.l1", "p.h0", "p.h1", "State.fin", "alpha.l", "alpha.h", "psi")
  
  #settings
  n.burnin <- 100
  n.iter <- 1000 + n.burnin
  n.chains <- 3
  n.thin <- 1
  
  for(s in 1:n.sims){
    my.data[[s]] <- list( #constants
                         n.sites = n.sites,
                         n.weeks = n.weeks,
                         n.occs = n.occs, 
                         neighbors = neighbors,

                         #data
                         yM= yM[,,,year,s],
                         site.char = site.char,
                         logeffort = logsearch.effort,
                         alpha = alpha,
                         rem.vec = rem.vec.dat[,,s],
                         removal.hours = removal.hours,
                         n.neighbors = n.neighbors,
                         
                         #priors
                         eps.l0.a = eps.l0.a[year,s], 
                         eps.l0.b = eps.l0.a[year,s], 
                         eps.l1.mean = eps.l1.mean[year,s],
                         eps.l1.sd= eps.l1.sd[year,s],
                         eps.h0.a = eps.h0.a[year,s],
                         eps.h0.b= eps.h0.b[year,s],
                         eps.h1.mean= eps.h1.mean[year,s],
                         eps.h1.sd= eps.h1.sd[year,s],
                         gamma.0.mean= gamma.0.mean[year,s],
                         gamma.0.sd= gamma.0.sd[year,s],
                         gamma.1.mean= gamma.1.mean[year,s],
                         gamma.1.sd= gamma.1.sd[year,s],
                         gamma.2.mean = gamma.2.mean[year,s],
                         gamma.2.sd = gamma.2.sd[year,s],
                         phi.lh.a = phi.lh.a[year,s],
                         phi.lh.b = phi.lh.b[year,s],
                         phi.lh1.mean = phi.lh1.mean[year,s],
                         phi.lh1.sd = phi.lh1.sd[year,s],
                         phi.hh.a = phi.hh.a[year,s],
                         phi.hh.b = phi.hh.b[year,s],
                         phi.hh1.mean = phi.hh1.mean[year,s],
                         phi.hh1.sd = phi.hh1.sd[year,s],
                         p.l0.a = p.l0.a[year,s],
                         p.l0.b = p.l0.b[year,s],
                         p.l1.mean = p.l1.mean[year,s],
                         p.l1.sd = p.l1.sd[year,s],
                         l.mean = l.mean[year,s], 
                         l.sd = l.sd[year,s],
                         p.h0.a = p.h0.a[year,s],
                         p.h0.b = p.h0.b[year,s],
                         p.h1.mean = p.h1.mean[year,s],
                         p.h1.sd = p.h1.sd[year,s],
                         h.mean = h.mean[year,s], 
                         h.sd = h.sd[year,s]
    )
  }
  
  ###### 2c. Run JAGS #####
  
  State.start <- array(NA, c(n.sites,n.weeks,n.sims)) #State initial values
  
  for(s in 1:n.sims){
    for(i in 1:n.sites){
      for(week in 1:n.weeks){
        if(rem.vec.dat[i,week,s] == 1){
          State.start[i,week,s] <- max(yM[i,,week,year,s], na.rm = T)
        }else{
        State.start[i,week,s] <- 2
        }
      }
    }
  }

  #Initial values
  for(s in 1:n.sims){
    initial.values[[s]] <- function()list(State = State.start[,,s])
  }
    
  #Running the model
  for(s in 1:n.sims){
    outs[s]<- paste("out", s, sep = "_")
    assign(outs[s],
           jagsUI::jags(data = my.data[[s]],inits = initial.values[[s]],
                        parameters.to.save = parameters.to.save, model.file = "Flower_multistate_datM.txt",
                        n.chains = n.chains, n.thin = n.thin, n.iter = n.iter , n.burnin = n.burnin))
  }
  
  #--------------------------------------------------------------------------------#
  #### 3. Decision for next year ####

  ###### 3a. Save data from MCMC  #####
  for(s in 1:n.sims){ 
    outputsfull[s]<- paste("outputfull", s, sep = "_")
    assign(outputsfull[s], 
           get(outs[s]))
    
    outputs[s]<- paste("output", s, sep = "_")
    assign(outputs[s], 
           as.data.frame((get(outputsfull[s]))$summary))
    
    assign(outputs[s], 
           cbind(get(outputs[s]), param = rownames(get(outputs[s]))))
  }
  
  #Save mcmcs
  for(s in 1:n.sims){
    mcmcs[s]<- paste("mcmc", s, sep = "_")
    assign(mcmcs[s], get(outs[s])$samples)
  }
  
 
  #----- EXTRACTING PARAMETERS FROM THE MODEL ----- #
  #-------State.fin -------# #estimated final state
  for(s in 1:n.sims){
    State.est[s]<- paste("State.est", s, sep = "_")
    assign(State.est[s], filter(get(outputs[s]), grepl("State", param)))
    
    sites[[s]] <- as.numeric(str_nth_number((get(State.est[s]))$param, n = 1))
    

    assign(State.est[s], 
           cbind(get(State.est[s]), site = sites[[s]])) #adding site column  
  
  # --- eps.l ---  eradication low state --- #
  #eps.l0 = base eradication at low state (beta distribution)

    eps.l0.est[s]<- paste("eps.l0", s, sep = "_")
    assign(eps.l0.est[s], filter(get(outputs[s]), grepl("eps.l0", param)))
    
    assign(eps.l0.est[s], 
           cbind(get(eps.l0.est[s]), cv = get(eps.l0.est[s])$sd/get(eps.l0.est[s])$mean
           ))  
  
  #eps.l1 = effect of eradication at low state (normal distribution)
    eps.l1.est[s]<- paste("eps.l1", s, sep = "_")
    assign(eps.l1.est[s], filter(get(outputs[s]), grepl("eps.l1", param)))
  
  
  # --- eps.h ---  eradication high state --- #
  #eps.h0 = base eradication at high state (beta distribution)
    eps.h0.est[s]<- paste("eps.h0", s, sep = "_")
    assign(eps.h0.est[s], filter(get(outputs[s]), grepl("eps.h0", param)))
    
    assign(eps.h0.est[s], 
           cbind(get(eps.h0.est[s]), cv = get(eps.h0.est[s])$sd/get(eps.h0.est[s])$mean
           ))  
  
  #eps.h1 = effect of eradication at high state (normal distribution)
    eps.h1.est[s]<- paste("eps.h1", s, sep = "_")
    assign(eps.h1.est[s], filter(get(outputs[s]), grepl("eps.h1", param)))
 
  # --- gamma ---  invasion --- #  
  #gamma.0 = intrinsic invasion (normal distribution)
    gamma.0.est[s]<- paste("gamma.0", s, sep = "_")
    assign(gamma.0.est[s], filter(get(outputs[s]), grepl("gamma.0", param)))
  
  #gamma.1 = effect of site characteristics (normal distribution)
    gamma.1.est[s]<- paste("gamma.1", s, sep = "_")
    assign(gamma.1.est[s], filter(get(outputs[s]), grepl("gamma.1", param)))

  #gamma.2 = effect of neighboring state (normal distribution)
    gamma.2.est[s]<- paste("gamma.2", s, sep = "_")
    assign(gamma.2.est[s], filter(get(outputs[s]), grepl("gamma.2", param)))
  
  # --- phi ---  transition rates --- #
  #phi.lh = transition low to high (beta distribution)
    phi0.lh.est[s]<- paste("phi0.lh", s, sep = "_")
    assign(phi0.lh.est[s], filter(get(outputs[s]), grepl("phi0.lh", param)))
    
    assign(phi0.lh.est[s], 
           cbind(get(phi0.lh.est[s]), cv = get(phi0.lh.est[s])$sd/get(phi0.lh.est[s])$mean
           ))  
    
    phi1.lh.est[s]<- paste("phi1.lh", s, sep = "_")
    assign(phi1.lh.est[s], filter(get(outputs[s]), grepl("phi1.lh", param)))
  
  #phi.hh = transition high to high (beta distribution)
     phi0.hh.est[s]<- paste("phi0.hh", s, sep = "_")
     assign(phi0.hh.est[s], filter(get(outputs[s]), grepl("phi0.hh", param)))
     
     assign(phi0.hh.est[s], 
            cbind(get(phi0.hh.est[s]), cv = get(phi0.hh.est[s])$sd/get(phi0.hh.est[s])$mean
            ))  
     
     phi1.hh.est[s]<- paste("phi1.hh", s, sep = "_")
     assign(phi1.hh.est[s], filter(get(outputs[s]), grepl("phi1.hh", param)))
  
  # --- p.l ---  detection low state --- #
  #p.l.0 = base detection low state (beta distribution)
    p.l0.est[s]<- paste("p.l0", s, sep = "_")
    assign(p.l0.est[s], filter(get(outputs[s]), grepl("p.l0", param)))
    
    assign(p.l0.est[s], 
           cbind(get(p.l0.est[s]), cv = get(p.l0.est[s])$sd/get(p.l0.est[s])$mean
           ))  

  #p.l.1 = effect of effort (normal distribution)
    p.l1.est[s]<- paste("p.l1", s, sep = "_")
    assign(p.l1.est[s], filter(get(outputs[s]), grepl("p.l1", param)))
  
  # --- alpha.l ---  difference between detection in data D and M --- #
  #alpha.l = effect of neighboring state (normal distribution)
    alpha.l.est[s]<- paste("alpha.l", s, sep = "_")
    assign(alpha.l.est[s], filter(get(outputs[s]), grepl("alpha.l", param)))

  # --- p.h ---  detection high state --- #
  #p.h.0 = base detection high state (beta distribution)
    p.h0.est[s]<- paste("p.h0", s, sep = "_")
    assign(p.h0.est[s], filter(get(outputs[s]), grepl("p.h0", param)))
    
    assign(p.h0.est[s], 
           cbind(get(p.h0.est[s]), cv = get(p.h0.est[s])$sd/get(p.h0.est[s])$mean
           ))  
  
  #p.h.1 = effect of effort (normal distribution)
    p.h1.est[s]<- paste("p.h1", s, sep = "_")
    assign(p.h1.est[s], filter(get(outputs[s]), grepl("p.h1", param)))
  
  # --- alpha.h ---  difference between detection in data D and M --- #
  #alpha.h = (normal distribution)
    alpha.h.est[s]<- paste("alpha.h", s, sep = "_")
    assign(alpha.h.est[s], filter(get(outputs[s]), grepl("alpha.h", param)))
    
  #save annual data
    assign(State.est[s], 
           cbind(get(State.est[s]), year = year))
    
    assign(eps.l0.est[s], 
           cbind(get(eps.l0.est[s]), year = year))
    
    assign(eps.l1.est[s], 
           cbind(get(eps.l1.est[s]), year = year))
    
    assign(eps.h0.est[s], 
           cbind(get(eps.h0.est[s]), year = year))
    
    assign(eps.h1.est[s], 
           cbind(get(eps.h1.est[s]), year = year))
    
    assign(gamma.0.est[s], 
           cbind(get(gamma.0.est[s]), year = year))
    
    assign(gamma.1.est[s], 
           cbind(get(gamma.1.est[s]), year = year))
    
    assign(gamma.2.est[s], 
           cbind(get(gamma.2.est[s]), year = year))
    
    assign(phi0.lh.est[s], 
           cbind(get(phi0.lh.est[s]), year = year))
    
    assign(phi1.lh.est[s], 
           cbind(get(phi1.lh.est[s]), year = year))
    
    assign(phi0.hh.est[s], 
           cbind(get(phi0.hh.est[s]), year = year))
    
    assign(phi1.hh.est[s], 
           cbind(get(phi1.hh.est[s]), year = year))
    
    assign(p.l0.est[s], 
           cbind(get(p.l0.est[s]), year = year))
    
    assign(p.l1.est[s], 
           cbind(get(p.l1.est[s]), year = year))
    
    assign(alpha.l.est[s], 
           cbind(get(alpha.l.est[s]), year = year))
    
    assign(p.h0.est[s], 
           cbind(get(p.h0.est[s]), year = year))
    
    assign(p.h1.est[s], 
           cbind(get(p.h1.est[s]), year = year))
    
    assign(alpha.h.est[s], 
           cbind(get(alpha.h.est[s]), year = year))
    
    
    all.State.est[s]<- paste("States.allsummary", s, sep = "_")
    all.eps.l0.est[s]<- paste("eps.l0.allsummary", s, sep = "_")
    all.eps.l1.est[s]<- paste("eps.l1.allsummary", s, sep = "_")
    all.eps.h0.est[s]<- paste("eps.h0.allsummary", s, sep = "_")
    all.eps.h1.est[s]<- paste("eps.h1.allsummary", s, sep = "_")
    all.gamma.0.est[s]<- paste("gamma.0.allsummary", s, sep = "_")
    all.gamma.1.est[s]<- paste("gamma.1.allsummary", s, sep = "_")
    all.gamma.2.est[s]<- paste("gamma.2.allsummary", s, sep = "_")
    
    all.phi0.lh.est[s]<- paste("phi0.lh.allsummary", s, sep = "_")
    all.phi1.lh.est[s]<- paste("phi1.lh.allsummary", s, sep = "_")
    all.phi0.hh.est[s]<- paste("phi0.hh.allsummary", s, sep = "_")
    all.phi1.hh.est[s]<- paste("phi1.hh.allsummary", s, sep = "_")
    
    all.p.l0.est[s]<- paste("p.l0.allsummary", s, sep = "_")
    all.p.l1.est[s]<- paste("p.l1.allsummary", s, sep = "_")
    all.alpha.l.est[s]<- paste("alpha.l.allsummary", s, sep = "_")
    all.p.h0.est[s]<- paste("p.h0.allsummary", s, sep = "_")
    all.p.h1.est[s]<- paste("p.h1.allsummary", s, sep = "_")
    all.alpha.h.est[s]<- paste("alpha.h.allsummary", s, sep = "_")
    
    
    #If year 1 we set summary data frame to itself
    if(year == 1){
      assign(all.State.est[s], 
             get(State.est[s]))
      
      assign(all.eps.l0.est[s], 
             get(eps.l0.est[s]))
      
      assign(all.eps.l1.est[s], 
             get(eps.l1.est[s]))
      
      assign(all.eps.h0.est[s], 
             get(eps.h0.est[s]))
      
      assign(all.eps.h1.est[s], 
             get(eps.h1.est[s]))
      
      assign(all.gamma.0.est[s], 
             get(gamma.0.est[s]))
      
      assign(all.gamma.1.est[s], 
             get(gamma.1.est[s]))
      
      assign(all.gamma.2.est[s], 
             get(gamma.2.est[s]))
      
      assign(all.phi0.lh.est[s], 
             get(phi0.lh.est[s]))
      
      assign(all.phi1.lh.est[s], 
             get(phi1.lh.est[s]))
      
      assign(all.phi0.hh.est[s], 
             get(phi0.hh.est[s]))
      
      assign(all.phi1.hh.est[s], 
             get(phi1.hh.est[s]))
      
      assign(all.p.l0.est[s], 
             get(p.l0.est[s]))
      
      assign(all.p.l1.est[s], 
             get(p.l1.est[s]))
      
      assign(all.alpha.l.est[s], 
             get(alpha.l.est[s]))
      
      assign(all.p.h0.est[s], 
             get(p.h0.est[s]))
      
      assign(all.p.h1.est[s], 
             get(p.h1.est[s]))
      
      assign(all.alpha.h.est[s], 
             get(alpha.h.est[s]))
      
      
      
    }else{ #if beyond first year, we append previous summary to new summary
      assign(all.State.est[s], 
             rbind(get(all.State.est[s]), get(State.est[s])))
      
      assign(all.eps.l0.est[s], 
             rbind(get(all.eps.l0.est[s]), get(eps.l0.est[s])))
      
      assign(all.eps.l1.est[s], 
             rbind(get(all.eps.l1.est[s]), get(eps.l1.est[s])))
      
      assign(all.eps.h0.est[s], 
             rbind(get(all.eps.h0.est[s]), get(eps.h0.est[s])))
      
      assign(all.eps.h1.est[s], 
             rbind(get(all.eps.h1.est[s]), get(eps.h1.est[s])))
      
      assign(all.gamma.0.est[s], 
             rbind(get(all.gamma.0.est[s]), get(gamma.0.est[s])))
      
      assign(all.gamma.1.est[s], 
             rbind(get(all.gamma.1.est[s]), get(gamma.1.est[s])))
      
      assign(all.gamma.2.est[s], 
             rbind(get(all.gamma.2.est[s]), get(gamma.2.est[s])))
    
      assign(all.phi0.lh.est[s], 
             rbind(get(all.phi0.lh.est[s]), get(phi0.lh.est[s])))
      
      assign(all.phi1.lh.est[s], 
             rbind(get(all.phi1.lh.est[s]), get(phi1.lh.est[s])))
      
      assign(all.phi0.hh.est[s], 
             rbind(get(all.phi0.hh.est[s]), get(phi0.hh.est[s])))
      
      assign(all.phi1.hh.est[s], 
             rbind(get(all.phi1.hh.est[s]), get(phi1.hh.est[s])))
      
      assign(all.p.l0.est[s], 
             rbind(get(all.p.l0.est[s]), get(p.l0.est[s])))
      
      assign(all.p.l1.est[s], 
             rbind(get(all.p.l1.est[s]), get(p.l1.est[s])))
      
      assign(all.alpha.l.est[s], 
             rbind(get(all.alpha.l.est[s]), get(alpha.l.est[s])))
      
      assign(all.p.h0.est[s], 
             rbind(get(all.p.h0.est[s]), get(p.h0.est[s])))      
      
      assign(all.p.h1.est[s], 
             rbind(get(all.p.h1.est[s]), get(p.h1.est[s])))
      
      assign(all.alpha.h.est[s], 
             rbind(get(all.alpha.h.est[s]), get(alpha.h.est[s])))
      
      
    }
  
  }
  
  
  #--------------------------------------------------------------------------------#
  ###### 3b. Make decision  #####
  # I assume that the state of the system stays the same after 48 weeks 
  # Thus, removal locations are selected by the estimated sites from the model
  
  for(s in 1:n.sims){
    #extracting and saving states data from the model
    States.mean.round[,year,s] <- round((get(State.est[s]))$mean)
    States.mean[,year,s] <- get(State.est[s])$mean #not rounded
  }
  
  States.mean.wide <- data.frame(t(States.mean[,year,]))
  colnames(States.mean.wide) <- seq(1:n.sites)
  States.mean.wide$sim <- seq(1:n.sims)  
  
  States.mean.long <- gather(States.mean.wide, site, state, all_of(n.sites), factor_key=TRUE)
  States.mean.long$site <- as.numeric(States.mean.long$site)
  
  States.mean.long$year <- year
  
  if(year == 1){
    States.mean.years <- States.mean.long
  }else{
    States.mean.years <- rbind(States.mean.years, States.mean.long)
  }
  
  ##---- making the decision for next year: ----##
  if(year < n.years){

    #Removal locations: rank sites by state
     for(s in 1:n.sims){
      sites.rem.M[,1,year+1,s] <- order(States.mean[,year,s], decreasing = T)
    }
    
  }else{
    #during the final year, we project the final state
    for(s in 1:n.sims){
      
      #simulating truth at the end
      for(i in 1:n.sites){
        State[i,1,year+1,s] <- rcat(1,TPM[State[i,4,year,s], i, 4, year, s,]) 
      }
      
      #estimating final state via results from estimation model
      for(i in 1:n.sites){
        D.est[i,s] <- sum(States.mean.round[neighbors[i,],year,s])/2 #state of neighbors
        
        gamma.est[i,s] <-invlogit(get(gamma.0.est[s])$mean + get(gamma.1.est[s])$mean*site.char[i] + get(gamma.2.est[s])$mean*D.est[i,s]) 
        eps.l.est[i,s] <- invlogit(get(eps.l0.est[s])$mean + get(eps.l1.est[s])$mean*rem.vec.dat[i,4,s]*removal.hours[2]) 
        eps.h.est[i,s] <- invlogit(get(eps.h0.est[s])$mean + get(eps.h1.est[s])$mean*rem.vec.dat[i,4,s]*removal.hours[3]) 
        phi.lh.est[i,s] <- invlogit(get(phi0.lh.est[s])$mean + get(phi1.lh.est[s])$mean*rem.vec.dat[i,4,s]*removal.hours[2])
        phi.hh.est[i,s] <- invlogit(get(phi0.hh.est[s])$mean + get(phi1.hh.est[s])$mean*rem.vec.dat[i,4,s]*removal.hours[3])
        
        
        TPM.est[1,i,s,1] <- 1-gamma.est[i,s]
        TPM.est[1,i,s,2] <- gamma.est[i,s]
        TPM.est[1,i,s,3] <- 0
        TPM.est[2,i,s,1] <- eps.l.est[i,s]
        TPM.est[2,i,s,2] <- (1- eps.l.est[i,s])*(1-phi.lh.est[i,s])
        TPM.est[2,i,s,3] <- (1- eps.l.est[i,s])*phi.lh.est[i,s] 
        TPM.est[3,i,s,1] <- eps.h.est[i,s]
        TPM.est[3,i,s,2] <- (1- eps.h.est[i,s])*(1-phi.hh.est[i,s])
        TPM.est[3,i,s,3] <- (1- eps.h.est[i,s])*phi.hh.est[i,s]
        
      }
      
      for(i in 1:n.sites){
      
        #if we visited the site for observation data
        if(!is.na(rem.vec[i,4,year,s])){ 
          S.end[i,s] <- max(yM[i,,4,year,s], na.rm = T) 
          
        }else{
          S.end[i,s] <- rcat(1,TPM.est[States.mean.round[i,year,s], i, s,]) 
        }
        
      } #sites loop
      
    } #sims loop
  } #final year loop

} #end adaptive management 

#################################################################################################
#### TIMING ####
end.time <- Sys.time()
time.taken <- end.time - start.time

#### Save True Data ####
#results for each sim
States.df <- adply(State, c(1,2,3,4))
colnames(States.df) <- c("site", "week", "year", "sim", "state")              

#mean across simulations
Mean.States.df <- aggregate(state ~ site+week+year,
                            data = as.data.frame(States.df), FUN = mean)

#observation data -multi
yM.df <- adply(yM, c(1,2,3,4,5))
colnames(yM.df) <- c("site", "occasion", "week", "year", "sim", "observed.state")              
rem.site.M.df <- yM.df %>% filter(observed.state > 1)

#### sites visited ####
sites.visit <- adply(rem.vec, c(1,2,4,3))
colnames(sites.visit) <- c("site", "week", "year", "sim", "rem.val")   
sites.visit <- sites.visit %>% filter(!is.na(rem.val))

#visit no remove
sites.visit.norem <- sites.visit %>% filter(rem.val == 0)
sites.visit.norem$rem.val <- 1
sites.visit.norem <- aggregate(rem.val ~ week + year + sim,
                               data = as.data.frame(sites.visit.norem), FUN = sum)


sites.visit.norem.avg <- aggregate(rem.val ~ week+ year,
                                   data = as.data.frame(sites.visit.norem), FUN = mean)

colnames(sites.visit.norem.avg)[3] <- "num.visit.norem"

#visit remove
sites.visit.rem <- sites.visit %>% filter(rem.val == 1)

sites.visit.rem<- aggregate(rem.val ~ week+ year + sim,
                            data = as.data.frame(sites.visit.rem), FUN = sum)

sites.visit.rem.avg <- aggregate(rem.val ~ week + year,
                                 data = as.data.frame(sites.visit.rem), FUN = mean)


colnames(sites.visit.rem.avg)[3] <- "num.visit.rem"

sites.df <- cbind(sites.visit.norem.avg, num.visit.rem = sites.visit.rem.avg$num.visit.rem)

#### Estimated Data ####
##### Estimated States ####
States.est.df <- States.mean.years %>% select(site,year,sim,state)

#mean across simulations
Mean.States.est.df <- aggregate(state ~ site+year,
                                data = as.data.frame(States.est.df), FUN = mean)

##### Estimated parameters ####
## --- eps.l0 -----------------------------------------------#
eps.l0s <- list()

for(s in 1:n.sims){
  assign(all.eps.l0.est[s], 
         cbind(get(all.eps.l0.est[s]), sim = s))
  
  eps.l0s[[s]] <- get(all.eps.l0.est[s])
}

eps.l0s.df <- do.call("rbind", eps.l0s)

## --- eps.l1 -----------------------------------------------#
eps.l1s <- list()

for(s in 1:n.sims){
  assign(all.eps.l1.est[s], 
         cbind(get(all.eps.l1.est[s]), sim = s))
  
  eps.l1s[[s]] <- get(all.eps.l1.est[s])
}

eps.l1s.df <- do.call("rbind", eps.l1s)

## --- eps.h0 -----------------------------------------------#
eps.h0s <- list()

for(s in 1:n.sims){
  assign(all.eps.h0.est[s], 
         cbind(get(all.eps.h0.est[s]), sim = s))
  
  eps.h0s[[s]] <- get(all.eps.h0.est[s])
}

eps.h0s.df <- do.call("rbind", eps.h0s)

## --- eps.h1 -----------------------------------------------#
eps.h1s <- list()

for(s in 1:n.sims){
  assign(all.eps.h1.est[s], 
         cbind(get(all.eps.h1.est[s]), sim = s))
  
  eps.h1s[[s]] <- get(all.eps.h1.est[s])
}

eps.h1s.df <- do.call("rbind", eps.h1s)

## --- gamma.0 -----------------------------------------------#
gamma.0s <- list()

for(s in 1:n.sims){
  assign(all.gamma.0.est[s], 
         cbind(get(all.gamma.0.est[s]), sim = s))
  
  gamma.0s[[s]] <- get(all.gamma.0.est[s])
}

gamma.0s.df <- do.call("rbind", gamma.0s)

## --- gamma.1 -----------------------------------------------#
gamma.1s <- list()

for(s in 1:n.sims){
  assign(all.gamma.1.est[s], 
         cbind(get(all.gamma.1.est[s]), sim = s))
  
  gamma.1s[[s]] <- get(all.gamma.1.est[s])
}

gamma.1s.df <- do.call("rbind", gamma.1s)

## --- gamma.2 -----------------------------------------------#
gamma.2s <- list()

for(s in 1:n.sims){
  assign(all.gamma.2.est[s], 
         cbind(get(all.gamma.2.est[s]), sim = s))
  
  gamma.2s[[s]] <- get(all.gamma.2.est[s])
}

gamma.2s.df <- do.call("rbind", gamma.2s)

## --- phi.lh -----------------------------------------------#
phi0.lhs <- list()
phi1.lhs <- list()

for(s in 1:n.sims){
  assign(all.phi0.lh.est[s], 
         cbind(get(all.phi0.lh.est[s]), sim = s))
  
  phi0.lhs[[s]] <- get(all.phi0.lh.est[s])
  
  assign(all.phi1.lh.est[s], 
         cbind(get(all.phi1.lh.est[s]), sim = s))
  
  phi1.lhs[[s]] <- get(all.phi1.lh.est[s])
}

phi0.lhs.df <- do.call("rbind", phi0.lhs)
phi1.lhs.df <- do.call("rbind", phi1.lhs)

## --- phi.hh -----------------------------------------------#
phi0.hhs <- list()
phi1.hhs <- list()

for(s in 1:n.sims){
  assign(all.phi0.hh.est[s], 
         cbind(get(all.phi0.hh.est[s]), sim = s))
  
  phi0.hhs[[s]] <- get(all.phi0.hh.est[s])
  
  assign(all.phi1.hh.est[s], 
         cbind(get(all.phi1.hh.est[s]), sim = s))
  
  phi1.hhs[[s]] <- get(all.phi1.hh.est[s])
}

phi0.hhs.df <- do.call("rbind", phi0.hhs)
phi1.hhs.df <- do.call("rbind", phi1.hhs)

## --- p.l0 -----------------------------------------------#
p.l0.s <- list()

for(s in 1:n.sims){
  assign(all.p.l0.est[s], 
         cbind(get(all.p.l0.est[s]), sim = s))
  
  p.l0.s[[s]] <- get(all.p.l0.est[s])
}

p.l0.s.df <- do.call("rbind", p.l0.s)

## --- p.l1 -----------------------------------------------#
p.l1.s <- list()

for(s in 1:n.sims){
  assign(all.p.l1.est[s], 
         cbind(get(all.p.l1.est[s]), sim = s))
  
  p.l1.s[[s]] <- get(all.p.l1.est[s])
}

p.l1.s.df <- do.call("rbind", p.l1.s)

## --- p.h0 -----------------------------------------------#
p.h0.s <- list()

for(s in 1:n.sims){
  assign(all.p.h0.est[s], 
         cbind(get(all.p.h0.est[s]), sim = s))
  
  p.h0.s[[s]] <- get(all.p.h0.est[s])
}

p.h0.s.df <- do.call("rbind", p.h0.s)

## --- p.h1 -----------------------------------------------#
p.h1.s <- list()

for(s in 1:n.sims){
  assign(all.p.h1.est[s], 
         cbind(get(all.p.h1.est[s]), sim = s))
  
  p.h1.s[[s]] <- get(all.p.h1.est[s])
}

p.h1.s.df <- do.call("rbind", p.h1.s)

#### QUICK RESULTS ####
Mean.States.df.fin <- Mean.States.df %>% filter(year == 11)
Mean.States.df.fin$state #distribution

mean(Mean.States.df.fin$state) #average final state

#check number of sites visited for removal on average each week
mean(sites.df$num.visit.norem)
mean(sites.df$num.visit.rem)

#correct results: true state
match <- array(NA, c(n.sites, n.weeks, n.years, n.sims))
match.dat <- array(NA, c(n.weeks, n.years, n.sims))

for(s in 1:n.sims){
  for(year in 1:n.years){
    for(week in 1:n.weeks){
      State.M <- State[,week,year,s]
      full.match <- (yM[,,week,year,s] == State.M)
      full.match [,1] <- as.numeric(full.match [,1])
      full.match [,2] <- as.numeric(full.match [,2])
      full.match [is.na(full.match )] <- 3 #replace NA with 3
      
      
      for(i in 1:n.sites){
        
        
        if(full.match[i,1] == 1 & full.match[i,1] == 3){ #true match first try
          match[i,week,year,s] <- 1
        }
        
        if(full.match[i,1] == 1 & full.match[i,2] == 1){ #true match
          match[i,week,year,s] <- 1
        }
        
        if(full.match[i,1] == 0 & full.match[i,2] == 1){ #true match on the second try
          match[i,week,year,s] <- 1
        }
        
        if(full.match[i,1] == 0 & full.match[i,2] == 0){ #not correct
          match[i,week,year,s] <- 0
        }
        
        if(full.match[i,1] == 3 & full.match[i,2] == 3){ #true match on the second try
          match[i,week,year,s] <- NA #not visited
        }
        
        
      } #sites
      
      match2 <- discard(match[,week,year,s], is.na)
      match.dat[week,year,s] <- sum(match2 == 1)/ length(match2) 
      
    } #weeks
  } #year
} #sims


mean(match.dat)

#correct results: detection/non-detection
match <- array(NA, c(n.sites, n.weeks, n.years, n.sims))
match.dat <- array(NA, c(n.weeks, n.years, n.sims))

for(s in 1:n.sims){
  for(year in 1:n.years){
    for(week in 1:n.weeks){
      State.D <- State[,week,year,s]
      State.D[State.D == 3] <- 2
      yM.2 <- yM[,,week,year,s]
      yM.2[yM.2 == 3] <- 2
      full.match <- (yM.2 == State.D)
      full.match [,1] <- as.numeric(full.match [,1])
      full.match [,2] <- as.numeric(full.match [,2])
      full.match [is.na(full.match )] <- 3 #replace NA with 3
      
      
      for(i in 1:n.sites){
        
        
        if(full.match[i,1] == 1 & full.match[i,1] == 3){ #true match first try
          match[i,week,year,s] <- 1
        }
        
        if(full.match[i,1] == 1 & full.match[i,2] == 1){ #true match
          match[i,week,year,s] <- 1
        }
        
        if(full.match[i,1] == 0 & full.match[i,2] == 1){ #true match on the second try
          match[i,week,year,s] <- 1
        }
        
        if(full.match[i,1] == 0 & full.match[i,2] == 0){ #not correct
          match[i,week,year,s] <- 0
        }
        
        if(full.match[i,1] == 3 & full.match[i,2] == 3){ #true match on the second try
          match[i,week,year,s] <- NA #not visited
        }
        
        
      } #sites
      
      match2 <- discard(match[,week,year,s], is.na)
      match.dat[week,year,s] <- sum(match2 == 1)/ length(match2) 
      
    } #weeks
  } #year
} #sims

mean(match.dat)

#### Final States average state ####
State.fins <- State[,4,n.years,]
State.fins.df <- adply(State.fins, c(1,2))
colnames(State.fins.df) <- c("site","sim", "state")

ggplot(State.fins.df)+
  geom_boxplot(mapping = aes(y = state, middle = mean(state)))

summary(State.fins.df$state)

#### site invasion ####
State.fins.avg <- aggregate(state ~ site, State.fins.df, mean)

head(State.fins.avg)

ggplot(State.fins.avg, aes(x = site, y = 1, fill = state)) +
  geom_tile()+
  theme_classic()

#### number of invaded sites ####
invasion <- rep(NA, n.sims)

for(s in 1:n.sims){
  df <- filter(State.fins.df, sim == s)
  invasion[s] <- sum(df$state == 1)
}

invasion.mean <- mean(invasion)


#percent of river uninvaded after 10 years
1- invasion.mean/n.sites

