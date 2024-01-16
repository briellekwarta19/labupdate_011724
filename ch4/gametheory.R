library(tidyverse)
library(plyr)
library(RColorBrewer) 
library(Recon)
library(rPref)

players <- 2 #Two players are the manager and neighbor
N.years <- 20

#### Biological Model ####
#harvest rate
hM <- seq(0,1,by = 0.5) 
hN <- seq(0,1,by = 0.5)

h <- expand.grid(hM,hN) #creating combinations of harvest rates (alternatives)
colnames(h) <- c("M", "N")
N.h <- length(h$M) #total number of management alternatives

X0 <- c(100,100) #initial population in manager, neighbor area respectively
X <- X.rem <- C<- G <- R <- array(0, dim = c(players, N.years, N.h)) #arrays
X[,1,] <- X0 #population at the first time step 

#Rates for biological model
r_M <- 1 #growth rate manager
r_N <- 1 #growth rate neighbor
K_M <- 300 #carrying capacity manager
K_N <- 300 #carrying capacity neighbor
b <- 0.5 #dispersal rate into other location

for(m in 1:N.h){
  for(t in 1:(N.years-1)){

    #First we calculate the population remaining after removals
    X.rem[1,t,m] <- X[1,t,m]- X[1,t,m]*h$M[m]
    X.rem[2,t,m] <- X[2,t,m] - X[2,t,m]*h$N[m]

    #Manager population:
    X[1,t+1,m] <- r_M*X.rem[1,t,m]*(1-(X.rem[1,t,m]/K_M)) + b*(X.rem[2,t,m]- X.rem[1,t,m]) + X.rem[1,t,m] 
    R[1,t+1,m] <- X[1,t,m]*h$M[m]  #manager removed

    #Neighbor population
    X[2,t+1,m] <- r_N*X.rem[2,t,m]*(1-(X.rem[2,t,m]/K_N)) - b*(X.rem[2,t,m]- X.rem[1,t,m]) + X.rem[2,t,m] 
    R[2,t+1,m] <- X[2,t,m]*h$N[m] #neighbor removed

  }
}


#Abundance plots
X.dat <- adply(X, c(1,2,3))
X.dat <- data.frame(X.dat)
colnames(X.dat) <- c("player", "year", "Management", "count")
X.dat$Management <- as.factor(X.dat$Management)

types <- c("1: (0,0)","2: (0.5,0)", "3: (1,0)", 
                      "4: (0,0.5)", "5: (0.5, 0.5)", "6: (1, 0.5)",
                      "7: (0,1)", "8: (0.5,1)", "9: (1,1)")

col <- brewer.pal(8, "Set1") 

colors <- colorRampPalette(col)(9)

ggplot(X.dat)+
  geom_point(aes(x = year, y = count, shape = player, color = Management), size = 2)+
  scale_shape_manual(values = c(16 ,1), labels = c("M", "N"), name = "Player") +
  geom_line(aes(x = year, y = count, color = Management, group = interaction(player)))+
  scale_color_manual(name = "Management rate (M,N)", labels = types, values = colors) + 
  facet_wrap(~Management, labeller = "label_both") +
  labs(title = "Abundance", x = "Year", y = "Abundance")

#### Cost Model ####
c_b <- 100 #base cost of removals
c_g <- 10 #income per removal 

C.total <- array(NA, c(players,N.h)) #total cost
G.total <- array(NA, c(players,N.h)) #total income
R.total <- array(NA, c(players,N.h)) #total removals

for(m in 1:N.h){
  for(t in 1:N.years){
    
    #Cost of removal
    C[1,t,m] <- log((X[1,t,m]*h$M[m])^2+1) + c_b
    C[2,t,m] <- log((X[2,t,m]*h$N[m])^2+1) + c_b
    
    #Income from removal (only neighbor has this)
    G[2,t,m] <- c_g*X[2,t,m]*h$N[m]
    
  }
  #Sum cost through time:
  C.total[1:2,m] <- c(sum(C[1,,m]), sum(C[2,,m]))
  
  #Sum income through time:
  G.total[1:2,m] <- c(0, sum(G[2,,m]))
  
  #Sum of removals through time
  R.total[1:2,m] <- c(sum(R[1,,m]), sum(R[2,,m]))
}

C.dat <- adply(C, c(1,2,3))
colnames(C.dat) <- c("player", "year", "Management", "cost")
C.dat$Management <- as.factor(C.dat$Management)

ggplot(C.dat)+
  geom_point(aes(x = year, y = cost, shape = player, color = Management), size = 2)+
  scale_shape_manual(values = c(16 ,1), labels = c("M", "N"), name = "Player") +
  geom_line(aes(x = year, y = cost, color = Management, group = interaction(player)))+
  scale_color_manual(name = "Management (M,N)", labels = types, values = colors) +
  facet_wrap(~Management, labeller = "label_both") +labs(title = "Cost", x = "Year", y = "Cost")

G.dat <- adply(G, c(1,2,3))
colnames(G.dat) <- c("player", "year", "Management", "income")
G.dat$Management <- as.factor(G.dat$Management)

ggplot(G.dat)+
  geom_point(aes(x = year, y = income, shape = player, color = Management), size = 2)+
  scale_shape_manual(values = c(16 ,1), labels = c("M", "N"), name = "Player") +
  geom_line(aes(x = year, y = income, color = Management, group = interaction(player)))+
  scale_color_manual(name = "Management (M,N)", labels = types, values = colors) +
  facet_wrap(~Management, labeller = "label_both") +
  labs(title = "Income", x = "Year", y = "Cost")

#### Payoffs Model ####
#The managers payoff model is always the same
#But we will define 4 different models for the neighbor

#-------- 1. Re-scale objective results --------#

#first I need to scale final population, accumulative removals, cost, and income objectives 
#to a 0 to 1 scale with a mean of 0 unit variance

#1.a: Final population:
X.fin <- X[1:2,N.years,1:N.h]

X.scaled <- rbind(c(scale(X.fin[1,])), #scaling X.M final population results mean 0, sd 1
                      c(scale(X.fin[2,]))) #scaling X.N final population results mean 0, sd 1

#1.b: Accumulative removals:
R.total[1:2,] <- R.total[1:2,]

R.scaled <- rbind(c(scale(R.total[1,])), 
                      c(scale(R.total[2,]))) 


#1.c: Accumulative Cost:
C.total[1:2,] <- C.total[1:2,]

C.scaled <- rbind(c(scale(C.total[1,])), 
                        c(scale(C.total[2,]))) 

#1.d: Accumulative Income:
G.total[1:2,] <- G.total[1:2,]

G.scaled <- rbind(rep(0,N.h), #Managers have no income
                        c(scale(G.total[2,]))) 

#-------- 2. Establish weights --------#
#two objectives are:
# 1. An invasive species objective (ex: Final population or accumulative removals)
# 2. An economic objective (ex: Accumulative cost or accumulative income)

#We need to establish weights for objectives they must sum to 1 but cant be 0 or 1
weight_grid <- data.frame(expand.grid(seq(0.1,0.9,0.1),seq(0.1,0.9,0.1)))
weight <- weight_grid[rowSums(weight_grid) == 1, ]
colnames(weight ) <- c("a", "b")
N.weights <- length(weight$a)

#-------- 3. Calculate Payoffs --------#

tot.players <- 5 #Manager, Environmentalist, Hobbyist, Profiteer, Bounty hunter

f <- array(NA, c(tot.players, N.h, N.weights))

for(m in 1:N.h){
  for(w in 1:N.weights){
    #Manager payoff
    f[1,m,w] <- -(weight$a[w]*(X.scaled[1,m]+X.scaled[2,m])) - (weight$b[w]*C.scaled[1,m])
    
    #Neighbor payoffs:
    # Environmentalist payoff
    f[2,m,w] <- -(weight$a[w]*(X.scaled[2,m])) - (weight$b[w]*C.scaled[2,m])
    
    # Hobbyist payoff
    f[3,m,w] <- (weight$a[w]*R.scaled[2,m]) - (weight$b[w]*C.scaled[2,m])
    
    # Hobbyist payoff
    f[4,m,w] <- (weight$a[w]*R.scaled[2,m]) + (weight$b[w]*G.scaled[2,m])
    
    # Bounty hunter payoff
    f[5,m,w] <- -(weight$a[w]*X.scaled[2,m]) + (weight$b[w]*G.scaled[2,m])

  }
}

#### MCDA result ####
#sum across weights
f.sum <- array(NA, dim = c(tot.players,N.h))

for(p in 1:tot.players){
  for(m in 1:N.h){
    f.sum[p,m] <- sum(f[p,m,]) #summing payoff across weights
  }
}



#Best alternative for Manager
h[which.max(f.sum[1,]),] 

#Best alternative for Environmentalist
h[which.max(f.sum[2,]),] 

#Best alternative for Hobbyist
h[which.max(f.sum[3,]),] 

#Best alternative for Profiteer
h[which.max(f.sum[4,]),] 

#Best alternative for Bounty Hunter
h[which.max(f.sum[5,]),] 

payoff.mcda <- adply(f.sum,c(1,2))
colnames(payoff.mcda) <- c("Player", "Management", "Payoff")

payoff.mcda$Player2 <- '1'
for(i in 1:length(payoff.mcda$Player)){
  if(as.numeric(payoff.mcda$Player[i]) > 1){
    payoff.mcda$Player2[i] <- '2'
  } 
}

#Plot of payoffs
ggplot(payoff.mcda)+
  geom_point(aes(x = Player, y = Payoff, color = Management, shape = Player2), size = 2)+
  scale_shape_manual(values = c(16 ,1), labels = c("M", "N"), name = "Player") +
  scale_color_manual(name = "Management (M,N)", labels = types, values = colors) + 
  facet_wrap(~Management, labeller = "label_both") +
  scale_x_discrete("Player", 
                   labels = c(expression(M), expression(N[E]),expression(N[H]),
                              expression(N[P]), expression(N[B])))+
  labs(title = "Payoff", x = "Player", y = "Payoff")


#### Game theory result ####
#expected value
f.exp<- array(NA, dim = c(tot.players,N.h))

for(p in 1:tot.players){
  for(m in 1:N.h){
    f.exp[p,m] <- mean(f[p,m,])
  }
}

#To make the comparison between mcda and game theory use the summed payoffs for game theory
f.game <- f.sum

M.matrix <- matrix(NA, nrow = 3, ncol = 3)
M.matrix[1,1:3] <- c(f.game[1,9],  f.game[1,6], f.game[1,3])
M.matrix[2,1:3] <- c(f.game[1,8], f.game[1,5], f.game[1,2])
M.matrix[3,1:3] <- c(f.game[1,7], f.game[1,4] , f.game[1,1])

N.matrix <- array(NA, dim = c(3,3,4))
N.matrix[1,1,] <- f.game[2:5,9] 
N.matrix[1,2,] <- f.game[2:5,6] 
N.matrix[1,3,] <- f.game[2:5,3] 
N.matrix[2,1,] <- f.game[2:5,8] 
N.matrix[2,2,] <- f.game[2:5,5] 
N.matrix[2,3,] <- f.game[2:5,2] 
N.matrix[3,1,] <- f.game[2:5,7] 
N.matrix[3,2,] <- f.game[2:5,4] 
N.matrix[3,3,] <- f.game[2:5,1] 

##### Nash Equilibrium #####
game <- matrix(NA, nrow = 3, ncol = 3)
game[1,] <- c("1, 1","1, 0.5", "1,0")
game[2,] <- c("0.5, 1","0.5, 0.5", "0.5,0")
game[3,] <- c("0, 1","0, 0.5", "0,0")

#Nash equilibrium manager vs environmentalist
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,1]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,1]))[2])]

#Nash equilibrium manager vs hobbyist
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,2]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,2]))[2])]

#Nash equilibrium manager vs profiteer
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,3]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,3]))[2])]

#Nash equilibrium manager vs bounty hunter
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,4]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,4]))[2])]

#above we set up matrix like:
#                    N
#    | 1.  (1,1)    4. (1,0.5)     7. (1,0)   |
# M  | 2. (0.5,1)   5. (0.5, 0.5)  8. (0.5,0) |
#    | 3. (0,0)     6. (0,0.5)     9. (0,0)   |

##### Pareto optimal #####
#Finding the pareto optimal solution
N <- 4 #Toggle this between 1,2,3,4 to find pareto optimal solutions against all neighbors

x <- as.vector(M.matrix)
y <- as.vector(N.matrix[,,N])
d = data.frame(x,y)
D = d[order(d$x,d$y,decreasing=TRUE),]
front = D[which(!duplicated(cummax(D$y))),]
front

p1 <- psel(d, high(x) * high(y))

# Plotting Pareto optimal solutions (bolded are pareto optimal)
ggplot(d, aes(x = x, y = y)) + 
  geom_point(shape = 21) +   
  geom_point(data = p1, size = 3)

po.title <- c(game[1,1], game[2,1], game[3,1], game[1,2], game[2,2], game[3,2], game[1,3], game[2,3], game[3,3])

#### Full result ####
##### Environmentalist #####
#Pareto optimal solutions
N <- 1 #Toggle this between 1,2,3,4 to find pareto optimal solutions against all neighbors
x <- as.vector(M.matrix)
y <- as.vector(N.matrix[,,N])
d = data.frame(x,y)
D = d[order(d$x,d$y,decreasing=TRUE),]
front = D[which(!duplicated(cummax(D$y))),]
p1 <- psel(d, high(x) * high(y))

#Best alternative for Manager
h[which.max(f.sum[1,]),] 

#Best alternative for Environmentalist
h[which.max(f.sum[2,]),] 

#nash equilibrium
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,1]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,1]))[2])]

#PO solution
po.title[as.numeric(rownames(p1))]

##### Hobbyist #####
#Pareto optimal solutions
N <- 2 #Toggle this between 1,2,3,4 to find pareto optimal solutions against all neighbors
x <- as.vector(M.matrix)
y <- as.vector(N.matrix[,,N])
d = data.frame(x,y)
D = d[order(d$x,d$y,decreasing=TRUE),]
front = D[which(!duplicated(cummax(D$y))),]
p1 <- psel(d, high(x) * high(y))

#Best alternative for Manager
h[which.max(f.sum[1,]),] 

#Best alternative for Environmentalist
h[which.max(f.sum[3,]),] 

#nash equilibrium
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,2]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,2]))[2])]


po.title[as.numeric(rownames(p1))]

##### Profiteer #####
#Pareto optimal solutions
N <- 3 #Toggle this between 1,2,3,4 to find pareto optimal solutions against all neighbors
x <- as.vector(M.matrix)
y <- as.vector(N.matrix[,,N])
d = data.frame(x,y)
D = d[order(d$x,d$y,decreasing=TRUE),]
front = D[which(!duplicated(cummax(D$y))),]
p1 <- psel(d, high(x) * high(y))

#Best alternative for Manager
h[which.max(f.sum[1,]),] 

#Best alternative for Profiteer
h[which.max(f.sum[4,]),] 

#nash equilibrium
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,3]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,3]))[2])]

po.title[as.numeric(rownames(p1))]

##### Bounty Hunter #####
#Pareto optimal solutions
N <- 4 #Toggle this between 1,2,3,4 to find pareto optimal solutions against all neighbors
x <- as.vector(M.matrix)
y <- as.vector(N.matrix[,,N])
d = data.frame(x,y)
D = d[order(d$x,d$y,decreasing=TRUE),]
front = D[which(!duplicated(cummax(D$y))),]
p1 <- psel(d, high(x) * high(y))

#Best alternative for Manager
h[which.max(f.sum[1,]),] 

#Best alternative for bounty hunter
h[which.max(f.sum[5,]),] 

#nash equilibrium
game[as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,4]))[1]),
     as.numeric(unlist(sim_nasheq(M.matrix, N.matrix[,,4]))[2])]

po.title[as.numeric(rownames(p1))]

