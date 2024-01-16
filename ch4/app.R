library(tidyverse)
library(plyr)
library(RColorBrewer) 
library(Recon)
library(rPref)


ui <- fluidPage(

    # Application title
    titlePanel("Game Theory App"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("r",
                        "Growth rate:",
                        min = -1,
                        max = 1,
                        value = 0.5,
                        step = 0.5),
            
            sliderInput("b",
                        "Dispersal rate:",
                        min = 0,
                        max = 1,
                        value = 0.2,
                        step = 0.2),
            
            sliderInput("K",
                        "Carrying Capacity:",
                        min = 100,
                        max = 500,
                        value = 300,
                        step = 100),
            
            sliderInput("M1",
                        "Initial Manager's population:",
                        min = 100,
                        max = 500,
                        value = 100,
                        step = 100),
            
            sliderInput("N1",
                        "Initial Manager's population:",
                        min = 100,
                        max = 500,
                        value = 100,
                        step = 100)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           tableOutput("table1")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
      players <- 2 #Two players are the manager and neighbor
      N.years <- 20
      
      hM <- seq(0,1,by = 0.5) 
      hN <- seq(0,1,by = 0.5)
      
      h <- expand.grid(hM,hN) #creating combinations of harvest rates (alternatives)
      colnames(h) <- c("M", "N")
      N.h <- length(h$M) #total number of management alternatives
      
      X0 <- c(input$M1,input$N1) #initial population in manager, neighbor area respectively
      X <- X.rem <- C<- G <- R <- array(0, dim = c(players, N.years, N.h)) #arrays
      X[,1,] <- X0 #population at the first time step 
      
      #Rates for biological model
      r_M <- input$r #growth rate manager
      r_N <- input$r #growth rate neighbor
      K_M <- input$K #carrying capacity manager
      K_N <- input$K #carrying capacity neighbor
      b <- input$b #dispersal rate into other location
      
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
      
    
      })
    

}

# Run the application 
shinyApp(ui = ui, server = server)
