library(shiny)
library(shinydashboard)
library(ggplot2)
library(rms)
library(splines)
library(huxtable)
library(DT)
ui <- fluidPage(
  titlePanel("Playing with restricted cubic splines"),
  sidebarLayout(
    sidebarPanel(
      h4("Tweak relationship with outcome:"),
      sliderInput("intercept", "Intercept:", min = -5, max = 5, value = -1),
      sliderInput("beta", "Beta:", min = -0.05, max = 0.05, value = 0.005),
      h4("Tweak simulation parameters:"),
      sliderInput("n", "Subject number to simulate:", min = 0, max = 10000, value = 1000),
      h4("Tweak simulated biomarker:"),
      sliderInput("mean", "Biomarker mean:", min = 0, max = 10000, value = 500),
      sliderInput("sd", "Biomarker standard deviation", min = 0, max = 10000, value = 250),
      selectInput(inputId = "knots", label = "Select number of knots", choices = c("knots3","knots4","knots5","knots6","knots7"),
                  selected = "knots3")
    ),
    mainPanel(plotOutput("ggplot"),
              plotOutput("plot2"),
               plotOutput("plot"),
               dataTableOutput("dataframe")
             
    )))


server <- shinyServer(function(input,output){
  biomarker_r <- reactive({
    round(rnorm(input$n,input$mean,input$sd),digits=2)
  })
  biomarkerlims<- reactive({
    range(biomarker_r())
  }) 
  #Generating Test Data
  
  biomarker.grid<- reactive({
    seq(from=biomarkerlims()[1], to = biomarkerlims()[2])
  })
  #true logit values for the xtest variable
  logit_l_r <- reactive({
    input$intercept + (biomarker_r() * input$beta)
  })
  #true probabilities
  # nonlinpred <- intercept + (xtest * beta) + (0.00001*xtest^2)+ (0.00000001*xtest^3)
  prob_l_r <- reactive({
    round(exp(logit_l_r())/(1 + exp(logit_l_r())),digits=2)
  })
  #prob_nl <- exp(nonlinpred)/(1 + exp(nonlinpred))
  # Simulate binary y to have Prob(y=1) = prob(linpred)]
  runis_r <- reactive({
    runis <- runif(input$n,0,1)
  })
  ytest_r <- reactive({
    ytest <- ifelse(runis_r() < prob_l_r(),1,0)
  })
  df_r <- reactive({
    data.frame(biomarker=biomarker_r(),event=ytest_r(),probability=prob_l_r())
  }) 
  output$dataframe = DT::renderDataTable({
    df_r()
  })
  output$plot<- renderPlot({
    hist(biomarker_r(),col='lightblue')
  })
  #Define knots locations
  knots3_r <- reactive({
  q3<-data.frame(quantile(biomarker_r(),probs=c(.10,.5,.90)))
  knots3 = c(q3[1,1],q3[2,1],q3[3,1])
  })
  
  knots4_r <- reactive({
  q4<-data.frame(quantile(biomarker_r(),probs=c(.05,.35,.65,.95)))
  knots4 = c(q4[1,1],q4[2,1],q4[3,1],q4[4,1])
  })
  
  knots5_r <- reactive({
  q5<-data.frame(quantile(biomarker_r(),probs=c(.05,.275,.5,.725,.95)))
  knots5 = c(q5[1,1],q5[2,1],q5[3,1],q5[4,1],q5[5,1])
  })
  
  knots6_r <- reactive({
  q6<-data.frame(quantile(biomarker_r(),probs=c(.05,.23,.41,.59,.77,.95)))
  knots6 = c(q6[1,1],q6[2,1],q6[3,1],q6[4,1],q6[5,1],q6[6,1])
  })
  
  knots7_r <- reactive({
  q7<-data.frame(quantile(biomarker_r(),probs=c(.025,.1833,.3417,.5,.6583,.8167,.975)))
  knots7 = c(q7[1,1],q7[2,1],q7[3,1],q7[4,1],q7[5,1],q7[6,1],q7[7,1])
  })
  
  output$plot2 <- renderPlot({
    biomarker_r <- biomarker_r()
    knots_s <- switch(input$knots, 
                      "knots3" = knots3_r(),
                      "knots4" = knots4_r(),
                      "knots5" = knots5_r(),
                      "knots6" = knots6_r(),
                      "knots7" = knots7_r()
    )
    fit3 <- glm(ytest_r() ~ bs(biomarker_r,knots = knots3_r()))
    fit4 <- glm(ytest_r() ~ bs(biomarker_r,knots = knots4_r()))
    fit5 <- glm(ytest_r() ~ bs(biomarker_r,knots = knots5_r()))
    fit6 <- glm(ytest_r() ~ bs(biomarker_r,knots = knots6_r()))
    fit7 <- glm(ytest_r() ~ bs(biomarker_r,knots = knots7_r()))
    
    plot(biomarker_r,ytest_r(),col="grey",xlab="Biomarker",ylab="Probability of Response")
    points(biomarker.grid(),predict(fit3,newdata = list(biomarker_r=biomarker.grid())),col="darkgreen",lwd=2,type="l") 
    points(biomarker.grid(),predict(fit4,newdata = list(biomarker_r=biomarker.grid())),col="darkred",lwd=2,type="l")
    points(biomarker.grid(),predict(fit5,newdata = list(biomarker_r=biomarker.grid())),col="darkblue",lwd=2,type="l")
    points(biomarker.grid(),predict(fit6,newdata = list(biomarker_r=biomarker.grid())),col="yellow",lwd=2,type="l")
    points(biomarker.grid(),predict(fit7,newdata = list(biomarker_r=biomarker.grid())),col="darkorange",lwd=2,type="l")
    
    abline(v=knots_s,lty=2,col="darkgreen")
  })
  output$ggplot <- renderPlot({
    biomarker_r <- biomarker_r()
    biomarker.grid <- biomarker.grid()
        knots_s <- switch(input$knots, 
                      "knots3" = knots3_r(),
                      "knots4" = knots4_r(),
                      "knots5" = knots5_r(),
                      "knots6" = knots6_r(),
                      "knots7" = knots7_r()
    )
    fit3 <- glm(ytest_r() ~ bs(biomarker_r,knots = knots3_r()))    
    #newdat = expand.grid(biomarker_r = seq(min(biomarker_r), max(biomarker_r), by = 1))
    newdat = data.frame(predict(fit3, newdata = biomarker.grid,se.fit=TRUE))
    newdat$upper <- newdat$fit + 1.96*newdat$se.fit
    newdat$lower <- newdat$fit - 1.96*newdat$se.fit
    newdat <-cbind(newdat,biomarker.grid)
    ggplot(df_r(), aes(x = biomarker_r, y = ytest_r())) +
      geom_point(alpha = .5) +
      geom_line(data = newdat, aes_string(x=biomarker,y = fit), size = 1)
  })
  })
shinyApp(ui, server)

