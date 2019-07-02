library(shiny)
library(shinydashboard)
library(ggplot2)
library(rms)
library(splines)
library(huxtable)
library(DT)
library(plotly)
theme_set(theme_minimal())
ui <- fluidPage(
  titlePanel("Playing with restricted cubic splines"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "knots", label = "Select number of knots", choices = c("3","4","5","6","7"),
                  selected = "knots3"),
      h4("Simulate non-linear relationship with outcome:"),
      sliderInput("intercept", "Intercept:", min = -5, max = 5, value = -1),
      sliderInput("beta", "Beta:", min = -0.05, max = 0.05, value = 0.005),
      sliderInput("a", "a:", min = 0.0000001, max = 0.0001, value = 0.000001),
      sliderInput("b", "b:", min = 1e-10, max = 1e-7, value = 0.00000001),
            h4("Simulation parameters:"),
      sliderInput("n", "Number of subjects:", min = 0, max = 10000, value = 1000),
      sliderInput("mean", "Biomarker mean:", min = 0, max = 10000, value = 500),
      sliderInput("sd", "Biomarker standard deviation", min = 0, max = 5000, value = 150)
    ),
    mainPanel(plotlyOutput("ggplot"),
              tableOutput("AIC_t")
              )))
server <- shinyServer(function(input,output){
  biomarker_r <- reactive({
    round(rnorm(input$n,input$mean,input$sd),digits=2)
  })
  #true logit values for the xtest variable
  logit_r <- reactive({
    input$intercept + (biomarker_r() * input$beta) + (input$a*biomarker_r()^2)+ (0.00000001*biomarker_r()^3)
  })
  #true probabilities
  prob_r <- reactive({
    exp(logit_r())/(1 + exp(logit_r()))
  })
  # Simulate binary y to have Prob(y=1) = prob(linpred)]
  runis_r <- reactive({
    runis <- runif(input$n,0,1)
  })
  ytest_r <- reactive({
    ytest <- ifelse(runis_r() < prob_r(),1,0)
  })
  df_r <- reactive({
    data.frame(biomarker=biomarker_r(),event=ytest_r(),trueprobability=prob_r())
  }) 
  output$AIC_t <- renderTable({
    biomarker_r <- biomarker_r()
    df_r <- df_r()
    ddist <<- datadist(df_r)
    options(datadist='ddist')
    knots_s <<- switch(input$knots, 
                       "knots3" = c(.10,.5,.90),
                       "knots4" = c(.05,.35,.65,.95),
                       "knots5" = c(.05,.275,.5,.725,.95),
                       "knots6" = c(.05,.23,.41,.59,.77,.95),
                       "knots7" = c(.025,.1833,.3417,.5,.6583,.8167,.975)
    )
    f <- lrm(ytest_r() ~ (biomarker_r))
    f_AIC <- data.frame(round(AIC(f),digits=0))
    f_rcs <- lrm(ytest_r() ~ rcs(biomarker_r,quantile(biomarker_r,knots_s)))  
    f_rcs_AIC <- data.frame(round(AIC(f_rcs),digits=0))
    AIC_t <- cbind(f_AIC,f_rcs_AIC)
    names(AIC_t) <- c("AIC without spline","AIC with spline")
    AIC_t
  })
  output$ggplot <- renderPlotly({
    biomarker_r <- biomarker_r()
    df_r <- df_r()
    ddist <<- datadist(df_r)
    options(datadist='ddist')
    knots_s <<- switch(input$knots, 
                      "3" = c(.10,.5,.90),
                      "4" = c(.05,.35,.65,.95),
                      "5" = c(.05,.275,.5,.725,.95),
                      "6" = c(.05,.23,.41,.59,.77,.95),
                      "7" = c(.025,.1833,.3417,.5,.6583,.8167,.975)
    )
    f <- lrm(ytest_r() ~ (biomarker_r))
    pred <- data.frame(Predict(f,biomarker_r=seq(10,1000,by=1),fun=plogis))
    f_rcs <- lrm(ytest_r() ~ rcs(biomarker_r,quantile(biomarker_r,knots_s)))  
    pred_rcs <- data.frame(Predict(f_rcs,biomarker_r=seq(10,1000,by=1),fun=plogis))
    p <- ggplot() + geom_point(data=df_r(),aes(x=biomarker,y=trueprobability),alpha=0.3,col="red") + annotate("text", x = 800, y = 0.4, label = "True probabilities",col="red")
    p <- p + geom_line(data=pred,aes(x=biomarker_r,y=yhat)) + annotate("text", x = 800, y = 0.3, label = "Predictions without spline")
    p<- p + geom_line(data=pred_rcs,aes(x=biomarker_r,y=yhat),col="blue") + annotate("text", x = 800, y = 0.2, label = "Predictions with spline",col="blue")
    p<-p + labs(x="Biomarker serum concentration",y="Probability of event") +theme(axis.title.x = element_text(size=12),
                                                                                axis.title.y = element_text(size=12),
                                                                                axis.text.x = element_text(size=10),
                                                                                axis.text.y = element_text(size=10))
    ggplotly(p)
  })
  })
shinyApp(ui, server)
