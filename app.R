library(shiny)
library(ggplot2)
library(rms)
library(plotly)
library(dplyr)
library(tidyr)
theme_set(theme_minimal())
ui <- fluidPage(
  titlePanel("Spliny: restricted cubic splines to model non-linear effects using logistic regression"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("intercept", "Intercept (baseline risk)", min = -2, max = 2, value = 2),
      sliderInput("beta","Coefficient 1 (α1)", min = -0.01, max = 0.002, value = -0.01),
      helpText("Tweaks shape of the relationship between the biomarker and the outcome variable"),
      sliderInput("a", "Coefficient 2 (α2)", min = -6, max = -5, value = -5),
      helpText("Tweaks shape of the relationship between the biomarker and the outcome variable"),
      sliderInput("n", "Number of subjects to simulate", min = 0, max = 1000, value = 500),
      sliderInput("sd", "Biomarker standard deviation (spread)", min = 10, max = 300, value = 200),
      selectInput(inputId = "knots", label = "Select number of knots for the restricted cubic spline", choices = c("3","4","5","6","7"),
                  selected = "knots3")
    ),
    mainPanel(plotOutput("ggplot",width = "100%",height = "500px",hover=TRUE),
              helpText("Shaded areas represent the 95% confidence bands of the predictions"),
              hr(),
              tableOutput("AIC_t"),
              textOutput("AIC_TEST"),
              helpText("NB: a lower AIC (Akaike Information Criterion) value indicates a better fit"),
              hr(),
              withMathJax(),
              helpText("The relationship between the biomarker (x) and the 'true' logit of the event is modeled using the following quadratic function:$$logit(x) = intercept + \\alpha_1 x^2 + 10^{\\alpha_2}x^2$$"),
              hr(),
              a(href="https://github.com/drjgauthier/splines/blob/master/app.R","R code")
              )))
server <- shinyServer(function(input,output){
  biomarker_r <- reactive({
    round(rnorm(input$n,500,input$sd),digits=0)
  })
  #true logit values for the xtest variable
  logit_r <- reactive({
    input$intercept + (biomarker_r() * input$beta) + ((10^input$a)*(biomarker_r()^2))
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
  
  knots_r <<- reactive({
    switch(input$knots, 
           "knots3" = c(.10,.5,.90),
           "knots4" = c(.05,.35,.65,.95),
           "knots5" = c(.05,.275,.5,.725,.95),
           "knots6" = c(.05,.23,.41,.59,.77,.95),
           "knots7" = c(.025,.1833,.3417,.5,.6583,.8167,.975)
    )
  })
  
  output$AIC_t <- renderTable({
    biomarker_r <- biomarker_r()
    df_r <- df_r()
    knots_r <- knots_r()
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
    f_AIC <- data.frame(round(AIC(f),digits=0))
    f_rcs <- lrm(ytest_r() ~ rcs(biomarker_r,quantile(biomarker_r,knots_s)))  
    f_rcs_AIC <- data.frame(round(AIC(f_rcs),digits=0))
    AIC_t <- cbind(f_AIC,f_rcs_AIC)
    names(AIC_t) <- c("AIC without spline","AIC with spline")
    AIC_t
    })
  output$ggplot <- renderPlot({
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
    #Fit without rcs
    f <- lrm(ytest_r() ~ biomarker_r)
    pred <- Predict(f,biomarker_r=seq(1,1000,by=1),fun=plogis,np = 1000)
    names(pred) <- c("biomarker","predicted_prob","upper","lower")
    #Fit with rcs
    f_rcs <- lrm(ytest_r() ~ rcs(biomarker_r,quantile(biomarker_r,knots_s)))  
    pred_rcs <- Predict(f_rcs,biomarker_r=seq(1,1000,by=1),fun=plogis,np = 1000)
    names(pred_rcs) <- c("biomarker","predicted_prob_s","upper_s","lower_s")
    #Data wrangling
    df <- left_join(df_r,pred,by="biomarker")
    df <- left_join(df,pred_rcs,by="biomarker")
    gather(df,key=prob_type,value=prediction,-biomarker,-event,-upper,-lower,-upper_s,-lower_s) -> df_g
    #Graph "true" probabilities and predictions
    p <- ggplot() + 
      geom_line(data=df_g,aes(x=biomarker,y=prediction,col=prob_type),alpha=0.8,size=1)+
      xlim(0,1000)+
      labs(x="Biomarker serum concentration",y="Probability of event")+
      scale_color_manual(values=c("red","blue","black"),labels=c("Predicted probabilities without spline","Predicted probabilities with spline","True probabilities"))+
      geom_point(data=subset(df_g,prob_type=="trueprobability"),aes(x=biomarker,y=prediction),alpha=0.1,size=4)+
      theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),axis.text.x = element_text(size=16,margin = margin(t = 0, r = 0, b = 10, l = 0)),axis.text.y = element_text(size=16),legend.title = element_blank(),legend.text = element_text(size=16))
    p <- p + geom_ribbon(data=df,aes(x=biomarker,ymin=lower,ymax=upper),fill="red",alpha=0.08)
    p + geom_ribbon(data=df,aes(x=biomarker,ymin=lower_s,ymax=upper_s),fill="blue",alpha=0.08)
  })
  output$math <- renderUI({
    withMathJax(helpText("The relationship between the biomarker (x) and the 'true' logit of the event is modeled using the following quadratic function: logit(x) = intercept + $$\\alpha1$$x + 10^α2*x^2) $$\\alpha1^2$$"))
  })
  })
shinyApp(ui, server)
