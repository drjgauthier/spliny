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
      selectInput("functiontype","Function modeling the non-linear relationship",
                  c("Hyperbolic tangent","Tangent","Sine")),
      conditionalPanel(
        condition="input.functiontype=='Tangent'",
        sliderInput("t","Tweak function shape (a)", min = 5, max = 20, value = 5),
        helpText("Increasing this parameter steepens the shape of the relationship")
              ),
      conditionalPanel(
        condition="input.functiontype=='Sine'",
        sliderInput("s","Tweak function shape (a)", min = 100, max = 500, value = 150)
        ),
      conditionalPanel(
        condition="input.functiontype=='Hyperbolic tangent'",
        sliderInput("ht","Tweak function shape (a)", min = 50, max = 400, value = 50),
        helpText("Increasing this parameter flattens the shape of the relationship")
        ),
      sliderInput("n", "Number of simulated subjects", min = 0, max = 1000, value = 500),
      selectInput(inputId = "knots", label = "Select number of knots for the restricted cubic spline", choices = c("3","4","5","6","7"),
                  selected = "5"),
      actionButton("do", "Simulate and plot probabilities")
    ),
    mainPanel(plotOutput("ggplot",width = "100%",height = "500px",hover=TRUE),
              helpText("Shaded areas represent the 95% confidence bands of the predictions"),
              hr(),
              tableOutput("AIC_t"),
              textOutput("AIC_TEST"),
              helpText("NB: a lower AIC (Akaike Information Criterion) value indicates a better fit"),
              hr(),
              conditionalPanel(
                condition="input.functiontype=='Hyperbolic tangent'",
                withMathJax(),
                helpText("The relationship between the biomarker and the 'true' simulated logit of the event is modeled using the hyperbolic tangent function: $$logit(biomarker) = \\frac{\\tanh(biomarker-500)}a$$")),
              conditionalPanel(
                condition="input.functiontype=='Tangent'",
                withMathJax(),
                helpText("The relationship between the biomarker and the 'true' simulated logit of the event is modeled using the trigonometric tangent function: $$logit(biomarker) = \\frac{\\frac{\\tan(biomarker-500)}{325}}a$$")),
              conditionalPanel(
                condition="input.functiontype=='Sine'",
                withMathJax(),
                helpText("The relationship between the biomarker and the 'true' simulated logit of the event is modeled using the trigonometric sine function: $$logit(biomarker) = \\frac{\\sin(biomarker)}a$$")),
              hr(),
              a(href="https://github.com/drjgauthier/spliny/blob/master/app.R","App R code on Github"),
              a(href="https://www.rdocumentation.org/packages/rms/versions/5.1-3.1","R documentation for the rms package")
              )))
server <- shinyServer(function(input,output){
  observeEvent(input$do,{
  biomarker_r <- reactive({
    round(runif(input$n,min=0,max=1000),digits=0)
  })
    #true logit values for the xtest variable
    logit_r <- reactive({
      if(input$functiontype=="Hyperbolic tangent") {
        tanh((biomarker_r()-500)/input$ht)}
      else 
        if(input$functiontype=="Sine") {
          sin(biomarker_r()/input$s)
        }
      else
        tan((biomarker_r()-500)/325)/input$t  
    })
    #simulated 'true' probabilities
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
           "3" = c(.10,.5,.90),
           "4" = c(.05,.35,.65,.95),
           "5" = c(.05,.275,.5,.725,.95),
           "6" = c(.05,.23,.41,.59,.77,.95),
           "7" = c(.025,.1833,.3417,.5,.6583,.8167,.975)
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
      #xlim(0,1000)+
      labs(x="Biomarker serum concentration",y="Probability of event")+
      scale_color_manual(values=c("red","blue","black"),labels=c("Predicted probabilities without spline","Predicted probabilities with spline","True probabilities"))+
      geom_point(data=subset(df_g,prob_type=="trueprobability"),aes(x=biomarker,y=prediction),alpha=0.1,size=4)+
      theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),axis.text.x = element_text(size=16,margin = margin(t = 0, r = 0, b = 10, l = 0)),axis.text.y = element_text(size=16),legend.title = element_blank(),legend.text = element_text(size=16))
    p <- p + geom_ribbon(data=df,aes(x=biomarker,ymin=lower,ymax=upper),fill="red",alpha=0.08)
    p + geom_ribbon(data=df,aes(x=biomarker,ymin=lower_s,ymax=upper_s),fill="blue",alpha=0.08)
  })
  })
})
shinyApp(ui, server)
