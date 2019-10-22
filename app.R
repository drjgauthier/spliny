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
      actionButton("do", "Simulate and plot probabilities"),
      helpText("Re-click to launch another simulation")
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
  biomarker_m_r <- reactive({
    cut2(biomarker_r(),g=2)
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
      biomarker_r <- biomarker_r()
      biomarker_m <- cut2(biomarker_r,g=2)
      data.frame(biomarker_r=biomarker_r(),biomarker_m,event=ytest_r(),trueprobability=prob_r())
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
    biomarker_m_r <- biomarker_m_r()
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
    f <- lrm(ytest_r() ~ biomarker_r)
    f_AIC <- data.frame(round(AIC(f),digits=0))
    stats <- data.frame(f$stats["Model L.R."],f$stats["C"],f$stats["Brier"],f_AIC)
    colnames(stats) <-c("Likelihood ratio","C-index","Brier Score","AIC")
    
    f_rcs <- lrm(ytest_r() ~ rcs(biomarker_r,quantile(biomarker_r,knots_s)))  
    f_rcs_AIC <- data.frame(round(AIC(f_rcs),digits=0))
    stats_rcs <- data.frame(f_rcs$stats["Model L.R."],f_rcs$stats["C"],f_rcs$stats["Brier"],f_rcs_AIC)
    colnames(stats_rcs) <-c("Likelihood ratio","C-index","Brier Score","AIC")
    
    f_m <- lrm(ytest_r() ~ biomarker_m_r)
    f_m_AIC <- data.frame(round(AIC(f_m),digits=0))
    stats_m <- data.frame(f_m$stats["Model L.R."],f_m$stats["C"],f_m$stats["Brier"],f_m_AIC)
    colnames(stats_m) <-c("Likelihood ratio","C-index","Brier Score","AIC")
    
    Model<-rbind("Biomarker dichotomized", "Biomarker modeled as continuous variable","Biomarker modeled with spline")
    
    AIC_t <- bind_rows(stats_m,stats,stats_rcs)
    stats_df <- data.frame(Model,AIC_t)
    colnames(stats_df) <-c("Model","Likelihood ratio","C-index","Brier Score","AIC")
    stats_df
  })
  output$ggplot <- renderPlot({
    df_r <- df_r()
    knots_s <<- switch(input$knots, 
                      "3" = c(.10,.5,.90),
                      "4" = c(.05,.35,.65,.95),
                      "5" = c(.05,.275,.5,.725,.95),
                      "6" = c(.05,.23,.41,.59,.77,.95),
                      "7" = c(.025,.1833,.3417,.5,.6583,.8167,.975)
    )
    #Fit with dichotomized model
    pred_m<-data.frame(predict(glm(ytest_r() ~ biomarker_m_r(),family="binomial"),se.fit=TRUE,type="response"))
    #Fit without rcs
    pred_c<-data.frame(predict(glm(ytest_r() ~ biomarker_r(),family="binomial"),se.fit=TRUE,type="response"))
    #Fit with rcs
    pred_s<-data.frame(predict(glm(ytest_r() ~ rcs(biomarker_r(),quantile(biomarker_r(),knots_s))),se.fit=TRUE,type="response"))
    df2 <- data.frame(biomarker_r(),prob_r(),pred_m[1:2],pred_c[1:2],pred_s[1:2])
    colnames(df2) <- c("biomarker_r","simulated probability","prob_m","prob_se_m","prob_c","prob_se_c",
                       "prob_s","prob_se_s")
    df2$prob_m_upper <- df2$prob_m + 1.96*df2$prob_se_m
    df2$prob_m_lower <- df2$prob_m - 1.96*df2$prob_se_m
    df2$prob_c_upper <- df2$prob_c + 1.96*df2$prob_se_c
    df2$prob_c_lower <- df2$prob_c - 1.96*df2$prob_se_c
    df2$prob_s_upper <- df2$prob_s + 1.96*df2$prob_se_s
    df2$prob_s_lower <- df2$prob_s - 1.96*df2$prob_se_s
    df2 %>% select(biomarker_r,
                   "simulated probability",
                   prob_m,
                   prob_m_upper,
                   prob_m_lower,
                   prob_c,
                   prob_c_upper,
                   prob_c_lower,
                   prob_s,
                   prob_s_upper,
                   prob_s_lower) %>% 
      gather(key=prob_type,value=prediction,-c("biomarker_r",
                                               "prob_m_upper",
                                               "prob_m_lower",
                                               "prob_c_upper",
                                               "prob_c_lower",
                                               "prob_s_upper",
                                               "prob_s_lower")) -> df2_g
    df2_g$prob_type <- as.factor(df2_g$prob_type)
    df2_g$prob_type = factor(df2_g$prob_type,levels(df2_g$prob_type)[c(4,2,1,3)])
    p <- df2_g %>% 
      ggplot() + 
      geom_line(data=df2_g,aes(x=biomarker_r,y=prediction,col=prob_type),alpha=0.8,size=1)+
      labs(x="Biomarker serum concentration",y="Probability of event")+
      scale_color_manual(values=c("red","blue","black","purple"),
                         labels=c("Simulated probabilities","Biomarker dichotomized","Biomarker continuous","Biomarker with spline"))+
      theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),axis.text.x = element_text(size=16,margin = margin(t = 0, r = 0, b = 10, l = 0)),axis.text.y = element_text(size=16),legend.title = element_blank(),legend.text = element_text(size=16))
    p <- p + geom_point(data=subset(df2_g,prob_type=="simulated probability"),aes(x=biomarker_r,y=prediction),col="red",alpha=0.1)
    p <- p + geom_ribbon(data=df2_g,aes(x=biomarker_r,ymin=prob_c_lower,ymax=prob_c_upper),fill="blue",alpha=0.08)
    p <- p + geom_ribbon(data=df2_g,aes(x=biomarker_r,ymin=prob_m_lower,ymax=prob_m_upper),fill="black",alpha=0.08)
    p + geom_ribbon(data=df2_g,aes(x=biomarker_r,ymin=prob_s_lower,ymax=prob_s_upper),fill="purple",alpha=0.08)
  })
  })
})
shinyApp(ui, server)
