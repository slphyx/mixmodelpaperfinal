library(shiny)
library(pROC)
library(mixtools)
library(MASS)
library(shinythemes)
library(shinydashboard)

ui <- fluidPage(
  theme = shinytheme("journal"),
  h2("Identify  artemisinin resistance from parasite clearance half-life data"),
  
  tabsetPanel(
    id="panels",
    tabPanel(title = "Introduction",
             br(),
             p("In the World Health Organization's ", a(href="http://www.who.int/malaria/publications/atoz/update-artemisinin-resistance-april2016/en/", 
                                                        "Update on artemisinin and ACT resistance - April 2016"),
               "the cut-off value of greater than 10% of patients with a half-life of the parasite clearance slope more than 5 hours after treatment with ACT or 
               artesunate monotherapy is used as one of the definitions of \"suspected endemic artemisinin resistance\"."),
            p("In the following two examples, the cut-off value will either miss or overestimate the artemisinin resistance.",
               "These examples assume that the parasite clearance half-lives are in log-normal distribution, and that the values for 
                sensitive and resistant populations each assume unimodal distribution."
             ),
             #p("Reactive buttons to change between Example 1 and Example 2"),
             p("Click on each button to populate the respective example.",
               "You can also use the parameters to change the outputs."),
             actionButton("eg1","Example 1"),
             actionButton("eg2","Example 2"),
             hr(),
             h4("Example"),
             p("In a sample of ",
               strong(textOutput("nnO",inline=T)),
               "individuals, when the sensitive population has a geomatric half-life mean of ", 
               strong(textOutput("senmuO", inline=T)), 
               "hours with the standard deviation of ",
               strong(textOutput("sensdO",inline=T)),
               "hours and the resistant population has a geomatric half-life mean of ",
               strong(textOutput("resmuO",inline=T)),
               "hours with the standard deviation of ",
               strong(textOutput("ressdO",inline=T)),
               "hours then the distribution will look like the following plot. If the percentage of 
resistant population is",
               strong(textOutput("prop_resistO",inline=T)),
               
               ", the cutoff value of ",
               strong(textOutput("cutoffO",inline=T)),
               "hours will detect ",
               strong(textOutput("capturedO", inline=T)),
               "percent of the population as resistant.", 
               "The use of cut-off value can be augmented by additional information.",#"The cut-off value then has to be adjusted.",
               br(),
               p("Here, we're providing a tool based on the model developed by", a(href="http://bit.ly/White-et-al-2015","White et al.(2015)"),
                "to anlayze the parasite clearance half-life data as 
               distributions of artemisinin-sensitive and artemisinin-resistant populations. 
               Please go to the", actionLink("link_to_MMpage", "next page"), "to use it.")
               #HTML("<a href='#histoplot1'>next page</a>."),
               
               #br(),
               #p("You can also use the parameters below to change the following plots.")
               #  "In fact ",
               # textOutput("overlayO", inline=T)
             ),
             
             fluidRow(
               column(5,
                      #plotOutput(outputId = "densityplot")
                      plotOutput(outputId = "histoplot1")
               ),
               column(5,
                      plotOutput(outputId = "ROC")
               )
               
             ),
             wellPanel(
               fluidRow(
                 column(4,
                        h3("Sensitive Distribution"),
                        sliderInput(inputId = "senmu",
                                    label = "Mean half-life (hours)",
                                    value = 3, min = 1, max = 6.5, step = .5
                        ),
                        sliderInput(inputId = "sensd",
                                    label = "SD (hours)",
                                    value = 1.26, min = 1, max = 2.1
                        )
                 ),
                 column(4,
                        h3("Resistant Distribution"),
                        sliderInput(inputId = "resmu",
                                    label = "Mean half-life (hours)",
                                    value = 6, min = 5, max = 10, step = .5
                        ),
                        sliderInput(inputId = "ressd",
                                    label = "SD (hours)",
                                    value = 1.22, min = 1, max = 2.1
                        ),
                        sliderInput(inputId = "prop_resist",
                                    label = "% of resistant population",
                                    value = 10, min = 0, max = 100
                        )
                 ),
                 column(4,
                        numericInput(inputId = "nn",
                                     label = "Sample Size:",
                                     value = 500
                        ),
                        checkboxInput(inputId = "showcutoff",
                                      label = "Show cutoff line in the histogram",
                                      value = TRUE
                        ),
                        sliderInput(inputId = "cutoff",
                                    label = "Cut-off half-life value",
                                    value = 5, min = 0, max = 10, step=.5
                        )
                 )
               )
             ),
             br() #,
             #p("This is the end of part 1.")
             
             ###test####
             #p(textOutput("genDataOut"))
    ),
    tabPanel(title = "Use the Mixture Model from White et al.",
             ############################
             ###Portions from MixModel###
             ############################
             br(),
             fluidRow(column(6,
                             h4("Before you start"),
                             p("Before using your own data to run the model, it might be worth checking out how
                               your data input should be like, since", strong("the model will not run if the data is not in the 
                                                                              correct format."), "The data input has to be a ", strong("csv"), "file with a", strong("single
                                                                                                                                                                     column of half-life clearance data."),"There must be",strong("no column names or no row names.")),  
                             #br(),
                             # "You can download our simulated default dataset and have a look.",
                             # "The following is the result of the model run using the default simulated dataset. 
                             #            You can download the default dataset here:",
                             downloadButton("defaultData", "Download default/template dataset"),
                             p("You can use our simulated dataset as a template to copy and paste (overwrite) the values of your data.
                               When saving, just keep the", strong("csv"), "format. Your uploaded data is used only for running the model, 
                               and it will not be stored.")
                             ),
                      column(5,
                             h4("Using your data"),
                             wellPanel(
                               fileInput(inputId = "file", label = "Your input file: ", accept = c(".csv"))
                             ),
                             #something about the default dataset
                             verbatimTextOutput("aboutDefault") #textOutput("aboutDefault")
                             )
                      ),
             hr(),
             #verbatimTextOutput("explain"),
             h4("Results"),
             textOutput("explain"),
             fluidRow(
               column(5,
                      plotOutput(outputId = "histoplot2")
               ),
               column(5,
                      plotOutput(outputId = "ROC2")
               )
               
             ),
             wellPanel(
               fluidRow(
                 column(4,
                        checkboxInput(inputId = "showcutoff2",
                                      label = "Show cutoff line in the histogram",
                                      value = FALSE
                        ),
                        sliderInput(inputId = "cutoff2",
                                    label = "Cut-off half-life value",
                                    value = 5, min = 0, max = 10, step=.5
                        )
                 )#,
                 # column(4,
                 #        ###other downloads maybe
                 # )
               )
             ),
             h4("Downloads"),
             downloadButton('downloadhistoplot2',"Download the histogram"),
             downloadButton('resultData',"Download the results in a table")
    ),
    tabPanel(title="Limitations",
             br(),
             h4("All the assumptions and limitations from the model of", a(href="http://bit.ly/White-et-al-2015","White et al.(2015)"),"are applied here."),
             tags$ul(tags$li("The clearance half-lives of infections with a particular sensitivity are assumed to follow unimodal distributions of log-normal type."),
                     tags$li("The maximum number of subpopulations the model can detect is 5."),
                     tags$li("As described in the", a(href="http://bit.ly/White-2015-S1","Supporting information 1 of White et al. (2015)"),", the model's ability to differentiate between subpopulations depends on means and standard deviations of the component distributions, sample size, and number of subpopulations. For instance, from a sample size of 50, the model will be able to differentiate between subpopulations of geometric mean half-lives with a difference of 3 or more hours. From a sample size of 1000, the model will be able to differentiate subpopulations whose geometric mean half-lives differ by only 0.5 hours. The model's prediction will also decrease with the increase in the true number of subpopulations. Eg., For a sample size of 1,000, the model will correctly predict 96%, 91%, 70%, 46% and 21% for the input mixture distributions of 1, 2, 3, 4 and 5 components respectively."),
                     tags$li("While using this web application, when the window of the browser is resized, the histogram will disappear. They will reappear when you change one of the parameters given for the histogram.")
             )),
    tabPanel(title="Related Resources",
             br(),
             tags$ul(
               tags$li("The original paper by White et al. (2015) on which this web application is based:", a(href="http://bit.ly/White-et-al-2015","Defining the In Vivo Phenotype of Artemisinin-Resistant Falciparum Malaria: A Modelling Approach")),
               tags$li("Exploration of the model's limitation:", a(href="http://bit.ly/White-2015-S1", "Supporting information 1, White et al. (2015)" )),
               tags$li(a(href="http://bit.ly/White-2015-code","Source codes of the original model by White et al. (2015)")),
               tags$li(a(href="http://bit.ly/2cJWpVr","WHO updates on artemisinin resistance")),
               tags$li(a(href="http://bit.ly/2d4hV7V","Parasite Clearance Estimator (PCE)")),
               tags$li("Published paper for this web application: submitted"),
               tags$li(a(href="http://bit.ly/White-2015-shiny-code","Source codes for this web application"))
             )
  )
)
)

server <- function(input, output, session) {
  #having the parameters reflect on the text description
  set.seed(3)
  rlnormR <- repeatable(rlnorm)
  dlnormR <- repeatable(dlnorm)
  
  output$senmuO <- renderText({
    as.character(input$senmu)
  })
  output$sensdO <- renderText({
    as.character(input$sensd)
  })
  output$resmuO <- renderText({
    as.character(input$resmu)
  })
  output$ressdO <- renderText({
    as.character(input$ressd)
  })
  output$prop_resistO <- renderText({
    as.character((input$prop_resist))
  })
  output$nnO <- renderText({
    as.character(input$nn)
  })
  output$cutoffO <- renderText({
    as.character(input$cutoff)
  })
  
  ####updateSliderInput####
  observeEvent(input$eg1,{
    updateSliderInput(session, "senmu", value=3)
    updateSliderInput(session, "sensd", value=1.22)
    updateSliderInput(session, "resmu", value=6)
    updateSliderInput(session, "ressd", value=1.22)
    updateSliderInput(session, "prop_resist", value=11)
    updateNumericInput(session, "nn", value=500)
    updateSliderInput(session, "cutoff", value=5)
  })
  observeEvent(input$eg2,{
    updateSliderInput(session, "senmu", value=3.5)
    updateSliderInput(session, "sensd", value=1.35)
    updateSliderInput(session, "resmu", value=6)
    updateSliderInput(session, "ressd", value=1.22)
    updateSliderInput(session, "prop_resist", value=2)
    updateNumericInput(session, "nn", value=1000)
    updateSliderInput(session, "cutoff", value=5)
  })
  
  #link to next page
  observeEvent(input$link_to_MMpage,{
    newvalue <- "Use the Mixture Model from White et al."
    updateTabItems(session, "panels", newvalue)
  })
  
  #something about the default dataset
  rvAboutDataset <- reactiveValues(text = "**Processing status** \nThe following is the output of the model \nusing the default dataset. This text will \nfade when the model is running.")
# The output will change once you've uploaded your data and the model approximation is completed. 
#                                    Note: After your data input if you are still seeing this message, 
#                                    the model is running in the background. It'll take sometime depending on the size of your data.")
  output$aboutDefault <- renderText(rvAboutDataset$text)
  observeEvent(input$file,{rvAboutDataset$text="*************************************\n* Model approximation is completed! *\n*************************************\nSee below for the results!"})
  
  rvExplain <- reactiveValues(text="")
  output$explain <- renderText(rvExplain$text)
  
  ####the functions for plotting####
  senmuR <- reactive({log(input$senmu)})
  sensdR <- reactive({log(input$sensd)})
  resmuR <- reactive({log(input$resmu)})
  ressdR <- reactive({log(input$ressd)})
  
  sen_popR <- reactive({rlnormR(input$nn*(1-(input$prop_resist/100)),senmuR(),sensdR())})
  res_popR <- reactive({rlnormR(input$nn*(input$prop_resist/100),resmuR(),ressdR())})
  
  
  genData <- reactive({
    
    sen_pop <- sen_popR() #sensitive population
    res_pop <- res_popR() #resistant population
    c(sen_pop,res_pop)
    #total_pop <- c(sen_pop,res_pop)
  })
  
  
  output$histoplot1 <- renderPlot({
    #for now it only works for 2 distributions (sensitive and resistant distributions)
    
    plam<-c(1-(input$prop_resist/100),(input$prop_resist/100))
    pmu<-c(senmuR(),resmuR())
    psig<-c(sensdR(),ressdR())
    
    hist(genData(),freq=FALSE,main = paste("Distribution of parasite clearance half lives","\n", "from the parameters"),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20) #taken out for shiny #,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12)
    
    x <- seq(0.1, max(genData()), length=1000)
    hx <- list()
    
    #casting multiple lines for different distributions
    lcolor <- c('blue','red','magenta','brown','green')
    for(k in 1:length(pmu)){
      hx[[k]]<-plam[k]*dlnormR(x,meanlog=(pmu[k]),sdlog=psig[k])
      lines(x,hx[[k]],col=lcolor[k], lwd=5)
    }
    if(input$showcutoff){abline(v=input$cutoff, col='red', lwd=3, lty=3)}
    
  })
  
  genData.DF <- reactive({
    cbind(genData(),c(rep(0,length(sen_popR())),rep(1,length(res_popR()))))
  })
  output$ROC <- renderPlot({
    popDF <- genData.DF()
    
    TPR <- sum(res_popR()>=input$cutoff)/length(res_popR())
    FPR <- sum(sen_popR()>=input$cutoff)/length(sen_popR())
    
    true_res <- sum(res_popR()>=input$cutoff)
    fal_res <- sum(sen_popR()>=input$cutoff)
    fal_sen <- sum(res_popR()<input$cutoff)
    true_sen <- sum(sen_popR()<input$cutoff)
    overlay <- paste(true_res," truly resistant, ", fal_res, " falsely resistant \n",fal_sen," falsely sensitive, ",true_sen," truly sensitive")
    #overlay2 <- paste(true_res," truly resistant, ", fal_res, " falsely resistant ",fal_sen," falsely sensitive, ",true_sen," truly sensitive")
    #for output in the text #comment out the line above and below to get them back
    output$capturedO <- renderText(as.character(round(sum(genData()>=input$cutoff)/length(genData())*100,1)))
    #output$overlayO <- renderText(overlay2)
    
    roc(popDF[,2], popDF[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, show.thres=TRUE, main="Receiver Operating Characteristic (ROC) Curve")
    points((1-FPR),TPR, col="red", pch=19)
    text(.5,.5,overlay, col="red")
  })
  
  #######################################################################
  ###For the users data, run the mixture model and draw the histogram####
  #######################################################################
  mixdatR <- reactive({
    read.csv(input$file$datapath)
  })
  
  MixModel <- reactive({
    mixdat <- mixdatR()
    
    N<-ncol(mixdat)
    M <- 5
    
    pval<-0.1
    nboot<-100 # number of iterations for bootstrap
    nsim<-1000 # number of iterations for creating probaility of resistance vs HL graphs
    P<-2 # use P or more samples to get geometric means and discard all other samples for permutation analysis
    T<-100 # number of permutations for permutation analysis
    smax=5000
    
    # create output matrices
    output.mu <- matrix(NA,nrow=M,ncol=N)
    output.sigma <- matrix(NA,nrow=M,ncol=N)
    output.lambda <- matrix(NA,nrow=M,ncol=N)
    output.loglik <- matrix(NA,nrow=M,ncol=N)
    output.mu.se <- matrix(NA,nrow=M,ncol=N)
    output.sigma.se <- matrix(NA,nrow=M,ncol=N)
    output.lambda.se <- matrix(NA,nrow=M,ncol=N)
    AIC<-matrix(0,nrow=M,ncol=N)
    AICdelta<-matrix(0,nrow=M,ncol=N)
    
    nb<-na.omit(mixdat[,N])
    
    # fit single component model
    for (i in 1:N){
      # 1 COMPONENT LOG NORMAL
      nmixdat<-na.omit(mixdat[,i])
      lmixdat<- log(nmixdat)
      xll<-fitdistr(lmixdat,"normal")
      output.loglik[1,i]<- xll$loglik
      output.mu[1,i]<-xll$estimate[1]
      output.lambda[1,i]<-1
      output.sigma[1,i]<-xll$estimate[2]
      output.mu.se[1,i]<-xll$sd[1]
      output.sigma.se[1,i]<-xll$sd[2]
      output.lambda.se[1,i]<-0
      AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
      AICdelta[1,i]<-0
    }
    
    # fit multiple component models sequentially
    for (i in 1:N){
      nmixdat<-na.omit(mixdat[,i])
      lmixdat<- log(nmixdat)
      # >=2 COMPONENTS LOG NORMAL
      j<-1
      # stop if j-component model is more parsimonious than (j-1)-compnent model
      while((j<=M-1) && AICdelta[j,i]<=pval){
        j<-j+1
        res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
        resboot <- boot.se(res, B = nboot)
        resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
        output.loglik[j,i]<-res$loglik
        AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
        AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
        if(AICdelta[j,i]<=pval){
          output.mu[1:j,i]<-res$mu
          output.sigma[1:j,i]<-res$sigma
          output.lambda[1:j,i]<-res$lambda
          output.mu.se[1:j,i]<-resboot$mu.se
          output.sigma.se[1:j,i]<-resboot$sigma.se
          output.lambda.se[1:j,i]<-resboot$lambda.se		
        }
      }
    }
    list(muR = na.omit(output.mu), sigmaR = na.omit(output.sigma), lambdaR = na.omit(output.lambda))
  })
  #pmu = c(1.392575, 1.947629),
  MixModelResult <- reactiveValues(Holder =  list(muR = c(1.392575, 1.947629),
                                                  sigmaR = c(0.4046619, 0.2105479),
                                                  lambdaR = c(0.231683825089153 ,0.768316174910847)))
                                     # list(muR = c(0.6931472, 1.791759),
                                     #             sigmaR = c(0.4046619, 0.2105479),
                                     #             lambdaR = c(0.768316174910847, 0.231683825089153)))
  observeEvent(input$file,{
    MixModelResult$Holder <- MixModel()
  })
  
  ############################################################
  ###Following has been put out of the previous renderPlot####
  ############################################################

  aboutDistribution <- function(){
    mm <- MixModelResult$Holder
    j <- length(mm$muR)
    k <- list()
    for (a in 1:j) {
      k[[a]] <- paste("Distribution",a,
                      "has a geometric mean of ", round(exp(mm$muR[a]),2),
                      "hours with SD ", round(exp(mm$sigmaR[a]),2),
                      "hours. Its contribution to composite distribution is ", round(mm$lambdaR[a]*100,2), "%."
      )
    }
    tmp <- capture.output(cat(unlist(k)))
    paste("From the current dataset, the model predicts ", j, 
          " component geometric mean half lives (hours).", tmp)
  }
  
  resultTable <- function(){
    mm <- MixModelResult$Holder
    j <- length(mm$muR)
    Distributions <- Mean_hours <- SD_hours <- Contribution_percent <- NA
    for (a in 1:j) {
      Distributions[a] <- paste("Distribution", a)
      Mean_hours[a] <- round(exp(mm$muR[a]),2)
      SD_hours[a] <- round(exp(mm$sigmaR[a]),2)
      Contribution_percent[a] <- round(mm$lambdaR[a]*100,2)
    }
    as.data.frame(cbind(Distributions, Mean_hours, SD_hours, Contribution_percent))
  }
  
  
  
  mixdatRHolder <- reactiveValues(Holder=read.csv("simulated_data_for_input.csv", header=F))
  observeEvent(input$file,{
    mixdatRHolder$Holder <- mixdatR()
  })
  
  histoplot2R <- reactive({
    nmixdat<-na.omit(mixdatRHolder$Holder[,1])
    
    plam<-MixModelResult$Holder$lambdaR
    pmu<-MixModelResult$Holder$muR
    psig<-MixModelResult$Holder$sigmaR
    hist(nmixdat,freq=FALSE,main = paste("Distribution of parasite clearance half lives","\n", "from the data input"),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20) #taken out for shiny #,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12)
    
    x <- seq(0.1, max(nmixdat), length=1000)
    hx <- list()
    
    #casting multiple lines for different distributions
    lcolor <- c('blue','red','magenta','brown','green')
    for(k in 1:length(pmu)){
      hx[[k]]<-plam[k]*dlnormR(x,meanlog=(pmu[k]),sdlog=psig[k])
      lines(x,hx[[k]],col=lcolor[k], lwd=5)
    }
    if(input$showcutoff2){abline(v=input$cutoff2, col='red', lwd=3, lty=3)}
  })
  
  
  output$histoplot2 <- renderPlot({
    histoplot2R()
  })
  
  #for downloading the histoplot2 (plot from the user data)
  histoplot2fun <- function(){
    #histoplot2R()
    nmixdat<-na.omit(mixdatRHolder$Holder[,1])
    
    plam<-MixModelResult$Holder$lambdaR
    pmu<-MixModelResult$Holder$muR
    psig<-MixModelResult$Holder$sigmaR
    hist(nmixdat,freq=FALSE,main = paste("Distribution of parasite clearance half lives","\n", "from the data input"),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20) #taken out for shiny #,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12)
    
    x <- seq(0.1, max(nmixdat), length=1000)
    hx <- list()
    
    #casting multiple lines for different distributions
    lcolor <- c('blue','red','magenta','brown','green')
    for(k in 1:length(pmu)){
      hx[[k]]<-plam[k]*dlnormR(x,meanlog=(pmu[k]),sdlog=psig[k])
      lines(x,hx[[k]],col=lcolor[k], lwd=5)
    }
    if(input$showcutoff2){abline(v=input$cutoff2, col='red', lwd=3, lty=3)}
  }
  
  ###downloadHandlers####
  output$downloadhistoplot2 <- downloadHandler(
    filename = function(){paste('histogram_',Sys.Date(),'.png',sep='')},
    content = function(file) {
      png(file) #if(...=="png"){png(file)} else if(...=="pdf"){pdf(file)}
      histoplot2fun()
      dev.off()
    })
  output$defaultData <- downloadHandler(
    filename = function(){paste('hl_data_', Sys.Date(),'.csv',sep='')},
    content = function(file){
      write.table(mixdatRHolder$Holder, file, col.names = FALSE, row.names = FALSE)
    }
  )
  output$resultData <- downloadHandler(
    filename = function(){paste('result_', Sys.Date(),'.csv',sep='')},
    content = function(file){
      write.csv(resultTable(), file) #, col.names = FALSE, row.names = FALSE)
    }
  )
  
  proportionsmm <- reactive({MixModelResult$Holder$lambdaR})
  nn.mixdatmm <- reactive({nrow(mixdatRHolder$Holder)})
  #################################
  ####plots from the cutoff app####
  ####the functions for plotting###
  #################################
  
  senmuRmm <- reactive({MixModelResult$Holder$muR[1]})
  sensdRmm <- reactive({MixModelResult$Holder$sigmaR[1]})    
  resmuRmm <- reactive({MixModelResult$Holder$muR[2]})
  ressdRmm <- reactive({MixModelResult$Holder$sigmaR[2]})
  
  ###test####
  output$genDataOut <- renderPrint({as.character(names(MixModelResult$Holder))})
  #output$genDataOut <- renderPrint({as.character(c(nn.mixdatmm(),MixModelResult$Holder$muR,proportionsmm(),senmuRmm(),sensdRmm()))})
  #output$genDataOut <- renderPrint({class(na.omit(mixdatRHolder$Holder))})
  # output$densityLine <- renderPlot({
  #   
  # })
  
  sen_popRmm <- reactive({rlnormR(nn.mixdatmm()*(1-proportionsmm()[2]),senmuRmm(),sensdRmm())})
  res_popRmm <- reactive({rlnormR(nn.mixdatmm()*proportionsmm()[2],resmuRmm(),ressdRmm())})
  
  
  genDatamm <- reactive({
    
    sen_pop <- sen_popRmm() #sensitive population
    res_pop <- res_popRmm() #resistant population
    c(sen_pop,res_pop)
    #total_pop <- c(sen_pop,res_pop)
  })
  
  
  
  genData.DFmm <- reactive({
    cbind(genDatamm(),c(rep(0,length(sen_popRmm())),rep(1,length(res_popRmm()))))
  })
  output$ROC2 <- renderPlot({
    if(length(MixModelResult$Holder$muR)==2)
    { #plot ROC only if the number of distributions is <=2!
      rvExplain$text <- aboutDistribution() #geometric_means_and_proportions() #"rvExplain$text for 2 distributions" #renderText({      })
      popDF <- genData.DFmm()
      
      TPR <- sum(res_popRmm()>=input$cutoff2)/length(res_popRmm())
      FPR <- sum(sen_popRmm()>=input$cutoff2)/length(sen_popRmm())
      
      true_res <- sum(res_popRmm()>=input$cutoff2)
      fal_res <- sum(sen_popRmm()>=input$cutoff2)
      fal_sen <- sum(res_popRmm()<input$cutoff2)
      true_sen <- sum(sen_popRmm()<input$cutoff2)
      overlay <- paste(true_res," truly resistant, ", fal_res, " falsely resistant \n",fal_sen," falsely sensitive, ",true_sen," truly sensitive")
      
      
      roc(popDF[,2], popDF[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, show.thres=TRUE, main="Receiver Operating Characteristic (ROC) Curve")
      points((1-FPR),TPR, col="red", pch=19)
      text(.5,.5,overlay, col="red")
    }
    else if(length(MixModelResult$Holder$muR)==1){
      rvExplain$text="The model predicts a single distribution!"
      frame()
      title(main="ROC curve can't be plotted \n since the model predicts a single distribution!")
    }
    else if(length(MixModelResult$Holder$muR)>2){
      rvExplain$text <- aboutDistribution()
      frame()
      title(main="ROC curve can't be plotted since \n the number of distributions is more than 2!")
    }
    else {
      frame()
      title(main="ROC curve can't be plotted \n either because there's no data or \n the data is in the wrong format!")
    }
  })
}

shinyApp(server = server, ui = ui)