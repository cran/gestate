## display options: align radio buttons in 2 columns
multicollab <- list(tags$head(tags$style(HTML(
  '.shiny-options-group{ 
  height: auto; width: auto; margin-top: 0px;
  -webkit-column-count: 2; -moz-column-count: 2; column-count: 2; 
  -webkit-column-fill: balance; -moz-column-fill: balance; column-fill: balance;
  } 
  .control-label{
  padding-bottom: 10px;
  }
  div.radio {
#  margin-top: 0px; margin-bottom: 0px; padding-bottom: 5px;
   margin-top: 0px; margin-bottom: 0px; padding-bottom: 0px;
  }'
))))

get_vignette_link <- function(...) {
  x <- vignette(...)
  if (nzchar(out <- x$PDF)) {
    ext <- tools::file_ext(out)
    port <- if (tolower(ext) == "html") 
      tools::startDynamicHelp(NA)
    else 0L
    if (port > 0L) {
      out <- sprintf("http://127.0.0.1:%d/library/%s/doc/%s", 
              port, basename(x$Dir), out)
      return(out)
    }
  }
  stop("no html help found")
}

## user interface
navbarPage(theme = shinythemes::shinytheme('flatly'), # layout
           paste("V",noquote(installed.packages()["gestate" ,"Version"]),sep=""), # title
           
           ## first tab: Analytical calculations & simulation
           tabPanel('Analytical Power Calculation & Simulation of Time-To-Event Clinical Trials',
                    sidebarLayout(
                      
                      ## input panel on the left
                      sidebarPanel( 
                        
                        ## input widgets
                        
                        fluidRow(column(12, tags$h3('Distributions'))),
                        br(),

                        fluidRow(column(12, tags$h4('Event Distributions'))),
                        fluidRow(column(12, tags$h4('Active Arm'))),
                        multicollab, 
                        fluidRow(column(12, radioButtons('event_active_distribution_ae', 'Distribution type', choices = list('Weibull'='weibull', 'Lognormal'='lognormal', 'Exponential'='exponential', 'Piecewise exponential'='pieceexponential', 'Log-Logistic'='loglogistic', 'Gompertz'='gompertz', 'Generalised Gamma'='ggamma'), selected='weibull'))),
                        fluidRow(column(12, uiOutput('event_active_type1_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, uiOutput('event_active_type2_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, tags$h4('Control Arm'))),
                        fluidRow(column(12, radioButtons('event_control_distribution_ae', 'Distribution type', choices = list('Weibull'='weibull', 'Lognormal'='lognormal', 'Exponential'='exponential', 'Piecewise exponential'='pieceexponential', 'Log-Logistic'='loglogistic', 'Gompertz'='gompertz', 'Generalised Gamma'='ggamma'), selected='weibull'))),
                        fluidRow(column(12, uiOutput('event_control_type1_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, uiOutput('event_control_type2_ae'))), # dynamic input field. See server.R
                        br(),
                        fluidRow(column(12, tags$h4('Censoring Distributions'))), 
                        fluidRow(column(12, tags$h4('Active Arm'))),
                        fluidRow(column(12, radioButtons('cens_active_distribution_ae', 'Distribution type', choices = list('Weibull'='weibull', 'Lognormal'='lognormal', 'Exponential'='exponential', 'Piecewise exponential'='pieceexponential', 'Log-Logistic'='loglogistic', 'Gompertz'='gompertz', 'Generalised Gamma'='ggamma', 'No censoring'='blank'), selected='blank'))),
                        fluidRow(column(12, uiOutput('cens_active_type1_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, uiOutput('cens_active_type2_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, tags$h4('Control Arm'))),
                        fluidRow(column(12, radioButtons('cens_control_distribution_ae', 'Distribution type', choices = list('Weibull'='weibull', 'Lognormal'='lognormal', 'Exponential'='exponential', 'Piecewise exponential'='pieceexponential', 'Log-Logistic'='loglogistic', 'Gompertz'='gompertz', 'Generalised Gamma'='ggamma', 'No censoring'='blank'), selected='blank'))),
                        fluidRow(column(12, uiOutput('cens_control_type1_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, uiOutput('cens_control_type2_ae'))), # dynamic input field. See server.R
                        br(),
                        fluidRow(column(12, tags$h4('Recruitment Distribution'))),
                        fluidRow(column(12, radioButtons('recruit_distribution_ae', 'Distribution type', choices = list('Instant recruitment' = 'instantRec', 'Linear recruitment'='linearRec', 'Piecewise recruitment'='pieceRec'), selected='linearRec'))),
                        fluidRow(column(12, uiOutput('recruit_type1_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, uiOutput('recruit_type2_ae'))), # dynamic input field. See server.R
                        
                        br(), br(),
                        
                        downloadButton('downloadPlots_dist', label='Download Distribution Plots'),
                        
                        br(), br(),
                        
                        fluidRow(column(12, tags$h3('Evaluation Parameters'))),
                        br(),
                        fluidRow(column(6 , numericInput('alpha1_ae', HTML('Type I error (1-sided, <br/>for power & RMST calculations)'), min=0.001, max=0.5, value=0.025, step=0.001)),
                                 column(6 , numericInput('power_ae', HTML('Power <br/>(for sample size calculation)'), min=0.100, max=0.999, value=0.900, step=0.001))),
                        fluidRow(column(12, tags$h4('RMST Calculations'))),
                        fluidRow(column(12, radioButtons('RMST_ae', 'Calculate RMST properties', choices = list('Yes'=TRUE, 'No'=FALSE), selected=FALSE))),
                        fluidRow(column(12, uiOutput('restrictiontime_ae'))), # dynamic input field. See server.R
                        fluidRow(column(12, tags$h4('Landmark Calculations'))),
                        fluidRow(column(12, radioButtons('landmark_ae', 'Calculate Landmark properties', choices = list('Yes'=TRUE, 'No'=FALSE), selected=FALSE))),
                        fluidRow(column(12, uiOutput('landmarktime_ae'))), # dynamic input field. See server.R
                        br(),
                        fluidRow(column(12, tags$h3('Analytic Approach'))),
                        br(),
                        fluidRow(column(6, numericInput('assess_ae', 'Maximum assessment time (in months)', min=0, max=300, value=50, step=1))),
                        br(),
                        actionButton('input_ae', 'Analyse'),
                        br(), br(),
                        downloadButton('downloadData_ae', label='Download Analytic Tables'),
                        downloadButton('downloadPlots_ae', label='Download Analytic Plots'),
                        
                        br(), br(),
                        
                        fluidRow(column(12, tags$h3('Simulation Approach'))),
                        br(),
                        fluidRow(column(6, numericInput('assess_sim_ae', 'Assessment time (in months)', min=0, max=300, value=24, step=1))),
                        br(),
                        fluidRow(column(6 , numericInput('iterations_sim', 'Number of iterations', min=10, max=100000, value=5000, step=1)),
                                 column(6 , numericInput('seed_sim', 'Seed number', min=1, max=100000, value=123456, step=1))),
                        fluidRow(column(12, radioButtons('create_events_sim', 'Fix event number', choices = list('Yes'=TRUE, 'No'=FALSE), selected=FALSE))),
                        fluidRow(column(12, uiOutput('event_creation_sim'))), # dynamic input field. See server.R
                        fluidRow(column(12, radioButtons('parallel_analysis_sim', 'Analysis by parallel processing', choices = list('Yes'=TRUE, 'No'=FALSE), selected=FALSE))),
                        fluidRow(column(12, uiOutput('core_number_sim'))), # dynamic input field. See server.R
                        fluidRow(column(12, uiOutput('LRCox_sim'))), # dynamic input field. See server.R
                        fluidRow(column(6, numericInput('displayData_sim', 'Display until iteration:', min=1, max=5000, value=1, step=1))),
                        br(), 
                        actionButton('input_sim', 'Perform Simulation'),
                        br(), br(),
                        downloadButton('downloadData_sim',label='Download Simulation Data'),
                        br(),
                        downloadButton('downloadAna_sim', label='Download Simulation Analysis'),
                        br(),
                        downloadButton('downloadSum_sim', label='Download Simulation Summary')
                        
                      ), 
                      
                      ## output panel on the right
                      mainPanel(
                        
                        ## output fields
                        
                        fluidRow(column(12, h4('Planned Survival Curves'), br(), plotOutput('KMCurvePlot', width="50%"), br())),
                        fluidRow(column(12, h4('Planned Censoring Curves'), br(), plotOutput('CensCurvePlot', width="50%"), br())),
                        fluidRow(column(12, h4('Planned Recruitment Plot'), br(), plotOutput('RecCurvePlot', width="50%"), br(), br())),
                        fluidRow(column(12, h3('Analysis'), br())),
                        fluidRow(column(12, p('Analytic results will appear in this section when the Analyse button is pressed.'), br())),      
                        fluidRow(column(12, h4('Values'), br())),
                        fluidRow(column(12, textOutput('SS_caveat_ae'), br())),
                        fluidRow(column(12, textOutput('text_ae'), br())),
                        fluidRow(column(10, tableOutput('Analysis_table'), br(), br())),
                        fluidRow(column(12, h4('Plots'), br(), plotOutput('ObservedEventsPlot', width="50%"),br())),
                        fluidRow(column(12, plotOutput('ExpectedLogHRPlot',width="50%"),br())),
                        fluidRow(column(12, plotOutput('PowerPlot',width="50%"),br())),
                        
                        fluidRow(column(12, h3('Simulation'), br())),
                        fluidRow(column(12, p('Simulation results will appear in this section when the Perform Simulation button is pressed. Summary contains the overall properties of the simulation following analysis. Analysis contains the analysis results for each simulation  requested to be displayed (all results may be found in the download). Values contains the patient-level data for each simulation requested to be displayed (all results may be found in the download). '), br())),
                        fluidRow(column(10, h4('Summary'), br(), tableOutput('summarised_analysis'), br())),
                        fluidRow(column(10, h4('Analysis'), br(), tableOutput('Simulation_analysis'), br())),
                        fluidRow(column(10, h4('Values'), br(), tableOutput('Simulation_table')))
                      )
                      
                    )
           ),
           
           
           ## second tab: about the Gestate Shiny app
           navbarMenu('Help',
                      
                      ## manual         
                      tabPanel('Manual', br(), h3('User Manual'), br(), a('Trial Planning User Manual', target='_self', href=get_vignette_link("trial_planning", package="gestate"))),
                      ## general information
                      tabPanel('About Gestate Shiny',
                               
                               br(), h3('R Shiny App for TTE calculation under NPH.'), br(),
                               h5(paste("Released as part of the gestate V",noquote(installed.packages()["gestate" ,"Version"]),' package.',sep="")), br(),
                               h5('Authors: James Bell (james.bell.ext@boehringer-ingelheim.com), Jasmin Ruehl'), br(), br(), br(), br(),
                               
                               h4('Details:'), br(),
                               p('The gestate package is a platform for design of generalised Time-To-Event trials and event prediction in ongoing trials.
                                  It supports a wide range of event, censoring and recruitment distributions through the use of Curve and RCurve objects. 
                                  The package is aimed at the design of trials with non-proportional hazards and/or less standard event distributions and/or complex censoring assumptions.
                                  This Shiny app is distributed with the gestate package and provides a GUI for most of the package trial-design functionality.
                                  The GUI is intended to give new users an introduction to the package or to allow more experienced users a means of casual exploration of possibilities by providing real-time updated plots of distribution functions.
                                  For thorough exploration of design space and/or integration into workflows, it is recommended to use the package functionality through R directly. '), br(), br(),
                         
                               h4('Shiny App Functionality:'), br(), br(),
                               p('Analytic estimation of HR and event numbers over time'),
                               p('Analytic power calculation (via Schoenfeld, Frontier methods)'),
                               p('Sample Size estimation under Schoenfeld (confirmation by subsequent power calculation recommended)'),  
                               p('Analytic RMST and Landmark expectations of point estimate and variance, with power calculation'),
                               p('2-arm Trial Simulation under assumption of Censoring Completely At Random'),
                               p('Automatic simulation analysis and result summarisation for log-rank, Cox, RMST and Landmark methods'),br(),
                               h4('Curve distribution parameterisations:'), br(),
                               p('Exponential:',br(), 'f(x) = lambda e^(- lambda x)', br(), 'F(x) = 1 - e^(- lambda x)'), br(),
                               p('Weibull:',br(), ' f(x) = (beta/alpha) (x/alpha)^(beta-1) exp(-(x/alpha)^beta)', br(), 'F(x)= 1 - exp(-(x/alpha)^beta)'), br(),    
                               p('Log-Normal:',br(), ' f(x) = 1/(v(2*pi) sigma  x) e^-((log x - mu)^2 / (2 sigma^2))', br(), 'F(x)=0.5(1 + erf((log(x)-mu)/(sigma sqrt(2))))'), br(),          
                               p('Log-Logistic:',br(), ' f(x) = eta (theta^beta) x^(eta-1) (theta^eta + x^eta)^-2 ', br(), ' F(x) = (x^eta) /(theta^eta+x^eta)'), br(),
                               p('Gompertz:',br(), ' f(x) = theta eta e^(eta + theta x - eta e^(theta x))', br(), ' F(x) = 1 - exp(eta - eta e^(theta x))'), br(),      
                               p('Generalised Gamma:',br(), ' f(x) = (rho x^((rho eta)-1) e^(-(x/theta)^rho) theta^(-rho eta) )/Gamma(eta) ', br(), ' F(x) = LPGamma(eta,(x/theta)^rho)/Gamma(eta)'), br(),    
                               p('Blank:', br(), 'f(x) = 0', br(), 'F(x) = 0'), br() 
                               
                      )
                      
           )
)
