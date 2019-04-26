## handle events and create output
function(input, output, session) {

  ### INPUT Analytical calculations & simulation #######################################
  ## distributions #####################################################################

  #Prespecify the defaults to avoid errors
  alpha_event_active_ae <- 1
  beta_event_active_ae <- 1
  mu_event_active_ae <- 0
  sigma_event_active_ae <- 1
  lambda_event_active_ae <- 1
  number_pieces_event_active_ae <- 1
  theta_ll_event_active_ae <- 1
  eta_ll_event_active_ae <- 1
  theta_gomp_event_active_ae <- 1
  eta_gomp_event_active_ae <- 1
  theta_gg_event_active_ae <- 1
  eta_gg_event_active_ae <- 2
  rho_gg_event_active_ae <- 1

  ## change event distribution parameter mask depending on distribution type (active arm)
  output$event_active_type1_ae <- renderUI({
    switch(input$event_active_distribution_ae,
           'weibull'          = fluidRow(column(6 , numericInput('alpha_event_active_ae', 'Scale parameter', min=0.001, max=1000, value=10, step=0.001)),
                                         column(6 , numericInput('beta_event_active_ae', 'Shape parameter',  min=0.001, max=1000, value=1, step=0.001))),
           'lognormal'        = fluidRow(column(6 , numericInput('mu_event_active_ae', 'Meanlog parameter', min=-1000, max=1000, value=2, step=0.01)),
                                         column(6 , numericInput('sigma_event_active_ae', 'Sdlog parameter',  min=0.001, max=1000, value=1, step=0.001))),
           'exponential'      = fluidRow(column(6 , numericInput('lambda_event_active_ae', 'Parameter', min=0.0001, max=1000, value=0.1, step=0.0001))),
           'pieceexponential' = fluidRow(column(6 , helpText('Start times must be in ascending order with no duplicates.')),
                                         column(6 , numericInput('number_pieces_event_active_ae', 'Number of intervals', min=1, max=20, value=1, step=1))),
           'loglogistic'      = fluidRow(column(6 , numericInput('theta_ll_event_active_ae', 'Scale parameter', min=0.001, max=1000, value=10, step=0.001)),
                                         column(6 , numericInput('eta_ll_event_active_ae', 'Shape parameter', min=0.001, max=1000, value=1, step=0.001))),
           'gompertz'         = fluidRow(column(6 , numericInput('theta_gomp_event_active_ae', 'Scale parameter', min=0.0001, max=1000, value=0.05, step=0.0001)),
                                         column(6 , numericInput('eta_gomp_event_active_ae', 'Shape parameter', min=0.001, max=1000, value=1, step=0.001))),
           'ggamma'           = fluidRow(column(6 , numericInput('theta_gg_event_active_ae', 'Scale parameter', min=0.001, max=1000, value=20, step=0.001)),
                                         column(6 , helpText('Shape parameter must be greater than 1.')))
    )
  })

 # if(input$alpha_event_active_ae < 0.001){updateNumericInput(session,input$alpha_event_active_ae,value=14)}

  ## change event distribution parameter mask depending on distribution type (active arm)
  output$event_active_type2_ae <- renderUI({
    if(input$event_active_distribution_ae == 'weibull'){
    }else if(input$event_active_distribution_ae == 'lognormal'){
    }else if(input$event_active_distribution_ae == 'exponential'){
    }else if(input$event_active_distribution_ae == 'pieceexponential'){
      expperiods <- ifelse(is.null(input$number_pieces_event_active_ae), 0, as.numeric(input$number_pieces_event_active_ae))
      if(expperiods >=1){
        lapply(1:expperiods, function(i){
          if(i > 1){
            fluidRow(column(6 , numericInput(inputId=paste0('start', i, '_event_active_ae'), label=paste0('Interval ', i, ': Start time'), min=1, max=300, value=i-1, step=1)),
                     column(6 , numericInput(inputId=paste0('lambda', i, '_event_active_ae'), label=paste0('Interval ', i, ': Parameter'), min=0.001, max=1000, value=0.1, step=0.001)))
          }else{
            fluidRow(column(6 , strong('Interval 1: Start time'), br(), h5('0')),
                     column(6 , numericInput(inputId='lambda1_event_active_ae', label='Interval 1: Parameter', min=0.001, max=1000, value=0.1, step=0.001)))
          }
        })
      }
    }else if(input$event_active_distribution_ae == 'loglogistic'){
    }else if(input$event_active_distribution_ae == 'gompertz'){
    }else if(input$event_active_distribution_ae == 'ggamma'){
      fluidRow(column(6 , numericInput('rho_gg_event_active_ae', 'Family parameter', min=0.001, max=1000, value=2, step=0.001)),
               column(6 , numericInput('eta_gg_event_active_ae', 'Shape parameter', min=1.001, max=1000, value=1.001, step=0.001)))
    }
  })


  #Prespecify the defaults to avoid errors
  alpha_event_control_ae <- 1
  beta_event_control_ae <- 1
  mu_event_control_ae <- 0
  sigma_event_control_ae <- 1
  lambda_event_control_ae <- 1
  number_pieces_event_control_ae <- 1
  theta_ll_event_control_ae <- 1
  eta_ll_event_control_ae <- 1
  theta_gomp_event_control_ae <- 1
  eta_gomp_event_control_ae <- 1
  theta_gg_event_control_ae <- 1
  eta_gg_event_control_ae <- 2
  rho_gg_event_control_ae <- 1

  ## change event distribution parameter mask depending on distribution type (control arm)

  output$event_control_type1_ae <- renderUI({
    switch(input$event_control_distribution_ae,
           'weibull'          = fluidRow(column(6 , numericInput('alpha_event_control_ae', 'Scale parameter', min=0.001, max=1000, value=10, step=0.001)),
                                         column(6 , numericInput('beta_event_control_ae', 'Shape parameter',  min=0.001, max=1000, value=1, step=0.001))),
           'lognormal'        = fluidRow(column(6 , numericInput('mu_event_control_ae', 'Meanlog parameter', min=-1000, max=1000, value=2, step=0.01)),
                                         column(6 , numericInput('sigma_event_control_ae', 'Sdlog parameter',  min=0.001, max=1000, value=1, step=0.001))),
           'exponential'      = fluidRow(column(6 , numericInput('lambda_event_control_ae', 'Parameter', min=0.0001, max=1000, value=0.1, step=0.0001))),
           'pieceexponential' = fluidRow(column(6 , helpText('Start times must be in ascending order with no duplicates.')),
                                         column(6 , numericInput('number_pieces_event_control_ae', 'Number of intervals', min=1, max=20, value=1, step=1))),
           'loglogistic'      = fluidRow(column(6 , numericInput('theta_ll_event_control_ae', 'Scale parameter', min=0.001, max=1000, value=10, step=0.001)),
                                         column(6 , numericInput('eta_ll_event_control_ae', 'Shape parameter', min=0.001, max=1000, value=1, step=0.001))),
           'gompertz'         = fluidRow(column(6 , numericInput('theta_gomp_event_control_ae', 'Scale parameter', min=0.0001, max=1000, value=0.05, step=0.0001)),
                                         column(6 , numericInput('eta_gomp_event_control_ae', 'Shape parameter', min=0.001, max=1000, value=1, step=0.001))),
           'ggamma'           = fluidRow(column(6 , numericInput('theta_gg_event_control_ae', 'Scale parameter', min=0.001, max=1000, value=20, step=0.001)),
                                         column(6 , helpText('Shape parameter must be greater than 1.')))
    )
  })

  ## change event distribution parameter mask depending on distribution type (control arm)
  output$event_control_type2_ae <- renderUI({
    if(input$event_control_distribution_ae == 'weibull'){
    }else if(input$event_control_distribution_ae == 'lognormal'){
    }else if(input$event_control_distribution_ae == 'exponential'){
    }else if(input$event_control_distribution_ae == 'pieceexponential'){
      expperiods <- ifelse(is.null(input$number_pieces_event_control_ae), 0, as.numeric(input$number_pieces_event_control_ae))
      lapply(1:expperiods, function(i){
        if(i > 1){
          fluidRow(column(6 , numericInput(inputId=paste0('start', i, '_event_control_ae'), label=paste0('Interval ', i, ': Start time'), min=1, max=300, value=i-1, step=1)),
                   column(6 , numericInput(inputId=paste0('lambda', i, '_event_control_ae'), label=paste0('Interval ', i, ': Parameter'), min=0.001, max=1000, value=0.1, step=0.001)))
        }else{
          fluidRow(column(6 , strong('Interval 1: Start time'), br(), h5('0')),
                   column(6 , numericInput(inputId='lambda1_event_control_ae', label='Interval 1: Parameter', min=0.001, max=1000, value=0.1, step=0.001)))
        }
      })
    }else if(input$event_control_distribution_ae == 'loglogistic'){
    }else if(input$event_control_distribution_ae == 'gompertz'){
    }else if(input$event_control_distribution_ae == 'ggamma'){
      fluidRow(column(6 , numericInput('rho_gg_event_control_ae', 'Family parameter', min=0.001, max=1000, value=2, step=0.001)),
               column(6 , numericInput('eta_gg_event_control_ae', 'Shape parameter', min=1.001, max=1000, value=1.001, step=0.001)))
    }
  })

  #Prespecify the defaults to avoid errors
  alpha_cens_active_ae <- 1
  beta_cens_active_ae <- 1
  mu_cens_active_ae <- 0
  sigma_cens_active_ae <- 1
  lambda_cens_active_ae <- 1
  number_pieces_cens_active_ae <- 1
  theta_ll_cens_active_ae <- 1
  eta_ll_cens_active_ae <- 1
  theta_gomp_cens_active_ae <- 1
  eta_gomp_cens_active_ae <- 1
  theta_gg_cens_active_ae <- 1
  eta_gg_cens_active_ae <- 2
  rho_gg_cens_active_ae <- 1

  ## change censoring distribution parameter mask depending on distribution type (active arm)
  output$cens_active_type1_ae <- renderUI({
    switch(input$cens_active_distribution_ae,
           'weibull'          = fluidRow(column(6 , numericInput('alpha_cens_active_ae', 'Scale parameter', min=0.001, max=1000, value=500, step=0.001)),
                                         column(6 , numericInput('beta_cens_active_ae', 'Shape parameter',  min=0.001, max=1000, value=0.3, step=0.001))),
           'lognormal'        = fluidRow(column(6 , numericInput('mu_cens_active_ae', 'Meanlog parameter', min=-1000, max=1000, value=5, step=0.01)),
                                         column(6 , numericInput('sigma_cens_active_ae', 'Sdlog parameter',  min=0.001, max=1000, value=2, step=0.001))),
           'exponential'      = fluidRow(column(6 , numericInput('lambda_cens_active_ae', 'Parameter', min=0.0001, max=1000, value=0.01, step=0.0001))),
           'pieceexponential' = fluidRow(column(6 , helpText('Start times must be in ascending order with no duplicates.')),
                                         column(6 , numericInput('number_pieces_cens_active_ae', 'Number of intervals', min=1, max=20, value=1, step=1))),
           'loglogistic'      = fluidRow(column(6 , numericInput('theta_ll_cens_active_ae', 'Scale parameter', min=0.001, max=1000, value=50, step=0.001)),
                                         column(6 , numericInput('eta_ll_cens_active_ae', 'Shape parameter', min=0.001, max=1000, value=2, step=0.001))),
           'gompertz'         = fluidRow(column(6 , numericInput('theta_gomp_cens_active_ae', 'Scale parameter', min=0.0001, max=1000, value=0.005, step=0.0001)),
                                         column(6 , numericInput('eta_gomp_cens_active_ae', 'Shape parameter', min=0.001, max=1000, value=2, step=0.001))),
           'ggamma'           = fluidRow(column(6 , numericInput('theta_gg_cens_active_ae', 'Scale parameter', min=0.001, max=1000, value=70, step=0.001)),
                                         column(6 , helpText('Shape parameter must be greater than 1.'))),
           'blank'            = NULL
    )
  })

  ## change censoring distribution parameter mask depending on distribution type (active arm)
  output$cens_active_type2_ae <- renderUI({
    if(input$cens_active_distribution_ae == 'weibull'){
    }else if(input$cens_active_distribution_ae == 'lognormal'){
    }else if(input$cens_active_distribution_ae == 'exponential'){
    }else if(input$cens_active_distribution_ae == 'pieceexponential'){
      expperiods <- ifelse(is.null(input$number_pieces_cens_active_ae), 0, as.numeric(input$number_pieces_cens_active_ae))
      lapply(1:expperiods, function(i){
        if(i > 1){
          fluidRow(column(6 , numericInput(inputId=paste0('start', i, '_cens_active_ae'), label=paste0('Interval ', i, ': Start time'), min=1, max=300, value=i-1, step=1)),
                   column(6 , numericInput(inputId=paste0('lambda', i, '_cens_active_ae'), label=paste0('Interval ', i, ': Parameter'), min=0.001, max=1000, value=0.1, step=0.001)))
        }else{
          fluidRow(column(6 , strong('Interval 1: Start time'), br(), h5('0')),
                   column(6 , numericInput(inputId='lambda1_cens_active_ae', label='Interval 1: Parameter', min=0.001, max=1000, value=0.1, step=0.001)))
        }
      })
    }else if(input$cens_active_distribution_ae == 'loglogistic'){
    }else if(input$cens_active_distribution_ae == 'gompertz'){
    }else if(input$cens_active_distribution_ae == 'ggamma'){
      fluidRow(column(6 , numericInput('rho_gg_cens_active_ae', 'Family parameter', min=0.001, max=100, value=2, step=0.001)),
               column(6 , numericInput('eta_gg_cens_active_ae', 'Shape parameter', min=1.001, max=100, value=1.001, step=0.001)))
    }else if(input$cens_active_distribution_ae == 'blank'){
    }
  })

    #Prespecify the defaults to avoid errors
  alpha_cens_control_ae <- 1
  beta_cens_control_ae <- 1
  mu_cens_control_ae <- 0
  sigma_cens_control_ae <- 1
  lambda_cens_control_ae <- 1
  number_pieces_cens_control_ae <- 1
  theta_ll_cens_control_ae <- 1
  eta_ll_cens_control_ae <- 1
  theta_gomp_cens_control_ae <- 1
  eta_gomp_cens_control_ae <- 1
  theta_gg_cens_control_ae <- 1
  eta_gg_cens_control_ae <- 2
  rho_gg_cens_control_ae <- 1

  ## change censoring distribution parameter mask depending on distribution type (control arm)
  output$cens_control_type1_ae <- renderUI({
    switch(input$cens_control_distribution_ae,
           'weibull'          = fluidRow(column(6 , numericInput('alpha_cens_control_ae', 'Scale parameter', min=0.001, max=1000, value=500, step=0.001)),
                                         column(6 , numericInput('beta_cens_control_ae', 'Shape parameter',  min=0.001, max=1000, value=0.3, step=0.001))),
           'lognormal'        = fluidRow(column(6 , numericInput('mu_cens_control_ae', 'Meanlog parameter', min=-1000, max=1000, value=5, step=0.01)),
                                         column(6 , numericInput('sigma_cens_control_ae', 'Sdlog parameter',  min=0.001, max=1000, value=2, step=0.001))),
           'exponential'      = fluidRow(column(6 , numericInput('lambda_cens_control_ae', 'Parameter', min=0.0001, max=1000, value=0.01, step=0.0001))),
           'pieceexponential' = fluidRow(column(6 , helpText('Start times must be in ascending order with no duplicates.')),
                                         column(6 , numericInput('number_pieces_cens_control_ae', 'Number of intervals', min=1, max=10, value=1, step=1))),
           'loglogistic'      = fluidRow(column(6 , numericInput('theta_ll_cens_control_ae', 'Scale parameter', min=0.001, max=1000, value=50, step=0.001)),
                                         column(6 , numericInput('eta_ll_cens_control_ae', 'Shape parameter', min=0.001, max=1000, value=2, step=0.001))),
           'gompertz'         = fluidRow(column(6 , numericInput('theta_gomp_cens_control_ae', 'Scale parameter', min=0.0001, max=1000, value=0.005, step=0.0001)),
                                         column(6 , numericInput('eta_gomp_cens_control_ae', 'Shape parameter', min=0.001, max=1000, value=2, step=0.001))),
           'ggamma'           = fluidRow(column(6 , numericInput('theta_gg_cens_control_ae', 'Scale parameter', min=0.001, max=1000, value=70, step=0.001)),
                                         column(6 , helpText('Shape parameter must be greater than 1.'))),
           'blank'            = NULL
    )
  })

  ## change censoring distribution parameter mask depending on distribution type (control arm)
  output$cens_control_type2_ae <- renderUI({
    if(input$cens_control_distribution_ae == 'weibull'){
    }else if(input$cens_control_distribution_ae == 'lognormal'){
    }else if(input$cens_control_distribution_ae == 'exponential'){
    }else if(input$cens_control_distribution_ae == 'pieceexponential'){
      expperiods <- ifelse(is.null(input$number_pieces_cens_control_ae), 0, as.numeric(input$number_pieces_cens_control_ae))
      lapply(1:expperiods, function(i){
        if(i > 1){
          fluidRow(column(6 , numericInput(inputId=paste0('start', i, '_cens_control_ae'), label=paste0('Interval ', i, ': Start time'), min=1, max=300, value=i-1, step=1)),
                   column(6 , numericInput(inputId=paste0('lambda', i, '_cens_control_ae'), label=paste0('Interval ', i, ': Parameter'), min=0.001, max=1000, value=0.1, step=0.001)))
        }else{
          fluidRow(column(6 , strong('Interval 1: Start time'), br(), h5('0')),
                   column(6 , numericInput(inputId='lambda1_cens_control_ae', label='Interval 1: Parameter', min=0.001, max=1000, value=0.1, step=0.001)))
        }
      })
    }else if(input$cens_control_distribution_ae == 'loglogistic'){
    }else if(input$cens_control_distribution_ae == 'gompertz'){
    }else if(input$cens_control_distribution_ae == 'ggamma'){
      fluidRow(column(6 , numericInput('rho_gg_cens_control_ae', 'Family parameter', min=0.001, max=1000, value=2, step=0.001)),
               column(6 , numericInput('eta_gg_cens_control_ae', 'Shape parameter', min=1.001, max=1000, value=1.001, step=0.001)))
    }else if(input$cens_control_distribution_ae == 'blank'){
    }
  })

  ## change recruitment distribution parameter mask depending on distribution type
  output$recruit_type1_ae <- renderUI({
    switch(input$recruit_distribution_ae,
           'instantRec' = fluidRow(column(6 , numericInput('Nactive_instant_ae', 'Number of patients (active arm)', min=0, max=25000, value=30, step=1)),
                                   column(6 , numericInput('Ncontrol_instant_ae', 'Number of patients (control arm)',  min=0, max=25000, value=30, step=1))),
           'linearRec'  = fluidRow(column(6 , numericInput('rlength_linear_ae', 'Length of recruitment (in months)', min=1, max=300, value=30, step=1))),
           'pieceRec'   = fluidRow(column(6 , numericInput('recruitperiods_piece_ae', 'Number of recruitment periods', min=1, max=30, value=1, step=1)),
                                   column(6 , numericInput('ratio_piece_ae', 'Recruitment ratio (active / control)', min=0.01, max=100, value=1, step=0.001)))
    )
  })

  ## change recruitment distribution parameter mask depending on distribution type
  output$recruit_type2_ae <- renderUI({
    if(input$recruit_distribution_ae == 'linearRec'){
      fluidRow(column(6 , numericInput('Nactive_linear_ae', 'Number of patients (active arm)', min=0, max=25000, value=30, step=1)),
               column(6 , numericInput('Ncontrol_linear_ae', 'Number of patients (control arm)',  min=0, max=25000, value=30, step=1)))
    }else if(input$recruit_distribution_ae == 'instantRec'){
    }else if(input$recruit_distribution_ae == 'pieceRec'){
      recruitperiods <- ifelse(is.null(input$recruitperiods_piece_ae), 0, as.numeric(input$recruitperiods_piece_ae))
      lapply(1:recruitperiods, function(i){
        fluidRow(column(6 , numericInput(inputId=paste0('periodlength', i, '_piece_ae'), label=paste0('Period ', i, ': length (in months)'), min=1, max=300, value=5, step=1)),
                 column(6 , numericInput(inputId=paste0('periodrate', i, '_piece_ae'), label=paste0('Period ', i, ': rate (in patient numbers)'), min=0, max=1000, value=10, step=1)))
      })
    }
  })

  ## evaluation ########################################################################

  ## change RMST calculation parameter mask if selected
  output$restrictiontime_ae <- renderUI({
    if(input$RMST_ae){
      fluidRow(column(12, numericInput('restrictionTime_ae', 'Restriction time (in months)', min=0, max=300, value=10, step=1)))
    }
  })

  ## change landmark calculation parameter mask if selected
  output$landmarktime_ae <- renderUI({
    if(input$landmark_ae){
      fluidRow(column(12, numericInput('landmarkTime_ae', 'Time for Landmark calculation (in months)', min=0, max=300, value=10, step=1)))
    }
  })

  ## simulation ########################################################################

  ## change event fixing parameter mask if selected
  output$event_creation_sim <- renderUI({
    if(input$create_events_sim){
      fluidRow(column(12, numericInput('force_events_sim', 'Number of events to force into each simulation', min=0, max=25000, value=100, step=1)))
    }
  })

  ## change parallel processing parameter mask if selected
  output$core_number_sim <- renderUI({
    if(input$parallel_analysis_sim){
      fluidRow(column(12, numericInput('cores_sim', 'Number of cores for parallel processing', min=1, max=16, value=4, step=1)))
    }
  })

  output$LRCox_sim <- renderUI({
    if(as.logical(input$RMST_ae) | as.logical(input$landmark_ae)){
      fluidRow(column(12, radioButtons('lrcox_sim', 'Perform Log-Rank/ Cox analysis', choices = list('Yes'=TRUE, 'No'=FALSE), selected=TRUE)))
    }else{
      fluidRow(column(12, radioButtons('lrcox_sim', 'Perform Log-Rank/ Cox analysis', choices = list('Yes'=TRUE), selected=TRUE)))
    }
  })

  ### OUTPUT Analytical calculations & simulation ######################################

  ## plot survival curves depending on assumptions
  output$KMCurvePlot <- renderPlot({
    ## set distributions
    ## piecewise exponential distribution
    if(input$event_active_distribution_ae == 'pieceexponential'){
      expperiods_active_e <- ifelse(is.null(input$number_pieces_event_active_ae), 0, as.numeric(input$number_pieces_event_active_ae))
      ## error handling
      if(expperiods_active_e <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_active_e  <- numeric(expperiods_active_e)
      lambda_active_e <- numeric(expperiods_active_e)
      for(i in 1:expperiods_active_e){
        result <- tryCatch({
          if(i == 1){
            start_active_e[i]  <- 0
          }else{
            start_active_e[i]  <- input[[paste0('start', i, '_event_active_ae')]]
          }
          lambda_active_e[i] <- input[[paste0('lambda', i, '_event_active_ae')]]
        }, error=function(e){
          start_active_e[i] <- i - 1
          lambda_active_e[i] <- 1
        })
      }
      if(is.unsorted(start_active_e, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }

    active_ecurve <- switch(input$event_active_distribution_ae,
                            'weibull'          = {if(is.null(input$beta_event_active_ae)){temp_beta <- 1}else{temp_beta <- input$beta_event_active_ae}
                                                  if(is.null(input$alpha_event_active_ae)){temp_alpha <- 1}else{temp_alpha <- input$alpha_event_active_ae}
                                                  Weibull(alpha=temp_alpha,beta=temp_beta)
                                                  },
                            'lognormal'        = {if(is.null(input$mu_event_active_ae)){temp_mu <- 0}else{temp_mu <- input$mu_event_active_ae}
                                                  if(is.null(input$sigma_event_active_ae)){temp_sigma <- 1}else{temp_sigma <- input$sigma_event_active_ae}
                                                  Lognormal(mu=temp_mu, sigma=temp_sigma)
                                                  },
                            'exponential'      = {if(is.null(input$lambda_event_active_ae)){temp_lambda <- 0}else{temp_lambda <- input$lambda_event_active_ae}
                                                  Exponential(lambda=temp_lambda)
                                                  },
                            'pieceexponential' = PieceExponential(start=start_active_e, lambda=lambda_active_e),
                            'loglogistic'      = LogLogistic(theta=input$theta_ll_event_active_ae, eta=input$eta_ll_event_active_ae),
                            'gompertz'         = Gompertz(theta=input$theta_gomp_event_active_ae, eta=input$eta_gomp_event_active_ae),
                            'ggamma'           = GGamma(theta=ifelse(is.null(input$theta_gg_event_active_ae),1,input$theta_gg_event_active_ae), eta=ifelse(is.null(input$eta_gg_event_active_ae),1.0001,input$eta_gg_event_active_ae), rho=ifelse(is.null(input$rho_gg_event_active_ae),1,input$rho_gg_event_active_ae))
    )

    ## piecewise exponential distribution
    if(input$event_control_distribution_ae == 'pieceexponential'){
      expperiods_control_e <- ifelse(is.null(input$number_pieces_event_control_ae), 0, as.numeric(input$number_pieces_event_control_ae))
      ## error handling
      if(expperiods_control_e <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_control_e  <- numeric(expperiods_control_e)
      lambda_control_e <- numeric(expperiods_control_e)
      for(i in 1:expperiods_control_e){
        result <- tryCatch({
          if(i == 1){
            start_control_e[i]  <- 0
          }else{
            start_control_e[i]  <- input[[paste0('start', i, '_event_control_ae')]]
          }
          lambda_control_e[i] <- input[[paste0('lambda', i, '_event_control_ae')]]
        }, error=function(e){
          start_control_e[i] <- i - 1
          lambda_control_e[i] <- 1
        })
      }
      if(is.unsorted(start_control_e, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }
    control_ecurve <- switch(input$event_control_distribution_ae,
                             'weibull'          = {if(is.null(input$beta_event_control_ae)){temp_beta <- 1}else{temp_beta <- input$beta_event_control_ae}
                               if(is.null(input$alpha_event_control_ae)){temp_alpha <- 1}else{temp_alpha <- input$alpha_event_control_ae}
                               Weibull(alpha=temp_alpha, beta=temp_beta)
                             },
                             'lognormal'        = {if(is.null(input$mu_event_control_ae)){temp_mu <- 0}else{temp_mu <- input$mu_event_control_ae}
                               if(is.null(input$sigma_event_control_ae)){temp_sigma <- 1}else{temp_sigma <- input$sigma_event_control_ae}
                               Lognormal(mu=temp_mu, sigma=temp_sigma)
                             },
                             'exponential'      = {if(is.null(input$lambda_event_control_ae)){temp_lambda <- 0}else{temp_lambda <- input$lambda_event_control_ae}
                               Exponential(lambda=temp_lambda)
                             },
                             'pieceexponential' = PieceExponential(start=start_control_e, lambda=lambda_control_e),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_event_control_ae, eta=input$eta_ll_event_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_event_control_ae, eta=input$eta_gomp_event_control_ae),
                             'ggamma'           = GGamma(theta=ifelse(is.null(input$theta_gg_event_control_ae),1,input$theta_gg_event_control_ae), eta=ifelse(is.null(input$eta_gg_event_control_ae),1.0001,input$eta_gg_event_control_ae), rho=ifelse(is.null(input$rho_gg_event_control_ae),1,input$rho_gg_event_control_ae))
    )

    result <- tryCatch({
      # Plot KM curves for events, taking the CDFs from those defined within the Curve objects
      plotSF(control_ecurve, overlay=FALSE, maxT=input$assess_ae, xlab='Patient Time (months)', main='Survival Curves for \nActive (blue) and Control (black) Arms')
      plotSF(active_ecurve, overlay=TRUE, maxT=input$assess_ae, col='steelblue3')
    }, error = function(e){
      x <- seq(from=0, to=input$assess_ae, by=0.1)
      plot(x, 1 - pweibull(x, 1, 1), type='l', lwd=1.8, xlim=c(0, input$assess_ae), ylim=c(0, 1), col=1, xlab='Patient Time (months)', ylab='S(t)', main='KM Curves for \nActive (blue) and Control (black) Arms')
      lines(x, 1 - pweibull(x, 1, 1), col='steelblue3', lwd=1.8)
    })

  })

  ## plot censoring curves depending on assumptions
  output$CensCurvePlot <- renderPlot({
    ## piecewise exponential distribution
    if(input$cens_active_distribution_ae == 'pieceexponential'){
      expperiods_active_d <- ifelse(is.null(input$number_pieces_cens_active_ae), 0, as.numeric(input$number_pieces_cens_active_ae))
      ## error handling
      if(expperiods_active_d <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_active_d  <- numeric(expperiods_active_d)
      lambda_active_d <- numeric(expperiods_active_d)
      for(i in 1:expperiods_active_d){
        result <- tryCatch({
          if(i == 1){
            start_active_d[i]  <- 0
          }else{
            start_active_d[i]  <- input[[paste0('start', i, '_cens_active_ae')]]
          }
          lambda_active_d[i] <- input[[paste0('lambda', i, '_cens_active_ae')]]
        }, error=function(e){
          start_active_d[i] <- i - 1
          lambda_active_d[i] <- 1
        })
      }
      if(is.unsorted(start_active_d, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }
    active_dcurve  <- switch(input$cens_active_distribution_ae,
                             'weibull'          = {if(is.null(input$beta_cens_active_ae)){temp_beta <- 1}else{temp_beta <- input$beta_cens_active_ae}
                               if(is.null(input$alpha_cens_active_ae)){temp_alpha <- 1}else{temp_alpha <- input$alpha_cens_active_ae}
                               Weibull(alpha=temp_alpha,beta=temp_beta)
                             },
                             'lognormal'        = {if(is.null(input$mu_cens_active_ae)){temp_mu <- 0}else{temp_mu <- input$mu_cens_active_ae}
                               if(is.null(input$sigma_cens_active_ae)){temp_sigma <- 1}else{temp_sigma <- input$sigma_cens_active_ae}
                               Lognormal(mu=temp_mu, sigma=temp_sigma)
                             },
                             'exponential'      = {if(is.null(input$lambda_cens_active_ae)){temp_lambda <- 0}else{temp_lambda <- input$lambda_cens_active_ae}
                               Exponential(lambda=temp_lambda)
                             },
                             'pieceexponential' = PieceExponential(start=start_active_d, lambda=lambda_active_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_active_ae, eta=input$eta_ll_cens_active_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_active_ae, eta=input$eta_gomp_cens_active_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_active_ae, eta=input$eta_gg_cens_active_ae, rho=input$rho_gg_cens_active_ae),
                             'blank'            = Blank()
    )

    ## piecewise exponential distribution
    if(input$cens_control_distribution_ae == 'pieceexponential'){
      expperiods_control_d <- ifelse(is.null(input$number_pieces_cens_control_ae), 0, as.numeric(input$number_pieces_cens_control_ae))
      ## error handling
      if(expperiods_control_d <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_control_d  <- numeric(expperiods_control_d)
      lambda_control_d <- numeric(expperiods_control_d)
      for(i in 1:expperiods_control_d){
        result <- tryCatch({
          if(i == 1){
            start_control_d[i]  <- 0
          }else{
            start_control_d[i]  <- input[[paste0('start', i, '_cens_control_ae')]]
          }
          lambda_control_d[i] <- input[[paste0('lambda', i, '_cens_control_ae')]]
        }, error=function(e){
          start_control_d[i] <- i - 1
          lambda_control_d[i] <- 1
        })
      }
      if(is.unsorted(start_control_d, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }
    control_dcurve <- switch(input$cens_control_distribution_ae,
                             'weibull'          = {if(is.null(input$beta_cens_control_ae)){temp_beta <- 1}else{temp_beta <- input$beta_cens_control_ae}
                               if(is.null(input$alpha_cens_control_ae)){temp_alpha <- 1}else{temp_alpha <- input$alpha_cens_control_ae}
                               Weibull(alpha=temp_alpha,beta=temp_beta)
                             },
                             'lognormal'        = {if(is.null(input$mu_cens_control_ae)){temp_mu <- 0}else{temp_mu <- input$mu_cens_control_ae}
                               if(is.null(input$sigma_cens_control_ae)){temp_sigma <- 1}else{temp_sigma <- input$sigma_cens_control_ae}
                               Lognormal(mu=temp_mu, sigma=temp_sigma)
                             },
                             'exponential'      = {if(is.null(input$lambda_cens_control_ae)){temp_lambda <- 0}else{temp_lambda <- input$lambda_cens_control_ae}
                               Exponential(lambda=temp_lambda)
                             },
                             'pieceexponential' = PieceExponential(start=start_control_d, lambda=lambda_control_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_control_ae, eta=input$eta_ll_cens_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_control_ae, eta=input$eta_gomp_cens_control_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_control_ae, eta=input$eta_gg_cens_control_ae, rho=input$rho_gg_cens_control_ae),
                             'blank'            = Blank()
    )

    result <- tryCatch({
      # Plot censoring curves, taking the CDFs from those defined within the Curve objects
      plotCDF(control_dcurve, overlay=FALSE, maxT=input$assess_ae, xlab='Patient Time (months)', main='Censoring Curves for \nActive (blue) and Control (black) Arms')
      plotCDF(active_dcurve, overlay=TRUE, maxT=input$assess_ae, col='steelblue3')
    }, error = function(e){
      x <- seq(from=0, to=input$assess_ae, by=0.1)
      plot(x, pmin(x), type='l', lwd=1.8, xlim=c(0, input$assess_ae), ylim=c(0,1), xlab='Patient Time (months)', main='Censoring Curves for \nActive (blue) and Control (black) Arms')
      lines(x, pmin(x), col='steelblue3', lwd=1.8)
    })
  })

  ## plot recruitment curve depending on assumptions
  output$RecCurvePlot <- renderPlot({
    ## set distribution
    ## piecewise recruitment
    if(input$recruit_distribution_ae == 'pieceRec'){
      recruitperiods <- ifelse(is.null(input$recruitperiods_piece_ae), 0, as.numeric(input$recruitperiods_piece_ae))
      ## error handling
      if(recruitperiods <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }else{
        recruitment <- matrix(numeric(2*recruitperiods), nrow=recruitperiods)
      }
      for(i in 1:recruitperiods){
        result1 <- tryCatch({
          recruitment[i,1] <- input[[paste0('periodlength', i, '_piece_ae')]]
          recruitment[i,2] <- input[[paste0('periodrate', i, '_piece_ae')]]
        }, error = function(e){
          recruitment[1,1] <- 5
          recruitment[1,2] <- 10
        })
      }
    }
    result2 <- tryCatch({
      rcurve <- switch(input$recruit_distribution_ae,
                       'linearRec'  = LinearR(rlength=input$rlength_linear_ae, Nactive=input$Nactive_linear_ae, Ncontrol=input$Ncontrol_linear_ae),
                       'instantRec' = InstantR(Nactive=input$Nactive_instant_ae, Ncontrol=input$Ncontrol_instant_ae),
                       'pieceRec'   = PieceR(recruitment=recruitment, ratio=input$ratio_piece_ae)
      )
    }, error = function(e){
      rcurve <- LinearR(rlength=12, Nactive=30, Ncontrol=30)
    })

    result3 <- tryCatch({
      # Plot recruitment curve, taking the CDFs from those defined within the Curve objects
      plotRecruitment(rcurve, overlay=FALSE, maxT=input$assess_ae, col=1, xlab='Trial Time (months)', ylab='Patients', main='Recruitment Plot')
    }, error = function(e){
      x <- seq(from=0, to=input$assess_ae, by=0.1)
      plot(x, getPatients(InstantR(0,0),x), type='l', lwd=1.8, xlim=c(0, input$assess_ae), ylim=c(0, 60), col=1, xlab='Time (months)', ylab='Patients', main='Recruitment Plot')
    })

  })
  ## Analytical Expectations ###########################################################

  ## create reactive result to be used for Analysis
  nph_curve_trajectories_eR <- eventReactive(input$input_ae, {
    ## set distributions
    ## piecewise exponential distribution
    if(input$event_active_distribution_ae == 'pieceexponential'){
      expperiods_active_e <- ifelse(is.null(input$number_pieces_event_active_ae), 0, as.numeric(input$number_pieces_event_active_ae))
      ## error handling
      if(expperiods_active_e <= 0){
        start_active_e <- -1
        lambda_active_e <- -1
      }else{
        start_active_e  <- numeric(expperiods_active_e)
        lambda_active_e <- numeric(expperiods_active_e)
        for(i in 1:expperiods_active_e){
          result <- tryCatch({
            if(i == 1){
              start_active_e[i]  <- 0
            }else{
              start_active_e[i]  <- input[[paste0('start', i, '_event_active_ae')]]
            }
            lambda_active_e[i] <- input[[paste0('lambda', i, '_event_active_ae')]]
          }, error=function(e){
            start_active_e[i] <- i - 1
            lambda_active_e[i] <- 1
          })
        }
      }
    }

    active_ecurve <- switch(input$event_active_distribution_ae,
                            'weibull'          = Weibull(alpha=input$alpha_event_active_ae,beta=input$beta_event_active_ae),
                            'lognormal'        = Lognormal(mu=input$mu_event_active_ae, sigma=input$sigma_event_active_ae),
                            'exponential'      = Exponential(lambda=input$lambda_event_active_ae),
                            'pieceexponential' = PieceExponential(start=start_active_e, lambda=lambda_active_e),
                            'loglogistic'      = LogLogistic(theta=input$theta_ll_event_active_ae, eta=input$eta_ll_event_active_ae),
                            'gompertz'         = Gompertz(theta=input$theta_gomp_event_active_ae, eta=input$eta_gomp_event_active_ae),
                            'ggamma'           = GGamma(theta=input$theta_gg_event_active_ae, eta=input$eta_gg_event_active_ae, rho=input$rho_gg_event_active_ae)
    )
    nae <- ifelse(is.null(input$number_pieces_event_active_ae), 0, input$number_pieces_event_active_ae)

    ## piecewise exponential distribution
    if(input$event_control_distribution_ae == 'pieceexponential'){
      expperiods_control_e <- ifelse(is.null(input$number_pieces_event_control_ae), 0, as.numeric(input$number_pieces_event_control_ae))
      ## error handling
      if(expperiods_control_e <= 0){
        start_control_e <- -1
        lambda_control_e <- -1
      }else{
        start_control_e  <- numeric(expperiods_control_e)
        lambda_control_e <- numeric(expperiods_control_e)
        for(i in 1:expperiods_control_e){
          result <- tryCatch({
            if(i == 1){
              start_control_e[i]  <- 0
            }else{
              start_control_e[i]  <- input[[paste0('start', i, '_event_control_ae')]]
            }
            lambda_control_e[i] <- input[[paste0('lambda', i, '_event_control_ae')]]
          }, error=function(e){
            start_control_e[i] <- i - 1
            lambda_control_e[i] <- 1
          })
        }
      }
    }
    control_ecurve <- switch(input$event_control_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_event_control_ae,beta=input$beta_event_control_ae),
                             'lognormal'        = Lognormal(mu=input$mu_event_control_ae, sigma=input$sigma_event_control_ae),
                             'exponential'      = Exponential(lambda=input$lambda_event_control_ae),
                             'pieceexponential' = PieceExponential(start=start_control_e, lambda=lambda_control_e),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_event_control_ae, eta=input$eta_ll_event_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_event_control_ae, eta=input$eta_gomp_event_control_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_event_control_ae, eta=input$eta_gg_event_control_ae, rho=input$rho_gg_event_control_ae)
    )
    nce <- ifelse(is.null(input$number_pieces_event_control_ae), 0, input$number_pieces_event_control_ae)

    ## piecewise exponential distribution
    if(input$cens_active_distribution_ae == 'pieceexponential'){
      expperiods_active_d <- ifelse(is.null(input$number_pieces_cens_active_ae), 0, as.numeric(input$number_pieces_cens_active_ae))
      ## error handling
      if(expperiods_active_d <= 0){
        start_active_d <- -1
        lambda_active_d <- -1
      }else{
        start_active_d  <- numeric(expperiods_active_d)
        lambda_active_d <- numeric(expperiods_active_d)
        for(i in 1:expperiods_active_d){
          result <- tryCatch({
            if(i == 1){
              start_active_d[i]  <- 0
            }else{
              start_active_d[i]  <- input[[paste0('start', i, '_cens_active_ae')]]
            }
            lambda_active_d[i] <- input[[paste0('lambda', i, '_cens_active_ae')]]
          }, error=function(e){
            start_active_d[i] <- i - 1
            lambda_active_d[i] <- 1
          })
        }
      }
    }
    active_dcurve  <- switch(input$cens_active_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_cens_active_ae,beta=input$beta_cens_active_ae),
                             'lognormal'        = Lognormal(mu=input$mu_cens_active_ae, sigma=input$sigma_cens_active_ae),
                             'exponential'      = Exponential(lambda=input$lambda_cens_active_ae),
                             'pieceexponential' = PieceExponential(start=start_active_d, lambda=lambda_active_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_active_ae, eta=input$eta_ll_cens_active_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_active_ae, eta=input$eta_gomp_cens_active_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_active_ae, eta=input$eta_gg_cens_active_ae, rho=input$rho_gg_cens_active_ae),
                             'blank'            = Blank()
    )
    nad <- ifelse(is.null(input$number_pieces_cens_active_ae), 0, input$number_pieces_cens_active_ae)

    ## piecewise exponential distribution
    if(input$cens_control_distribution_ae == 'pieceexponential'){
      expperiods_control_d <- ifelse(is.null(input$number_pieces_cens_control_ae), 0, as.numeric(input$number_pieces_cens_control_ae))
      ## error handling
      if(expperiods_control_d <= 0){
        start_control_d <- -1
        lambda_control_d <- -1
      }else{
        start_control_d  <- numeric(expperiods_control_d)
        lambda_control_d <- numeric(expperiods_control_d)
        for(i in 1:expperiods_control_d){
          result <- tryCatch({
            if(i == 1){
              start_control_d[i]  <- 0
            }else{
              start_control_d[i]  <- input[[paste0('start', i, '_cens_control_ae')]]
            }
            lambda_control_d[i] <- input[[paste0('lambda', i, '_cens_control_ae')]]
          }, error=function(e){
            start_control_d[i] <- i - 1
            lambda_control_d[i] <- 1
          })
        }
      }
    }
    control_dcurve <- switch(input$cens_control_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_cens_control_ae,beta=input$beta_cens_control_ae),
                             'lognormal'        = Lognormal(mu=input$mu_cens_control_ae, sigma=input$sigma_cens_control_ae),
                             'exponential'      = Exponential(lambda=input$lambda_cens_control_ae),
                             'pieceexponential' = PieceExponential(start=start_control_d, lambda=lambda_control_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_control_ae, eta=input$eta_ll_cens_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_control_ae, eta=input$eta_gomp_cens_control_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_control_ae, eta=input$eta_gg_cens_control_ae, rho=input$rho_gg_cens_control_ae),
                             'blank'            = Blank()
    )
    ncd <- ifelse(is.null(input$number_pieces_cens_control_ae), 0, input$number_pieces_cens_control_ae)
    nec <- c(nae, nce, nad, ncd)

    ## set recruitment distribution
    ## piecewise recruitment
    if(input$recruit_distribution_ae == 'pieceRec'){
      recruitperiods <- ifelse(is.null(input$recruitperiods_piece_ae), 0, as.numeric(input$recruitperiods_piece_ae))
      ## error handling
      if(recruitperiods <= 0){
        recruitment <- matrix(c(-1,-1), nrow=2)
      }else{
        recruitment <- matrix(numeric(2*recruitperiods), nrow=recruitperiods)
      }
      for(i in 1:recruitperiods){
        recruitment[i,1] <- input[[paste0('periodlength', i, '_piece_ae')]]
        recruitment[i,2] <- input[[paste0('periodrate', i, '_piece_ae')]]
      }
    }
    rCurve <- switch(input$recruit_distribution_ae,
                     'linearRec'  = LinearR(rlength=input$rlength_linear_ae, Nactive=input$Nactive_linear_ae, Ncontrol=input$Ncontrol_linear_ae),
                     'instantRec' = InstantR(Nactive=input$Nactive_instant_ae, Ncontrol=input$Ncontrol_instant_ae),
                     'pieceRec'   = PieceR(recruitment=recruitment, ratio=input$ratio_piece_ae))
    nr <- input$recruitperiods_piece_ae

    ## set missing values
    landmarkTime    <- ifelse(is.null(input$landmarkTime_ae), Inf, input$landmarkTime_ae)
    restrictionTime <- ifelse(is.null(input$restrictionTime_ae), NA, input$restrictionTime_ae)

    landmark <- if(as.logical(input$landmark_ae)){landmarkTime} else {NULL}
    RMST <- if(as.logical(input$RMST_ae)){restrictionTime} else {NULL}

    result <- tryCatch({
    }, error=function(e){
    })

    ## conduct calculations with specified distributions

    outputtable <- tryCatch({
      nph_traj(active_ecurve=active_ecurve, control_ecurve=control_ecurve,active_dcurve=active_dcurve,control_dcurve=control_dcurve,  rcurve=rCurve, max_assessment=input$assess_ae,
                             alpha1=input$alpha1_ae, landmark=landmark, RMST=RMST,required_power=input$power_ae)$Summary
    }, error=function(e){
      NULL
    })


msg <- " "
    # currently missing input$power_ae as a final argument
    #Check for theoretical PH in Exponential/Weibull case
    PH <- FALSE
    PHa <- -1234567
    PHb <- -9876543

    if(getType(active_ecurve) == "Exponential"){
      PHa <- 1
      lambdaA <- getParams(active_ecurve)[[1]]
    } else if(getType(active_ecurve) == "Weibull"){
      PHa <- getParams(active_ecurve)[[2]]
      lambdaA <- getParams(active_ecurve)[[1]]^-getParams(active_ecurve)[[2]]
    }
    if(getType(control_ecurve) == "Exponential"){
      PHb <- 1
      lambdaC <- getParams(control_ecurve)[[1]]
    } else if(getType(control_ecurve) == "Weibull"){
      PHb <- getParams(control_ecurve)[[2]]
      lambdaC <- getParams(control_ecurve)[[1]]^-getParams(control_ecurve)[[2]]
    }
    if(PHa == PHb){
      PH <- TRUE
      HRtemp <- lambdaA/lambdaC
      msg=paste("Proportional Hazards detected: Using exact HR point estimate of",HRtemp)
    }

    ## perform checks of inputs
    if(getType(control_ecurve) != 'Lognormal' && getType(control_ecurve) != 'PieceExponential'){
      error1 <- ifelse(any(c(do.call('cbind', getParams(control_ecurve))) <= 0), 'Misspecified parameter value for event distribution (control arm).', '')
      if(getType(control_ecurve) == 'GGamma' && getParams(control_ecurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma event distribution (control arm): Must be greater than 1.'}
    }else if(getType(control_ecurve) == 'Lognormal'){
      error1 <- ifelse(getParams(control_ecurve)[[2]] <= 0, 'Misspecified parameter value for event distribution (control arm).', '')
    }else if(getType(control_ecurve) == 'PieceExponential'){
      if(nec[2] <= 0){
        error1 <- 'Invalid number of intervals for piecewiese exponential event distribution (control arm).'
      }else if(is.unsorted(getParams(control_ecurve)[[1]], strictly=TRUE)){
        error1 <- 'Start times of piecewise exponential event distribution (control arm) must be in ascending order with no duplicates'
      }else if(any(getParams(control_ecurve)[[1]] < 0) || any(getParams(control_ecurve)[[2]] <= 0)){
        error1 <- 'Start times and parameter values of piecewise exponential event distribution (control arm) must be positive.'
      }else{
        error1 <- ''
      }
    }
    if(getType(active_ecurve) != 'Lognormal' && getType(active_ecurve) != 'PieceExponential'){
      error2 <- ifelse(any(c(do.call('cbind', getParams(active_ecurve))) <= 0), 'Misspecified parameter value for event distribution (active arm).', '')
      if(getType(active_ecurve) == 'GGamma' && getParams(active_ecurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma event distribution (active arm): Must be greater than 1.'}
    }else if(getType(active_ecurve) == 'Lognormal'){
      error2 <- ifelse(getParams(active_ecurve)[[2]] <= 0, 'Misspecified parameter value for event distribution (active arm).', '')
    }else if(getType(active_ecurve) == 'PieceExponential'){
      if(nec[1] <= 0){
        error2 <- 'Invalid number of intervals for piecewise exponential event distribution (active arm).'
      }else if(is.unsorted(getParams(active_ecurve)[[1]], strictly=TRUE)){
        error2 <- 'Start times of piecewise exponential event distribution (active arm) must be in ascending order with no duplicates'
      }else if(any(getParams(active_ecurve)[[1]] < 0) || any(getParams(active_ecurve)[[2]] <= 0)){
        error2 <- 'Start times and parameter values of piecewise exponential event distribution (active arm) must be positive.'
      }else{
        error2 <- ''
      }
    }
    if(getType(control_dcurve) != 'Lognormal' && getType(control_dcurve) != 'PieceExponential' && getType(control_dcurve) != 'Blank'){
      error3 <- ifelse(any(c(do.call('cbind', getParams(control_dcurve))) <= 0), 'Misspecified parameter value for censoring distribution (control arm).', '')
      if(getType(control_dcurve) == 'GGamma' && getParams(control_dcurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma censoring distribution (control arm): Must be greater than 1.'}
    }else if(getType(control_dcurve) == 'Lognormal'){
      error3 <- ifelse(getParams(control_dcurve)[[2]] <= 0, 'Misspecified parameter value for censoring distribution (control arm).', '')
    }else if(getType(control_dcurve) == 'PieceExponential'){
      if(nec[4] <= 0){
        error3 <- 'Invalid number of intervals for piecewiese exponential censoring distribution (control arm).'
      }else if(is.unsorted(getParams(control_dcurve)[[1]], strictly=TRUE)){
        error3 <- 'Start times of piecewise exponential censoring distribution (control arm) must be in ascending order with no duplicates'
      }else if(any(getParams(control_dcurve)[[1]] < 0) || any(getParams(control_dcurve)[[2]] <= 0)){
        error3 <- 'Start times and parameter values of piecewise exponential censoring distribution (control arm) must be positive.'
      }else{
        error3 <- ''
      }
    }else{
      error3 <- ''
    }
    if(getType(active_dcurve) != 'Lognormal' && getType(active_dcurve) != 'PieceExponential' && getType(active_dcurve) != 'Blank'){
      error4 <- ifelse(any(c(do.call('cbind', getParams(active_dcurve))) <= 0), 'Misspecified parameter value for censoring distribution (active arm).', '')
      if(getType(active_dcurve) == 'GGamma' && getParams(active_dcurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma censoring distribution (active arm): Must be greater than 1.'}
    }else if(getType(active_dcurve) == 'Lognormal'){
      error4 <- ifelse(getParams(active_dcurve)[[2]] <= 0, 'Misspecified parameter value for censoring distribution (active arm).', '')
    }else if(getType(active_dcurve) == 'PieceExponential'){
      if(nec[3] <= 0){
        error4 <- 'Invalid number of intervals for piecewiese exponential censoring distribution (active arm).'
      }else if(is.unsorted(getParams(active_dcurve)[[1]], strictly=TRUE)){
        error4 <- 'Start times of piecewise exponential censoring distribution (active arm) must be in ascending order with no duplicates'
      }else if(any(getParams(active_dcurve)[[1]] < 0) || any(getParams(active_dcurve)[[2]] <= 0)){
        error4 <- 'Start times and parameter values of piecewise exponential censoring distribution (active arm) must be positive.'
      }else{
        error4 <- ''
      }
    }else{
      error4 <- ''
    }
    error5 <- ''
    if(getType(rCurve) == 'LinearR'){
      error5 <- ifelse(getRatio(rCurve) <= 0, 'Recruitment ratio must be positive.', '')
    }else if(getType(rCurve) == 'PieceR'){
      if(nr <= 0){
        error5 <- 'invalid number of intervals for piecewise recruitment.'
      }else if(any(c(do.call('cbind', getParams(rCurve))) < 0) || any(getParams(rCurve)[[2]] %% 1 != 0)){
        error5 <- 'Period length and numbers specified for piecewise recruitment must be positive. Period numbers must be integer values.'
      }
    }
    if(any(c(getNactive(rCurve), getNcontrol(rCurve)) <= 0)){
      error5 <- 'Numbers of patients specified for recruitment distribution must be positive integer values.'
    }
    error6 <- ifelse(input$alpha1_ae <= 0 | input$alpha1_ae > 0.5, 'Please specify a positive one-sided alpha less than 0.5', '')
    error7 <- ifelse(input$power_ae <= 0 | input$power_ae > 1, 'Please specify a positive power less than 1', '')
    error8 <- ifelse((input$assess_ae %% 1 != 0 | input$assess_ae < 0), 'Please specify a positive integer for the maximum assessment time', '')
    error9 <- ifelse(as.logical(input$RMST_ae) && restrictionTime < 0, 'Please specify a positive integer for the restriction time.', '')
    error10 <- ifelse(as.logical(input$landmark_ae) && landmarkTime < 0, 'Please specify a positive integer for the landmark time.', '')

    output <- list(c(error1, error2, error3, error4, error5, error6, error7, error8, error9, error10), msg, outputtable)
    return(output)
  })


  ## output reactive message to explain and warn about SS calculation
  output$SS_caveat_ae <- renderText({
    ae <- nph_curve_trajectories_eR()
    msg <- ''
    if(ae[[1]][1] == ''){msg <- 'Below is a table with estimated values for several key trial parameters. Note that the estimated required SS column provides an estimate at each assessment time of the sample size required to reach the pre-specified power if all parameters other than patient numbers are kept the same (based on Schoenfeld). This sample size should be used only as a tool to guide future runs of GESTATE.'}
    msg
  })

  ## output reactive message if PH trial is given
  output$text_ae <- renderText({
    ae <- nph_curve_trajectories_eR()
    ae[[2]]
  })

  ## output reactive analysis values table
  output$Analysis_table <- renderTable({
    ae <- nph_curve_trajectories_eR()
    ## error handling
    error <- FALSE
    for(i in 1:10){
      if(ae[[1]][i] != '') error <- TRUE
    }
    if(error){
      table <- matrix(numeric(10), nrow=10)
      for(i in 1:10){
        table[i,1] <- ae[[1]][i]
        colnames(table) <- c('ERROR')
      }
      table
      ## table output
    }else{
      table <- ae[[3]]
      out <- matrix(rep('', 11 * nrow(table)), ncol=11)
      out[,1]  <- round(table$Time, 0)
      out[,2]  <- round(table$Patients, 0)
      out[,3]  <- round(table$Events_Control, 3)
      out[,4]  <- round(table$Events_Active, 3)
      out[,5]  <- round(table$Events_Total, 3)
      out[,6]  <- round(table$HR, 4)
      out[,7]  <- round(table$LogHR, 4)
      out[,8]  <- round(table$LogHR_SE, 4)
      out[,9]  <- round(table$Schoenfeld_Power, 4)
      out[,10] <- round(table$Frontier_Power, 4)
      out[,11] <- round(table$Estimated_SS, 0)
      colnames(out) <- c('Assessment Time', 'Patients Recruited', 'Events (Control)', 'Events (Active)', 'Events (Total)', 'Hazard Ratio',
                         'Log(HR)', 'Log(HR) SE', 'Power (Schoenfeld)', 'Power (Frontier)','Estimated Required SS')

      if('RMST_Control' %in% colnames(table)){
        RMST <- matrix(rep('', 6 * nrow(table)), ncol=6)
        result1 <- tryCatch({
          RMST[,1] <- round(table$RMST_Control, 4)
          RMST[,2] <- round(table$RMST_Active, 4)
          RMST[,3] <- round(table$RMST_Delta, 4)
          RMST[,4] <- round(table$RMST_SE, 5)
          RMST[,5] <- round(table$RMST_Power, 4)
          RMST[,6] <- round(table$RMST_Failure, 4)
        }, error = function(e){

        })
        colnames(RMST) <- c('RMST (Control)', 'RMST (Active)', 'RMST Delta', 'RMST SE', 'RMST Power', 'RMST Failure')
        out <- cbind(out, RMST)
      }

      if('LM_Control' %in% colnames(table)){
        landmark <- matrix(rep('', 5 * nrow(table)), ncol=5)
        result2 <- tryCatch({
          landmark[,1] <- round(table$LM_Control, 4)
          landmark[,2] <- round(table$LM_Active, 4)
          landmark[,3] <- round(table$LM_Delta, 4)
          landmark[,4] <- round(table$LM_D_SE, 4)
          landmark[,5] <- round(table$LM_Power, 4)
        }, error = function(e){

        })
        colnames(landmark) <- c('Landmark (Control)', 'Landmark (Active)', 'Landmark Delta', 'Landmark Delta SE', 'Landmark Power')
        out <- cbind(out, landmark)
      }
      out
    }
  })

  ## plot reactive analysis results: Observed events over time
  output$ObservedEventsPlot <- renderPlot({
    ae <- nph_curve_trajectories_eR()
    ## error handling
    error <- FALSE
    for(i in 1:10){
      if(ae[[1]][i] != '') error <- TRUE
    }
    if(error){
      plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
      ## create plots for key variables against time
    }else{
      output <- ae[[3]]
      plot(output[,'Time'], output[,'Events_Total'], type='l', lwd=1.5, xlab='Trial Time (months)', ylab='Observed Events', main='Observed Events over Time')
    }
  })

  ## plot reactive analysis results: Expected log(HR) over time
  output$ExpectedLogHRPlot <- renderPlot({
    ae <- nph_curve_trajectories_eR()
    ## error handling
    error <- FALSE
    for(i in 1:10){
      if(ae[[1]][i] != '') error <- TRUE
    }
    if(error){
      plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
      ## Create plots for key variables against time
    }else{
      output <- ae[[3]]
      plot(output[,'Time'], output[,'LogHR'], type='l', lwd=1.5, xlab='Trial Time (months)', ylab='log(HR)', main='Expected log(HR) over Time', ylim=c(min(-1, output[,'LogHR']), max(1, output[,'LogHR'])))
      lines(output[,'Time'], output[,'LogHR']+qnorm(1-input$alpha1_ae)*output[,'LogHR_SE'], type='l', lwd=1.5, col='steelblue3', lty=3)
      lines(output[,'Time'], output[,'LogHR']-qnorm(1-input$alpha1_ae)*output[,'LogHR_SE'], type='l', lwd=1.5, col='steelblue3', lty=3)
      lines(output[,'Time'], output[,'Time']*0, type='l', lwd=1.5, col='darkslategray2')
    }
  })

  ## plot reactive analysis results: Power over time
  output$PowerPlot <- renderPlot({
    ae <- nph_curve_trajectories_eR()
    ## error handling
    error <- FALSE
    for(i in 1:10){
      if(ae[[1]][i] != '') error <- TRUE
    }
    if(error){
      plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
      ## Create plots for key variables against time
    }else{
      output <- ae[[3]]
      plot(output[,'Time'], output[,'Schoenfeld_Power'], type='l', lwd=1.5, xlab='Trial Time (months)', ylab='Power', main='Power over Time', ylim=c(0,1))
      lines(output[,'Time'], output[,'Frontier_Power'], type='l', lwd=1.5, col='steelblue3')
      legend('bottomright', xjust=1, yjust=1, legend=c('Schoenfeld', 'Frontier'), col=c('black', 'steelblue3'), pch=15, bty='n', horiz=FALSE, text.width=input$assess_ae*0.2)
      #text.width=input$assess_ae*0.2*c(10/29, 14/29))
    }
  })

  ## downloads #########################################################################

  ## download data (sep=',' in .csv will separate the columns for excel)
  output$downloadData_ae <- downloadHandler(
    filename = function() {paste0('analytical-expectations-', Sys.time(),'.csv')},
    content = function(file){
      ae <- nph_curve_trajectories_eR()
      ## error handling
      error <- FALSE
      for(i in 1:10){
        if(ae[[1]][i] != '') error <- TRUE
      }
      if(!error){
        table <- ae[[3]]
      }else{
        ## error output
        table <- matrix(numeric(9), nrow=9)
        for(i in 1:10){
          table[i,1] <- ae[[1]][i]
          colnames(table) <- c('ERROR')
        }
      }
      write.table(as.matrix(table), file, sep=',',row.names=FALSE)
    })

  ## download figures
  makeplot_dist <- function(){
    ## set distribution
    ## piecewise exponential distribution
    if(input$event_active_distribution_ae == 'pieceexponential'){
      expperiods_active_e <- ifelse(is.null(input$number_pieces_event_active_ae), 0, as.numeric(input$number_pieces_event_active_ae))
      ## error handling
      if(expperiods_active_e <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_active_e  <- numeric(expperiods_active_e)
      lambda_active_e <- numeric(expperiods_active_e)
      for(i in 1:expperiods_active_e){
        result <- tryCatch({
          if(i == 1){
            start_active_e[i]  <- 0
          }else{
            start_active_e[i]  <- input[[paste0('start', i, '_event_active_ae')]]
          }
          lambda_active_e[i] <- input[[paste0('lambda', i, '_event_active_ae')]]
        }, error=function(e){
          start_active_e[i] <- i - 1
          lambda_active_e[i] <- 1
        })
      }
      if(is.unsorted(start_active_e, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }
    active_ecurve <- switch(input$event_active_distribution_ae,
                            'weibull'          = Weibull(alpha=input$alpha_event_active_ae,beta=input$beta_event_active_ae),
                            'lognormal'        = Lognormal(mu=input$mu_event_active_ae, sigma=input$sigma_event_active_ae),
                            'exponential'      = Exponential(lambda=input$lambda_event_active_ae),
                            'pieceexponential' = PieceExponential(start=start_active_e, lambda=lambda_active_e),
                            'loglogistic'      = LogLogistic(theta=input$theta_ll_event_active_ae, eta=input$eta_ll_event_active_ae),
                            'gompertz'         = Gompertz(theta=input$theta_gomp_event_active_ae, eta=input$eta_gomp_event_active_ae),
                            'ggamma'           = GGamma(theta=input$theta_gg_event_active_ae, eta=input$eta_gg_event_active_ae, rho=input$rho_gg_event_active_ae)
    )

    ## piecewise exponential distribution
    if(input$event_control_distribution_ae == 'pieceexponential'){
      expperiods_control_e <- ifelse(is.null(input$number_pieces_event_control_ae), 0, as.numeric(input$number_pieces_event_control_ae))
      ## error handling
      if(expperiods_control_e <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_control_e  <- numeric(expperiods_control_e)
      lambda_control_e <- numeric(expperiods_control_e)
      for(i in 1:expperiods_control_e){
        result <- tryCatch({
          if(i == 1){
            start_control_e[i]  <- 0
          }else{
            start_control_e[i]  <- input[[paste0('start', i, '_event_control_ae')]]
          }
          lambda_control_e[i] <- input[[paste0('lambda', i, '_event_control_ae')]]
        }, error=function(e){
          start_control_e[i] <- i - 1
          lambda_control_e[i] <- 1
        })
      }
      if(is.unsorted(start_control_e, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }
    control_ecurve <- switch(input$event_control_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_event_control_ae,beta=input$beta_event_control_ae),
                             'lognormal'        = Lognormal(mu=input$mu_event_control_ae, sigma=input$sigma_event_control_ae),
                             'exponential'      = Exponential(lambda=input$lambda_event_control_ae),
                             'pieceexponential' = PieceExponential(start=start_control_e, lambda=lambda_control_e),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_event_control_ae, eta=input$eta_ll_event_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_event_control_ae, eta=input$eta_gomp_event_control_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_event_control_ae, eta=input$eta_gg_event_control_ae, rho=input$rho_gg_event_control_ae)
    )

    ## piecewise exponential distribution
    if(input$cens_active_distribution_ae == 'pieceexponential'){
      expperiods_active_d <- ifelse(is.null(input$number_pieces_cens_active_ae), 0, as.numeric(input$number_pieces_cens_active_ae))
      ## error handling
      if(expperiods_active_d <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_active_d  <- numeric(expperiods_active_d)
      lambda_active_d <- numeric(expperiods_active_d)
      for(i in 1:expperiods_active_d){
        result <- tryCatch({
          if(i == 1){
            start_active_d[i]  <- 0
          }else{
            start_active_d[i]  <- input[[paste0('start', i, '_cens_active_ae')]]
          }
          lambda_active_d[i] <- input[[paste0('lambda', i, '_cens_active_ae')]]
        }, error=function(e){
          start_active_d[i] <- i - 1
          lambda_active_d[i] <- 1
        })
      }
      if(is.unsorted(start_active_d, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }
    active_dcurve  <- switch(input$cens_active_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_cens_active_ae,beta=input$beta_cens_active_ae),
                             'lognormal'        = Lognormal(mu=input$mu_cens_active_ae, sigma=input$sigma_cens_active_ae),
                             'exponential'      = Exponential(lambda=input$lambda_cens_active_ae),
                             'pieceexponential' = PieceExponential(start=start_active_d, lambda=lambda_active_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_active_ae, eta=input$eta_ll_cens_active_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_active_ae, eta=input$eta_gomp_cens_active_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_active_ae, eta=input$eta_gg_cens_active_ae, rho=input$rho_gg_cens_active_ae),
                             'blank'            = Blank()
    )

    ## piecewise exponential distribution
    if(input$cens_control_distribution_ae == 'pieceexponential'){
      expperiods_control_d <- ifelse(is.null(input$number_pieces_cens_control_ae), 0, as.numeric(input$number_pieces_cens_control_ae))
      ## error handling
      if(expperiods_control_d <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
      start_control_d  <- numeric(expperiods_control_d)
      lambda_control_d <- numeric(expperiods_control_d)
      for(i in 1:expperiods_control_d){
        result <- tryCatch({
          if(i == 1){
            start_control_d[i]  <- 0
          }else{
            start_control_d[i]  <- input[[paste0('start', i, '_cens_control_ae')]]
          }
          lambda_control_d[i] <- input[[paste0('lambda', i, '_cens_control_ae')]]
        }, error=function(e){
          start_control_d[i] <- i - 1
          lambda_control_d[i] <- 1
        })
      }
      if(is.unsorted(start_control_d, strictly=TRUE)){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }
    }
    control_dcurve <- switch(input$cens_control_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_cens_control_ae, beta=input$beta_cens_control_ae),
                             'lognormal'        = Lognormal(mu=input$mu_cens_control_ae, sigma=input$sigma_cens_control_ae),
                             'exponential'      = Exponential(lambda=input$lambda_cens_control_ae),
                             'pieceexponential' = PieceExponential(start=start_control_d, lambda=lambda_control_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_control_ae, eta=input$eta_ll_cens_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_control_ae, eta=input$eta_gomp_cens_control_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_control_ae, eta=input$eta_gg_cens_control_ae, rho=input$rho_gg_cens_control_ae),
                             'blank'            = Blank()
    )

    result1 <- tryCatch({
      # Plot KM curves for events, taking the CDFs from those defined within the Curve objects
      plotSF(control_ecurve, overlay=FALSE, maxT=input$assess_ae, xlab='Patient Time (months)', main='Survival Curves for \nActive (blue) and Control (black) Arms')
      plotSF(active_ecurve, overlay=TRUE, maxT=input$assess_ae, col='steelblue3')
    }, error = function(e){
      x <- seq(from=0, to=input$assess_ae, by=0.1)
      plot(x, 1 - pweibull(x, 1, 1), type='l', lwd=1.8, xlim=c(0, input$assess_ae), ylim=c(0, 1), col=1, xlab='Patient Time (months)', ylab='S(t)', main='KM Curves for \nActive (blue) and Control (black) Arms')
      lines(x, 1 - pweibull(x, 1, 1), col='steelblue3', lwd=1.8)
    })

    result2 <- tryCatch({
      # Plot censoring curves, taking the CDFs from those defined within the Curve objects
      plotCDF(control_dcurve, overlay=FALSE, maxT=input$assess_ae, xlab='Patient Time (months)', main='Censoring Curves for \nActive (blue) and Control (black) Arms')
      plotCDF(active_dcurve, overlay=TRUE, maxT=input$assess_ae, col='steelblue3')
    }, error = function(e){
      x <- seq(from=0, to=input$assess_ae, by=0.1)
      plot(x, pmin(x), type='l', lwd=1.8, xlim=c(0, input$assess_ae), ylim=c(0,1), col=1, xlab='Patient Time (months)', ylab='P(t)', main='Censoring Curves for \nActive (blue) and Control (black) Arms')
      lines(x, pmin(x), col='steelblue3', lwd=1.8)
    })

    if(input$recruit_distribution_ae == 'pieceRec'){
      recruitperiods <- ifelse(is.null(input$recruitperiods_piece_ae), 0, as.numeric(input$recruitperiods_piece_ae))
      ## error handling
      if(recruitperiods <= 0){
        plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
        return()
      }else{
        recruitment <- matrix(numeric(2*recruitperiods), nrow=recruitperiods)
      }
      for(i in 1:recruitperiods){
        result3 <- tryCatch({
          recruitment[i,1] <- input[[paste0('periodlength', i, '_piece_ae')]]
          recruitment[i,2] <- input[[paste0('periodrate', i, '_piece_ae')]]
        }, error = function(e){
          recruitment[1,1] <- 5
          recruitment[1,2] <- 10
        })
      }
    }
    result4 <- tryCatch({
      rcurve <- switch(input$recruit_distribution_ae,
                       'linearRec'  = LinearR(rlength=input$rlength_linear_ae, Nactive=input$Nactive_linear_ae, Ncontrol=input$Ncontrol_linear_ae),
                       'instantRec' = InstantR(Nactive=input$Nactive_instant_ae, Ncontrol=input$Ncontrol_instant_ae),
                       'pieceRec'   = PieceR(recruitment=recruitment, ratio=input$ratio_piece_ae)
      )
    }, error = function(e){
      rcurve <- LinearR(rlength=input$rlength_linear_ae, Nactive=30, Ncontrol=30)
    })

    result5 <- tryCatch({
      # Plot recruitment curve, taking the CDFs from those defined within the Curve objects
      plotRecruitment(rcurve, overlay=FALSE, maxT=input$assess_ae, xlab='Time (months)', ylab='Patients', main='Recruitment Plot')
    }, error = function(e){
      x <- seq(from=0, to=input$assess_ae, by=0.1)
      plot(x, getPatients(InstantR(0,0),x), type='l', lwd=1.8, xlim=c(0, input$assess_ae), ylim=c(0, 60), xlab='Time (months)', ylab='Patients', main='Recruitment Plot')
    })
  }

  output$downloadPlots_dist <- downloadHandler(
    filename = function() {paste0('distribution-plots-', Sys.time(), '.pdf')},
    content = function(file_dist_out){
      pdf(file_dist_out)
      makeplot_dist()
      dev.off()
    }
  )

  makeplot_ae <- function(){
    ae <- nph_curve_trajectories_eR()
    ## error handling
    error <- FALSE
    for(i in 1:10){
      if(ae[[1]][i] != '') error <- TRUE
    }
    if(error){
      plot(0, 0, type='n', xlim=c(-1, 1), ylim=c(0, 0), xaxt='n', yaxt='n', xlab=' ', ylab=' ', main=' ')
      ## create plots for key variables against time
    }else{
      output <- ae[[3]]
      plot(output[,'Time'], output[,'Events_Total'], type='l', lwd=1.5, xlab='Trial Time (months)', ylab='Observed Events', main='Observed Events over Time')
      plot(output[,'Time'], output[,'LogHR'], type='l', lwd=1.5, xlab='Trial Time (months)', ylab='log(HR)', main='Expected log(HR) over Time', ylim=c(min(-1, output[,'LogHR']), max(1, output[,'LogHR'])))
      lines(output[,'Time'], output[,'LogHR']+qnorm(1-input$alpha1_ae)*output[,'LogHR_SE'], type='l', lwd=1.5, col='steelblue3', lty=3)
      lines(output[,'Time'], output[,'LogHR']-qnorm(1-input$alpha1_ae)*output[,'LogHR_SE'], type='l', lwd=1.5, col='steelblue3', lty=3)
      lines(output[,'Time'], output[,'Time']*0, type='l', lwd=1.5, col='darkslategray2')
      plot(output[,'Time'], output[,'Schoenfeld_Power'], type='l', lwd=1.5, xlab='Trial Time (months)', ylab='Power', main='Power over Time', ylim=c(0,1))
      lines(output[,'Time'], output[,'Frontier_Power'], type='l', lwd=1.5, col='steelblue3')
      legend('bottomright', xjust=1, yjust=1, legend=c('Schoenfeld', 'Frontier'), col=c('black', 'steelblue3'), pch=15, bty='n', horiz=FALSE, text.width=input$assess_ae*0.2)

    }
  }

  output$downloadPlots_ae <- downloadHandler(
    filename = function() {paste0('analytical-expectations-plots-', Sys.time(), '.pdf')},
    content = function(file_ae_out){
      pdf(file_ae_out)
      makeplot_ae()
      dev.off()
    }
  )

  ## Simulation ########################################################################

  ## create reactive simulation
  simulate_trials_eR <- eventReactive(input$input_sim, {
    ## set distributions
    ## piecewise exponential distribution
    if(input$event_active_distribution_ae == 'pieceexponential'){
      expperiods_active_e <- ifelse(is.null(input$number_pieces_event_active_ae), 0, as.numeric(input$number_pieces_event_active_ae))
      ## error handling
      if(expperiods_active_e <= 0){
        start_active_e <- -1
        lambda_active_e <- -1
      }else{
        start_active_e  <- numeric(expperiods_active_e)
        lambda_active_e <- numeric(expperiods_active_e)
        for(i in 1:expperiods_active_e){
          result <- tryCatch({
            if(i == 1){
              start_active_e[i]  <- 0
            }else{
              start_active_e[i]  <- input[[paste0('start', i, '_event_active_ae')]]
            }
            lambda_active_e[i] <- input[[paste0('lambda', i, '_event_active_ae')]]
          }, error=function(e){
            start_active_e[i] <- i - 1
            lambda_active_e[i] <- 1
          })
        }
      }
    }
    active_ecurve <- switch(input$event_active_distribution_ae,
                            'weibull'          = Weibull(alpha=input$alpha_event_active_ae, beta=input$beta_event_active_ae),
                            'lognormal'        = Lognormal(mu=input$mu_event_active_ae, sigma=input$sigma_event_active_ae),
                            'exponential'      = Exponential(lambda=input$lambda_event_active_ae),
                            'pieceexponential' = PieceExponential(start=start_active_e, lambda=lambda_active_e),
                            'loglogistic'      = LogLogistic(theta=input$theta_ll_event_active_ae, eta=input$eta_ll_event_active_ae),
                            'gompertz'         = Gompertz(theta=input$theta_gomp_event_active_ae, eta=input$eta_gomp_event_active_ae),
                            'ggamma'           = GGamma(theta=input$theta_gg_event_active_ae, eta=input$eta_gg_event_active_ae, rho=input$rho_gg_event_active_ae)
    )
    nae <- ifelse(is.null(input$number_pieces_event_active_ae), 0, input$number_pieces_event_active_ae)

    ## piecewise exponential distribution
    if(input$event_control_distribution_ae == 'pieceexponential'){
      expperiods_control_e <- ifelse(is.null(input$number_pieces_event_control_ae), 0, as.numeric(input$number_pieces_event_control_ae))
      ## error handling
      if(expperiods_control_e <= 0){
        start_control_e <- -1
        lambda_control_e <- -1
      }else{
        start_control_e  <- numeric(expperiods_control_e)
        lambda_control_e <- numeric(expperiods_control_e)
        for(i in 1:expperiods_control_e){
          result <- tryCatch({
            if(i == 1){
              start_control_e[i]  <- 0
            }else{
              start_control_e[i]  <- input[[paste0('start', i, '_event_control_ae')]]
            }
            lambda_control_e[i] <- input[[paste0('lambda', i, '_event_control_ae')]]
          }, error=function(e){
            start_control_e[i] <- i - 1
            lambda_control_e[i] <- 1
          })
        }
      }
    }
    control_ecurve <- switch(input$event_control_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_event_control_ae,beta=input$beta_event_control_ae),
                             'lognormal'        = Lognormal(mu=input$mu_event_control_ae, sigma=input$sigma_event_control_ae),
                             'exponential'      = Exponential(lambda=input$lambda_event_control_ae),
                             'pieceexponential' = PieceExponential(start=start_control_e, lambda=lambda_control_e),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_event_control_ae, eta=input$eta_ll_event_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_event_control_ae, eta=input$eta_gomp_event_control_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_event_control_ae, eta=input$eta_gg_event_control_ae, rho=input$rho_gg_event_control_ae)
    )
    nce <- ifelse(is.null(input$number_pieces_event_control_ae), 0, input$number_pieces_event_control_ae)

    ## piecewise exponential distribution
    if(input$cens_active_distribution_ae == 'pieceexponential'){
      expperiods_active_d <- ifelse(is.null(input$number_pieces_cens_active_ae), 0, as.numeric(input$number_pieces_cens_active_ae))
      ## error handling
      if(expperiods_active_d <= 0){
        start_active_d <- -1
        lambda_active_d <- -1
      }else{
        start_active_d  <- numeric(expperiods_active_d)
        lambda_active_d <- numeric(expperiods_active_d)
        for(i in 1:expperiods_active_d){
          result <- tryCatch({
            if(i == 1){
              start_active_d[i]  <- 0
            }else{
              start_active_d[i]  <- input[[paste0('start', i, '_cens_active_ae')]]
            }
            lambda_active_d[i] <- input[[paste0('lambda', i, '_cens_active_ae')]]
          }, error=function(e){
            start_active_d[i] <- i - 1
            lambda_active_d[i] <- 1
          })
        }
      }
    }
    active_dcurve  <- switch(input$cens_active_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_cens_active_ae,beta=input$beta_cens_active_ae),
                             'lognormal'        = Lognormal(mu=input$mu_cens_active_ae, sigma=input$sigma_cens_active_ae),
                             'exponential'      = Exponential(lambda=input$lambda_cens_active_ae),
                             'pieceexponential' = PieceExponential(start=start_active_d, lambda=lambda_active_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_active_ae, eta=input$eta_ll_cens_active_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_active_ae, eta=input$eta_gomp_cens_active_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_active_ae, eta=input$eta_gg_cens_active_ae, rho=input$rho_gg_cens_active_ae),
                             'blank'            = Blank()
    )
    nad <- ifelse(is.null(input$number_pieces_cens_active_ae), 0, input$number_pieces_cens_active_ae)

    ## piecewise exponential distribution
    if(input$cens_control_distribution_ae == 'pieceexponential'){
      expperiods_control_d <- ifelse(is.null(input$number_pieces_cens_control_ae), 0, as.numeric(input$number_pieces_cens_control_ae))
      ## error handling
      if(expperiods_control_d <= 0){
        start_control_d <- -1
        lambda_control_d <- -1
      }else{
        start_control_d  <- numeric(expperiods_control_d)
        lambda_control_d <- numeric(expperiods_control_d)
        for(i in 1:expperiods_control_d){
          result <- tryCatch({
            if(i == 1){
              start_control_d[i]  <- 0
            }else{
              start_control_d[i]  <- input[[paste0('start', i, '_cens_control_ae')]]
            }
            lambda_control_d[i] <- input[[paste0('lambda', i, '_cens_control_ae')]]
          }, error=function(e){
            start_control_d[i] <- i - 1
            lambda_control_d[i] <- 1
          })
        }
      }
    }
    control_dcurve <- switch(input$cens_control_distribution_ae,
                             'weibull'          = Weibull(alpha=input$alpha_cens_control_ae,beta=input$beta_cens_control_ae),
                             'lognormal'        = Lognormal(mu=input$mu_cens_control_ae, sigma=input$sigma_cens_control_ae),
                             'exponential'      = Exponential(lambda=input$lambda_cens_control_ae),
                             'pieceexponential' = PieceExponential(start=start_control_d, lambda=lambda_control_d),
                             'loglogistic'      = LogLogistic(theta=input$theta_ll_cens_control_ae, eta=input$eta_ll_cens_control_ae),
                             'gompertz'         = Gompertz(theta=input$theta_gomp_cens_control_ae, eta=input$eta_gomp_cens_control_ae),
                             'ggamma'           = GGamma(theta=input$theta_gg_cens_control_ae, eta=input$eta_gg_cens_control_ae, rho=input$rho_gg_cens_control_ae),
                             'blank'            = Blank()
    )
    ncd <- ifelse(is.null(input$number_pieces_cens_control_ae), 0, input$number_pieces_cens_control_ae)
    nec <- c(nae, nce, nad, ncd)

    ## set recruitment distribution
    ## piecewise recruitment
    if(input$recruit_distribution_ae == 'pieceRec'){
      recruitperiods <- ifelse(is.null(input$recruitperiods_piece_ae), 0, as.numeric(input$recruitperiods_piece_ae))
      ## error handling
      if(recruitperiods <= 0){
        recruitment <- matrix(c(-1,-1), nrow=2)
      }else{
        recruitment <- matrix(numeric(2*recruitperiods), nrow=recruitperiods)
      }
      for(i in 1:recruitperiods){
        recruitment[i,1] <- input[[paste0('periodlength', i, '_piece_ae')]]
        recruitment[i,2] <- input[[paste0('periodrate', i, '_piece_ae')]]
      }
    }
    rCurve <- switch(input$recruit_distribution_ae,
                     'linearRec'  = LinearR(rlength=input$rlength_linear_ae, Nactive=input$Nactive_linear_ae, Ncontrol=input$Ncontrol_linear_ae),
                     'instantRec' = InstantR(Nactive=input$Nactive_instant_ae, Ncontrol=input$Ncontrol_instant_ae),
                     'pieceRec'   = PieceR(recruitment=recruitment, ratio=input$ratio_piece_ae))
    nr <- input$recruitperiods_piece_ae

    ## carry out simulation with specified distributions

    ## fix event number if selected
    fixevents <- NULL
    if(as.logical(input$create_events_sim)){
      if(input$force_events_sim %% 1 != 0 || input$force_events_sim < 0){
        output[[1]][8] <- 'Number of events to fix in each simulation must be a positive integer less than or equal to the number of patients.'
        error <- TRUE
        fixevents <- NULL
      }else{
        fixevents <- input$force_events_sim
      }
    }

    sim <- tryCatch({
        simulate_trials(active_ecurve=active_ecurve, control_ecurve=control_ecurve, active_dcurve=active_dcurve,control_dcurve=control_dcurve,rcurve=rCurve,assess=input$assess_sim_ae,fix_events=fixevents, iterations=input$iterations_sim, seed=input$seed_sim)
      }, error=function(e){
      NULL
    })

    ## perform checks of inputs
    if(getType(control_ecurve) != 'Lognormal' && getType(control_ecurve) != 'PieceExponential'){
      error1 <- ifelse(any(c(do.call('cbind', getParams(control_ecurve))) <= 0), 'Misspecified parameter value for event distribution (control arm).', '')
      if(getType(control_ecurve) == 'GGamma' && getParams(control_ecurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma event distribution (control arm): Must be greater than 1.'}
    }else if(getType(control_ecurve) == 'Lognormal'){
      error1 <- ifelse(getParams(control_ecurve)[[2]] <= 0, 'Misspecified parameter value for event distribution (control arm).', '')
    }else if(getType(control_ecurve) == 'PieceExponential'){
      if(nec[2] <= 0){
        error1 <- 'Invalid number of intervals for piecewiese exponential event distribution (control arm).'
      }else if(is.unsorted(getParams(control_ecurve)[[1]], strictly=TRUE)){
        error1 <- 'Start times of piecewise exponential event distribution (control arm) must be in ascending order with no duplicates'
      }else if(any(getParams(control_ecurve)[[1]] < 0) || any(getParams(control_ecurve)[[2]] <= 0)){
        error1 <- 'Start times and parameter values of piecewise exponential event distribution (control arm) must be positive.'
      }else{
        error1 <- ''
      }
    }
    if(getType(active_ecurve) != 'Lognormal' && getType(active_ecurve) != 'PieceExponential'){
      error2 <- ifelse(any(c(do.call('cbind', getParams(active_ecurve))) <= 0), 'Misspecified parameter value for event distribution (active arm).', '')
      if(getType(active_ecurve) == 'GGamma' && getParams(active_ecurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma event distribution (active arm): Must be greater than 1.'}
    }else if(getType(active_ecurve) == 'Lognormal'){
      error2 <- ifelse(getParams(active_ecurve)[[2]] <= 0, 'Misspecified parameter value for event distribution (active arm).', '')
    }else if(getType(active_ecurve) == 'PieceExponential'){
      if(nec[1] <= 0){
        error2 <- 'Invalid number of intervals for piecewiese exponential event distribution (active arm).'
      }else if(is.unsorted(getParams(active_ecurve)[[1]], strictly=TRUE)){
        error2 <- 'Start times of piecewise exponential event distribution (active arm) must be in ascending order with no duplicates'
      }else if(any(getParams(active_ecurve)[[1]] < 0) || any(getParams(active_ecurve)[[2]] <= 0)){
        error2 <- 'Start times and parameter values of piecewise exponential event distribution (active arm) must be positive.'
      }else{
        error2 <- ''
      }
    }
    if(getType(control_dcurve) != 'Lognormal' && getType(control_dcurve) != 'PieceExponential' && getType(control_dcurve) != 'Blank'){
      error3 <- ifelse(any(c(do.call('cbind', getParams(control_dcurve))) <= 0), 'Misspecified parameter value for censoring distribution (control arm).', '')
      if(getType(control_dcurve) == 'GGamma' && getParams(control_dcurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma censoring distribution (control arm): Must be greater than 1.'}
    }else if(getType(control_dcurve) == 'Lognormal'){
      error3 <- ifelse(getParams(control_dcurve)[[2]] <= 0, 'Misspecified parameter value for censoring distribution (control arm).', '')
    }else if(getType(control_dcurve) == 'PieceExponential'){
      if(nec[4] <= 0){
        error3 <- 'Invalid number of intervals for piecewiese exponential censoring distribution (control arm).'
      }else if(is.unsorted(getParams(control_dcurve)[[1]], strictly=TRUE)){
        error3 <- 'Start times of piecewise exponential censoring distribution (control arm) must be in ascending order with no duplicates'
      }else if(any(getParams(control_dcurve)[[1]] < 0) || any(getParams(control_dcurve)[[2]] <= 0)){
        error3 <- 'Start times and parameter values of piecewise exponential censoring distribution (control arm) must be positive.'
      }else{
        error3 <- ''
      }
    }else{
      error3 <- ''
    }
    if(getType(active_dcurve) != 'Lognormal' && getType(active_dcurve) != 'PieceExponential' && getType(active_dcurve) != 'Blank'){
      error4 <- ifelse(any(c(do.call('cbind', getParams(active_dcurve))) <= 0), 'Misspecified parameter value for censoring distribution (active arm).', '')
      if(getType(active_dcurve) == 'GGamma' && getParams(active_dcurve)[[2]] <= 1){error1 <- 'Misspecified shape parameter for the generalised gamma censoring distribution (active arm): Must be greater than 1.'}
    }else if(getType(active_dcurve) == 'Lognormal'){
      error4 <- ifelse(getParams(active_dcurve)[[2]] <= 0, 'Misspecified parameter value for censoring distribution (active arm).', '')
    }else if(getType(active_dcurve) == 'PieceExponential'){
      if(nec[3] <= 0){
        error4 <- 'Invalid number of intervals for piecewiese exponential censoring distribution (active arm).'
      }else if(is.unsorted(getParams(active_dcurve)[[1]], strictly=TRUE)){
        error4 <- 'Start times of piecewise exponential censoring distribution (active arm) must be in ascending order with no duplicates'
      }else if(any(getParams(active_dcurve)[[1]] < 0) || any(getParams(active_dcurve)[[2]] <= 0)){
        error4 <- 'Start times and parameter values of piecewise exponential censoring distribution (active arm) must be positive.'
      }else{
        error4 <- ''
      }
    }else{
      error4 <- ''
    }
    error5 <- ''
    if(getType(rCurve) == 'LinearR'){
      error5 <- ifelse(getRatio(rCurve) < 0, 'Recruitment ratio must be positive.', '')
    }else if(getType(rCurve) == 'PieceR'){
      if(nr <= 0){
        error5 <- 'invalid number of intervals for piecewise recruitment.'
      }else if(any(c(do.call('cbind', getParams(rCurve))) < 0) || any(getParams(rCurve)[[2]] %% 1 != 0)){
        error5 <- 'Period length and numbers specified for piecewise recruitment must be positive. Period numbers must be integer values.'
      }
    }
    if(any(c(getNactive(rCurve), getNcontrol(rCurve)) <= 0)){
      error5 <- 'Numbers of patients specified for recruitment distribution must be positive integer values.'
    }
    error6 <- ifelse((input$iterations_sim %% 1 != 0 | input$iterations_sim < 0), 'Please specify a positive integer for the number of simulations', '')
    error7 <- ifelse((input$assess_ae %% 1 != 0 | input$assess_ae < 0), 'Please specify a positive integer for the assessment time', '')

    output <- list(c(error1, error2, error3, error4, error5, error6, error7, ''), sim)


    ## error handling
    error <- FALSE
    for(i in 1:7){
      if(output[[1]][i] != '') {error <- TRUE}
    }



    return(output)
  })

  ## output reactive values table based on simulation
  output$Simulation_table <- renderTable({
    sim <- simulate_trials_eR()
    ## error handling
    error <- FALSE
    for(i in 1:8){
      if(sim[[1]][i] != '') error <- TRUE
    }
    if(error){
      table <- matrix(numeric(8), nrow=8)
      for(i in 1:8){
        table[i,1] <- sim[[1]][i]
        colnames(table) <- c('ERROR')
      }
      table
      ## table output
    }else{
      table <- sim[[2]]

      ## handle display options

      disp <- ifelse(is.null(input$displayData_sim) | is.na(input$displayData_sim), 1, input$displayData_sim)
      disp <- round(disp, 0)
      if(disp < 1){disp <- 1}
      if(disp > input$iterations_sim){disp <- input$iterations_sim}
      l <- nrow(table[table[,'Iter'] <= disp,])
      width <- 4
      out <- matrix(rep('', width * l), ncol=width)
      out[,1] <- round(table[table[,'Iter'] <= disp, 'Time'], 4)
      out[,2] <- round(table[table[,'Iter'] <= disp, 'Censored'], 0)
      out[,2] <- ifelse(out[,2] == '1', 'Yes', 'No')
      out[,3] <- round(table[table[,'Iter'] <= disp, 'Trt'], 0)
      out[,3] <- ifelse(out[,3] == '2', 'Active', 'Control')
      out[,4] <- round(table[table[,'Iter'] <= disp, 'Iter'], 0)
      colnames(out) <- c('Time','Censored','Treatment','Simulation')
      out
    }
  })

  ## create reactive analysis based on simulation values
  sim_analysis_eR <- eventReactive(input$input_sim, {
    sim <- simulate_trials_eR()
    ## error handling
    error <- FALSE
    for(i in 1:8){
      if(sim[[1]][i] != '') error <- TRUE
    }
    if(error){
      table <- matrix(c(''), ncol=1)
      colnames(table) <- c('ERROR')
      ## table
    }else{
      sim <- sim[[2]]
      ## conduct analysis by parallel processing if selected
      ## set missing values
      cores            <- ifelse(is.null(input$cores_sim), 1, input$cores_sim)
      #RMST            <- ifelse(is.null(input$restrictionTime_ae), NULL, input$restrictionTime_ae)
      if(is.null(input$restrictionTime_ae)){RMST <- NULL} else{RMST <- input$restrictionTime_ae}
      #landmark        <- ifelse(is.null(input$landmarkTime_ae), NULL, input$landmarkTime_ae)
      if(is.null(input$restrictionTime_ae)){landmark <- NULL} else{landmark <- input$restrictionTime_ae}
      alpha2           <- ifelse(is.null(input$alpha1_ae), 0.05, 2 * input$alpha1_ae)
      ##error handling
      error1 <- ifelse(as.logical(input$parallel_analysis_sim) && (cores %% 1 != 0 || cores < 1), 'Number of cores for parallel processing must be a positive integer.', '')
      error2 <- ifelse(as.logical(input$landmark_ae) && landmark < 0, 'Time for Landmark calculation must be positive.', '')
      error3 <- ifelse(as.logical(input$RMST_ae) && RMST < 0, 'Restriction time must be positive.', '')

      if(any(c(error1, error2, error3) != '')){
        table <- matrix(c(error1, error2, error3), nrow=3)
        colnames(table) <- c('ERROR')
      }else{
        table <- analyse_sim(data=sim, LR=as.logical(input$lrcox_sim),RMST=RMST,landmark=landmark,stratum="", parallel_cores=cores)
        Sim_Number <- c(1:nrow(table))
        table <- cbind(Sim_Number,table)
      }
    }

    return(table)
  })

  ## output reactive analysis table based on simulation values
  output$Simulation_analysis <- renderTable({
    table <- sim_analysis_eR()
    if(length(table[1,] > 1)){
      ## handle display options
      disp <- ifelse(is.null(input$displayData_sim) || is.na(input$displayData_sim), 1, input$displayData_sim)
      disp <- round(disp, 0)
      if(disp < 1){disp <- 1}
      if(disp > input$iterations_sim){disp <- input$iterations_sim}
      if('LogHR' %in% colnames(table)){
        res <- matrix(rep('', 9 * disp), ncol=9)
        res[,1] <- round(table[1:disp, 'LogHR'], 4)
        res[,2] <- round(table[1:disp, 'HR'], 4)
        res[,3] <- round(table[1:disp, 'LogHR_SE'], 5)
        res[,4] <- round(table[1:disp, 'HR_Z'], 4)
        res[,5] <- round(table[1:disp, 'LR_Z'], 4)
        res[,6] <- round(table[1:disp, 'HR_P'], 4)
        res[,7] <- round(table[1:disp, 'LR_P'], 4)
        res[,8] <- round(table[1:disp, 'Events_Active'], 4)
        res[,9] <- round(table[1:disp, 'Events_Control'], 4)

        colnames(res) <- c('Log(HR)', 'HR', 'Log(HR) SE', 'HR Z-Score', 'Log-Rank Z-Score', 'HR p-Value', 'Log-Rank p-Value',
                           'Observed Events (Active)', 'Observed Events (Control)')
      }
      if('RMST_Control' %in% colnames(table)){
        RMST <- matrix(rep('', 6 * disp), ncol=6)
        RMST[,1] <- round(table[1:disp, 'RMST_Active'], 4)
        RMST[,2] <- round(table[1:disp, 'RMST_A_SE'], 5)
        RMST[,3] <- round(table[1:disp, 'RMST_Control'], 4)
        RMST[,4] <- round(table[1:disp, 'RMST_C_SE'], 5)
        RMST[,5] <- round(table[1:disp, 'RMST_Delta'], 4)
        RMST[,6] <- round(table[1:disp, 'RMST_D_SE'], 5)
        colnames(RMST) <- c('RMST (Active)', 'RMST (Active) SE', 'RMST (Control)', 'RMST (Control) SE', 'RMST Difference',
                            'RMST Difference SE')
        res <- if('LogHR' %in% colnames(table)){cbind(res, RMST)} else {RMST}
      }

      if('LM_Active' %in% colnames(table)){
        LM <- matrix(rep('', 8 * disp), ncol=8)
        LM[,1] <- round(table[1:disp, 'LM_Active'], 4)
        LM[,2] <- round(table[1:disp, 'LM_A_SE'], 5)
        LM[,3] <- round(table[1:disp, 'LM_Control'], 4)
        LM[,4] <- round(table[1:disp, 'LM_C_SE'], 5)
        LM[,5] <- round(table[1:disp, 'LM_Delta'], 4)
        LM[,6] <- round(table[1:disp, 'LM_D_SE'], 5)
        LM[,7] <- round(table[1:disp, 'LM_Z'], 4)
        LM[,8] <- round(table[1:disp, 'LM_P'], 4)
        colnames(LM) <- c('Landmark Survival (Active)', 'Landmark Survival (Active) SE', 'Landmark Survival (Control)', 'Landmark Survival (Control) SE',
                          'Landmark Survival Difference', 'Landmark Survival Difference SE', 'Landmark Z-Score', 'Landmark p-Value')
        res <- if('LogHR'%in% colnames(table) | 'RMST_Control' %in% colnames(table)){cbind(res, LM)} else{LM}
      }
      table <- res
    }
    table
  })

  ## output summary of analysis based on simulation values
  output$summarised_analysis <- renderTable({
    ana <- sim_analysis_eR()
    if(length(ana[1,]) > 1){ # if an error has occured, the table has only one column and can't be summarised
      ## set missing values
      alpha1 <- ifelse(is.null(input$alpha1_ae), 0.025, input$alpha1_ae)
      table <- summarise_analysis(ana, alpha1)
        if('LogHR' %in% colnames(table)){
        res <- matrix(rep('', 13), ncol=13)
        res[,1]  <- round(table[,'LogHR'], 4)
        res[,2]  <- round(table[,'HR'], 4)
        res[,3]  <- round(table[,'LogHR_SE'], 5)
        res[,4]  <- round(table[,'HR_Z'], 4)
        res[,5]  <- round(table[,'HR_P'], 4)
        res[,6]  <- round(table[,'LR_Z'], 4)
        res[,7]  <- round(table[,'LR_P'], 4)
        res[,8]  <- round(table[,'Events_Active'], 4)
        res[,9]  <- round(table[,'Events_Control'], 4)
        res[,10] <- round(table[,'Events_Total'], 4)
        res[,11] <- round(table[,'LR_Failed'], 4)
        res[,12] <- round(table[,'LR_Power'], 4)
        res[,13] <- round(table[,'Simulations'], 0)
        colnames(res) <- c('Log(HR)', 'HR', 'Log(HR) SE', 'Cox Z-value', 'Cox P-Value', 'Log-Rank Z-Value', 'Log-Rank P-Value',
                           'Observed Events (Active)', 'Observed Events (Control)', 'Observed Events (Total)', 'Failed (Log-Rank)', 'Log-Rank Power','Simulations')
      }

      if('RMST_Control' %in% colnames(table)){
        RMST <- matrix(rep('', 8), ncol=8)
        RMST[,1]  <- round(table[,'RMST_Active'], 4)
        RMST[,2]  <- round(table[,'RMST_A_SE'], 5)
        RMST[,3]  <- round(table[,'RMST_Control'], 4)
        RMST[,4]  <- round(table[,'RMST_C_SE'], 5)
        RMST[,5]  <- round(table[,'RMST_Delta'], 4)
        RMST[,6]  <- round(table[,'RMST_D_SE'], 5)
        RMST[,7]  <- round(table[,'RMST_Failed'], 4)
        RMST[,8]  <- round(table[,'RMST_Power'], 4)
        colnames(RMST) <- c('RMST (Active)', 'RMST (Active) SE', 'RMST (Control)', 'RMST (Control) SE', 'RMST Difference', 'RMST Difference SE',
                            'Failed (RMST)', 'RMST Power')
        res <- if('LogHR' %in% colnames(table)){cbind(res, RMST)} else{RMST}
      }

      if('LM_Active' %in% colnames(table)){
        LM <- matrix(rep('', 8 ), ncol=8)
        LM[,1] <- round(table[,'LM_Active'], 4)
        LM[,2] <- round(table[,'LM_A_SE'], 5)
        LM[,3] <- round(table[,'LM_Control'], 4)
        LM[,4] <- round(table[,'LM_C_SE'], 5)
        LM[,5] <- round(table[,'LM_Delta'], 4)
        LM[,6] <- round(table[,'LM_D_SE'], 5)
        LM[,7] <- round(table[,'LM_Failed'], 4)
        LM[,8] <- round(table[,'LM_Power'], 4)
        colnames(LM) <- c('Landmark Survival (Active)', 'Landmark Survival (Active) SE', 'Landmark Survival (Control)', 'Landmark Survival (Control) SE',
                          'Landmark Survival Difference', 'Landmark Survival Delta SE', 'Failed (Landmark)', 'Landmark Power')
        res <- if('LogHR'%in% colnames(table) | 'RMST_Control' %in% colnames(table)){cbind(res, LM)} else{LM}
      }
      table <- res

    }else{
      table <- matrix(c(''), ncol=1)
      colnames(table) <- c('ERROR')
    }
    table
  })

  ## downloads #########################################################################

  ## download data (sep=',' in .csv will separate the columns for excel)
  output$downloadData_sim <- downloadHandler(
    filename = function() {paste0('simulation-data-', Sys.time(),'.csv')},
    content = function(file){
      sim <- simulate_trials_eR()
      error <- FALSE
      for(i in 1:8){
        if(sim[[1]][i] != '') error <- TRUE
      }
      if(error){
        table <- matrix(numeric(8), nrow=8)
        for(i in 1:8){
          table[i,1] <- sim[[1]][i]
          colnames(table) <- c('ERROR')
        }
        ## table output
      }else{
        table <- sim[[2]]
      }
      write.table(as.matrix(table), file, sep=',',row.names=FALSE)
    })

  ## download data (sep=',' in .csv will separate the columns for excel)
  output$downloadAna_sim <- downloadHandler(
    filename = function() {paste0('simulation-analysis-', Sys.time(),'.csv')},
    content = function(file){
      table <- sim_analysis_eR()
      write.table(as.matrix(table), file, sep=',',row.names=FALSE)
    })

  ## download data (sep=',' in .csv will separate the columns for excel)
  output$downloadSum_sim <- downloadHandler(
    filename = function() {paste0('simulation-summary-', Sys.time(),'.csv')},
    content = function(file){
      ana <- sim_analysis_eR()
      if(length(ana[1,]) > 1){ # if an error has occured, the table has only one column and can't be summarised
        ## set missing values
        alpha1 <- ifelse(is.null(input$alpha1_ae), 0.025, input$alpha1_ae)
        table <- summarise_analysis(ana, alpha1)
      }else{
        table <- matrix(c(''), ncol=1)
        colnames(table) <- c('ERROR')
      }
      write.table(as.matrix(table), file, sep=',',row.names=FALSE)
    })
}
