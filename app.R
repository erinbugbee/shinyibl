# Install Packages (Run these lines in this console once, removing the #)
#install.packages("shiny")
#install.packages("ggvis")
#install.packages("dplyr")
#install.packages("reshape2")
#install.packages("tidyr")
#install.packages("rhandsontable")
#install.packages("data.table")

# Load Packages
library(shiny)
library(ggvis)
library(dplyr)
library(reshape2)
library(tidyr)
library(rhandsontable)
library(data.table)
library(shinylive)
library(ggplot2)

# Initial values and names for input table
init <- data.table(V1 = c(0, 10),
                   V2 = c(100, 20),
                   p1 = c(.1, .5))
init$p2 <- 1 - init$p1
init$ev <- init$V1 * init$p1 + init$V2 * init$p2
rows <- c("A", "B")
name_out <- c("O1", "O2", "p(O1)", "p(O2)", "EV")
name_in <- c("V1", "V2", "p1", "p2", "ev")
rownames(init) <- rows


# Define UI
ui <- fluidPage(
  titlePanel("ShinyIBL"),
  # following fluidRow implements the loading message while simulation is running
  fluidRow(
    tags$style(type="text/css", "
               #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
               }
               "),
    ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Instructions:"),
      tags$h6("This is a GUI to simulate choices from the instance-based learning model. To use it:"),
      tags$h6("1. Enter the binary choice problem that you want to simulate."),
      tags$h6("2. Define the simulation settings"),

      tags$h6("3. Define the IBL model parameters."),
      tags$h6("4. Run the simulation."),
      h4("Define Gamble Values:"),
      rHandsontableOutput("hot"),
      h4("Define Simulation Settings:"),
      sliderInput("subj",
                  "Number of Subjects:",
                  min = 1,
                  max = 200,
                  value = 10),
      sliderInput("trial",
                  "Number of Trials:",
                  min = 5,
                  max = 300,
                  value = 20),
      h4("Define Model Parameters:"),
      sliderInput("decay",
                  "Value of Decay Parameter:",
                  min = -5,
                  max = 5,
                  value = 0.75,
                  step = .25),
      sliderInput("sigma",
                  "Value of Noise Parameter:",
                  min = 0.1,
                  max = 5,
                  value = 0.5),
      p(actionButton("go", "Run Simulation", icon("play"))),

     # tags$h5("Select which plots you would like to see:"),
    #  checkboxInput('checkp', 'p(A) Plot'),
    #  checkboxInput('checkbv', 'Blended Values Plot'),
    #  checkboxInput('checkprob', 'Probability of Retrieval Plot'),
    #  checkboxInput('checkact', 'Activation Plot')
    ),
    mainPanel(
      fluidRow(
        column(6,
               tags$div(style = "border: 2px solid #ddd; padding: 15px; border-radius: 10px; background-color: #f9f9f9; margin-bottom: 20px;",
                        tags$h4("Simulation Settings", style = "font-weight: bold; color: #333;"),
                        tags$p("Number of Subjects: ", strong(textOutput("subj_var", inline = TRUE))),
                        tags$p("Number of Trials: ", strong(textOutput("trial_var", inline = TRUE)))
               )
        ),
        column(6,
               tags$div(style = "border: 2px solid #ddd; padding: 15px; border-radius: 10px; background-color: #f9f9f9; margin-bottom: 20px;",
                        tags$h4("Model Parameters", style = "font-weight: bold; color: #333;"),
                        tags$p("Decay Parameter: ", strong(textOutput("decay_var", inline = TRUE))),
                        tags$p("Noise Parameter: ", strong(textOutput("noise_var", inline = TRUE)))
               )
        )
      ),
      
      fluidRow(
        column(6, plotOutput("iblPlot")),
        column(6, plotOutput("bvPlot"))
      ),
      fluidRow(
        column(6, plotOutput("probPlot")),
        column(6, plotOutput("actsPlot"))
      )
    )
  )
)

# Define server
server <- function(input, output) {
  v <- reactiveValues(data = NULL)
  
  # IBL calculation function
  iblCalc <- function(idx, v_a1, v_a2, v_b1, v_b2, pa, pb, decay, trial, sigma, timer) {
    # get IBl predictions for current gamble/participant sampling period
    tau <- sigma * sqrt(2)
    # set up memory structure and sets pre-populated values as 10% larger than largest outcome
    pval <- c(v_a1, v_a2, v_b1, v_b2)
    pp <- max(pval[c(pa, 1 - pa, pb, 1 - pb) > 0])
    pp <- pp * 1.1
    mem <- setNames(list(c(), c(), c(), c(), c(), c()), 
                    c(paste0("1_1_", as.character(v_a1)),
                      paste0("1_2_", as.character(v_a2)),
                      as.character(pp),
                      paste0("2_1_", as.character(v_b1)),
                      paste0("2_2_", as.character(v_b2)),
                      as.character(pp)))
    mem[names(mem) == pp] <- 1

    ns <- matrix(runif((trial+1) * 6, min = 0.1), ncol = 6) # add uniform samples to generate logistic noise
    
    ### prespecify output data structure
    out <- rep(NA, trial+1) # simplex probability of choosing option a
    act_out <- matrix(rep(-10, 4 * (trial+1)), ncol = 4) # activations
    bv_out <- matrix(rep(NA, 2 * (trial+1)), ncol = 2) # blended values
    pr_out <- matrix(rep(NA, 4 * (trial+1)), ncol = 4) # probability of recovery
    out[1] <- 1 / 2
    act_out[1,] <- -10
    for (t in 2:(trial+1)) {
      # calculate activations (decayed traces plus noise)
      acts <- sapply(1:6, function(x) log(sum((t - mem[[x]]) ^ (-decay), na.rm = TRUE))) + sigma * log((1 -  ns[t, ]) / ns[t, ])
      acts <- ifelse(acts == -Inf, -10, acts)
     # acts <- ifelse(c(pa, (1 - pa), 1, pb, (1 - pb), 1) == 0, -Inf, acts)
      act_out[t,] <- acts[c(1:2, 4:5)] # only save non-prepopulated value activations
      # calculate cognitive probabilities
      ps <- c(exp(acts[1:3] / tau) / sum(exp(acts[1:3] / tau)),  exp(acts[4:6] / tau) / sum(exp(acts[4:6] / tau)))
      pr_out[t,] <- ps[c(1:2,4:5)] # only save non-prepopulated value probabilities
      
      # calculate option-wise blended values
      vals <- c(ps[1:3] %*% c(v_a1, v_a2, pp), ps[4:6] %*% c(v_b1, v_b2, pp))
      bv_out[t,] <- vals

      out[t] <- mean(vals[1] > vals[2]) # Calculate P(A) as proportion when BV of A > BV of B

      # Fixes error with NAs
      safe_softmax <- function(x) {
        x <- x - max(x, na.rm = TRUE)  # Normalize for numerical stability
        exp_x <- exp(x)
        exp_x / sum(exp_x, na.rm = TRUE)
      }
      
      probs <- safe_softmax(vals)
      
      # play relevant gamble and record observed outcome into memory
      tmp <- ifelse(
        sample(1:2, 1, prob = probs) == 1, 
        as.character(sample(paste0("1_", c(paste0("1_", v_a1), paste0("2_", v_a2))), 1, prob = c(pa, 1 - pa))), 
        as.character(sample(paste0("2_", c(paste0("1_", v_b1), paste0("2_", v_b2))), 1, prob = c(pb, 1 - pb)))
      )
      mem[names(mem) == tmp][[1]] <- c(mem[names(mem) == tmp][[1]], t)
    }
    tmp <- data.frame(cbind(idx, out, v_a1, v_a2, v_b1, v_b2, pa, pb, act_out, bv_out, 1:(trial+1), pr_out, max(trial)))
    names(tmp) <- c("idx", "out", "v_a1", "v_a2", "v_b1", "v_b2", "pa", "pb", "aa_1", "aa_2", "ab_1", "ab_2", "bv_a", "bv_b", "trial", "pa_1", "pa_2", "pb_1", "pb_2", "maxTrial")
    tmp
  }
  

  values <- reactiveValues(hot = init) # table setup based on initial values at top
  
  # Initialized values for first plot
  output$subj_var <- renderText({ 
    isolate(input$subj)
  })
  
  output$trial_var <- renderText({ 
    isolate(input$trial)
  })
  
  output$decay_var <- renderText({ 
    isolate(input$decay)
  })
  
  output$noise_var <- renderText({ 
    isolate(input$sigma)
  })
  
  # Update values when Run Simulation is clicked
  observeEvent(input$go, {
    output$subj_var <- renderText({ 
      isolate(input$subj)
    })
    
    output$trial_var <- renderText({ 
      isolate(input$trial)
    })
    
    output$decay_var <- renderText({ 
      isolate(input$decay)
    })
    
    output$noise_var <- renderText({ 
      isolate(input$sigma)
    })
    
    })
  
  
  # following chunk handles manual table updating
  output$hot = renderRHandsontable({
    DT = NULL
    if (!is.null(input$hot)) {
      DT = setDT(hot_to_r(input$hot))
      values[["hot"]] = DT
    } else if (!is.null(values[["hot"]])) {
      DT = values[["hot"]]
    } 
    if (!is.null(DT)){
      names(DT) <- name_in
    DT$p1 <- ifelse(DT$p1 > 1 | DT$p1 < 0, .5, DT$p1)
    DT$p2 <- 1 - DT$p1
    DT$ev <- DT$V1 * DT$p1 + DT$V2 * DT$p2
    names(DT) <- name_out
    rownames(DT) <- rows

    }

    rhandsontable(DT) %>%
      hot_col(col = "O1") %>%
      hot_col(col = "O2") %>%
      hot_col(col = "p(O1)") %>%
      hot_col(col = "p(O2)", readOnly = T) %>%
      hot_col(col = "EV", readOnly = T)
    
  })
  
  ibl <- reactiveValues(plota = NULL)
  iblreset <- reactiveValues(plota = NULL)

  # p_dat contains the main simulation data
  p_dat <- reactive({
    input$go
    isolate(port <- values[["hot"]])
    isolate(names(port) <- name_in)
    isolate(rownames(port) <- rows)
    isolate(do.call("rbind", lapply(1:input$subj, iblCalc, port$V1[1], port$V2[1], port$V1[2], port$V2[2], port$p1[1], port$p1[2], input$decay, input$trial, input$sigma, input$subj
    )))
  })
  
  # Probability of Choosing Option Plot
  output$iblPlot <- renderPlot({
    raw_data <- p_dat() %>%
      dplyr::filter(trial > 1) %>%
      mutate(pB = 1 - out) %>%  # Compute p(B)
      pivot_longer(cols = c(out, pB), names_to = "Option", values_to = "Probability") %>%
      mutate(Option = recode(Option, "out" = "A", "pB" = "B"))  # Rename legend labels
    
    avg_data <- raw_data %>%
      group_by(Option, trial) %>%
      summarise(Probability = mean(Probability), .groups = "drop")  # Compute averages
    
    ggplot() +
      geom_point(data = raw_data, aes(x = trial - 1, y = Probability, color = Option),
                 alpha = 0.1, size = 2) +  
      
      geom_line(data = avg_data, aes(x = trial - 1, y = Probability, color = Option),
                linewidth = 1) +  
      
      ylim(0, 1) +
      labs(x = "Trial", y = "Probability", title = "Probability of Choice", color = "Option") +  
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "left",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
      ) +
      annotate("text", x = max(raw_data$trial) - 2, y = 0.75, 
               label = paste0("p(A): ", sprintf("%.2f", mean(raw_data$Probability[raw_data$Option == 'A'])), 
                              "\n p(B): ", sprintf("%.2f", mean(raw_data$Probability[raw_data$Option == 'B']))),
               hjust = 1, vjust = 1, size = 6, color = "black", fontface = "bold")
  })
  

  # Aggregate data for Blended Value Plot
  bv_dat <- reactive({
    p_dat() %>%
      dplyr::select(bv_a, bv_b, idx, trial) %>%
      dplyr::filter(trial > 1) %>%
      pivot_longer(cols = c("bv_a", "bv_b"), names_to = "Option", values_to = "bv") %>%
      mutate(Option = ifelse(Option == "bv_a", "A", "B"))
  })
  
  # Blended Value Plot
  output$bvPlot <- renderPlot({
    raw_data <- bv_dat()
    
    avg_data <- raw_data %>%
      group_by(Option, trial) %>%
      summarise(bv = mean(bv), .groups = "drop")
    
    ggplot() +
      geom_point(data = raw_data, aes(x = trial - 1, y = bv, color = Option),
                 alpha = 0.1, size = 2) +  
      
      geom_line(data = avg_data, aes(x = trial - 1, y = bv, color = Option),
                linewidth = 1) +  
      
      labs(x = "Trial", y = "Blended Value", title = "Blended Values") +
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "left",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
      )
  })
  
  
  # Aggregate data for Probability of Retrieval Plot
  pr_dat <- reactive({
    p_dat() %>%
      dplyr::select(idx, pa_1:pb_2, trial, v_a1, v_a2, v_b1, v_b2, pa, pb) %>%
      dplyr::filter(trial > 1) %>%
      pivot_longer(cols = c(pa_1, pa_2, pb_1, pb_2), names_to = "opts", values_to = "acts") %>%
      filter(acts > -Inf) %>%
      mutate(opts = case_when(
        opts == "pa_1" ~ paste("Option A, Outcome: ", v_a1, ", Probability: ", pa, sep = ""),
        opts == "pa_2" ~ paste("Option A, Outcome: ", v_a2, ", Probability: ", 1 - pa, sep = ""),
        opts == "pb_1" ~ paste("Option B, Outcome: ", v_b1, ", Probability: ", pb, sep = ""),
        opts == "pb_2" ~ paste("Option B, Outcome: ", v_b2, ", Probability: ", 1 - pb, sep = "")
      ))
  })
  
  # Probability of Retrieval Plot
  output$probPlot <- renderPlot({
    raw_data <- pr_dat()
    
    avg_data <- raw_data %>%
      group_by(opts, trial) %>%
      summarise(acts = mean(acts), .groups = "drop")  # Compute averages
    
    ggplot() +
      geom_point(data = raw_data, aes(x = trial - 1, y = acts, color = opts),
                 alpha = 0.1, size = 2) +  
      
      geom_line(data = avg_data, aes(x = trial - 1, y = acts, color = opts),
                linewidth = 1) +  
      
      labs(x = "Trial", y = "Probability of Retrieval", title = "Probability of Retrieval", color = "Option and Outcome/Probability") +
      scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +  
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "left",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
      )
  })
  
  
  # Aggregate data for Activation Plot
  a_dat <- reactive({
    p_dat() %>%
      dplyr::select(idx, aa_1:ab_2, trial, v_a1, v_a2, v_b1, v_b2, pa, pb) %>%
      dplyr::filter(trial > 1) %>%
      pivot_longer(cols = c(aa_1, aa_2, ab_1, ab_2), names_to = "opts", values_to = "acts") %>%
      filter(acts > -Inf) %>%
      mutate(opts = case_when(
        opts == "aa_1" ~ paste("Option A, Outcome: ", v_a1, ", Probability: ", pa, sep = ""),
        opts == "aa_2" ~ paste("Option A, Outcome: ", v_a2, ", Probability: ", 1 - pa, sep = ""),
        opts == "ab_1" ~ paste("Option B, Outcome: ", v_b1, ", Probability: ", pb, sep = ""),
        opts == "ab_2" ~ paste("Option B, Outcome: ", v_b2, ", Probability: ", 1 - pb, sep = "")
      ))
  })
  
  
    # Activation Plot
  output$actsPlot <- renderPlot({
    raw_data <- a_dat()  # Get raw data (all points)
    
    avg_data <- raw_data %>%
      group_by(opts, trial) %>%
      summarise(acts = mean(acts), .groups = "drop")  # Compute averages
    
    ggplot() +
      geom_point(data = raw_data, aes(x = trial - 1, y = acts, color = opts),
                 alpha = 0.1, size = 2) +  
      
      geom_line(data = avg_data, aes(x = trial - 1, y = acts, color = opts),
                size = 1) +  
      
      labs(x = "Trial", y = "Activation", title = "Activation", color = "Option and Outcome/Probability") +
      scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +  
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "left",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
      )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)