library(shiny)
library(ggvis)
library(dplyr)
library(reshape2)
library(tidyr)
library(rhandsontable)
library(data.table)

#cat("Got here", file=stderr())

init <- data.table(V1 = c(0, 10),
                   V2 = c(100, 20),
                   p1 = c(.1, .5))
init$p2 <- 1 - init$p1
init$ev <- init$V1 * init$p1 + init$V2 * init$p2
rows <- c("A", "B")
name_out <- c("O1", "O2", "p(O1)", "p(O2)", "EV")
name_in <- c("V1", "V2", "p1", "p2", "ev")
rownames(init) <- rows

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Shiny IBL"),
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
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div("Loading...",id="loadmessage")
    )
    ),
  sidebarBottomPage(
    sidebarBottomPanel(
      h4("Instructions:"),
      tags$h6("This is a GUI to simulate choices from the instance-based learning model. To use it:"),
      tags$h6("1. Enter the binary choice problem that you want to simulate"),
      tags$h6("2. Define the simulation parameters"),
      tags$h6("3. Define the IBL model parameters"),
      tags$h6("4. Run simulation"),
      h4("Define gamble values:"),
      rHandsontableOutput("hot"),
      h4("Define simulation settings:"),
      sliderInput("subj",
                  "Number of subjects:",
                  min = 1,
                  max = 200,
                  value = 10),
      sliderInput("trial",
                  "Maximum number of trials:",
                  min = 3,
                  max = 300,
                  value = 20),
      h4("Define model parameters:"),
      sliderInput("decay",
                  "Value of decay parameter:",
                  min = -5,
                  max = 5,
                  value = 0.75,
                  step = .25),
      sliderInput("sigma",
                  "Value of noise parameter:",
                  min = 0.1,
                  max = 5,
                  value = 0.5),
      p(actionButton("go", "Run simulation", icon("random")))
    ),
    mainTopPanel(
      fluidRow(
        column(10, ggvisOutput("iblPlot"))
      ),
      fluidRow(
        column(10, ggvisOutput("bvPlot"))
      ),
      fluidRow(
        column(12, ggvisOutput("actsPlot"))
      ),
      fluidRow(
        column(12, ggvisOutput("probPlot"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # IBL calculation
  iblCalc <- function(idx, v_a1, v_a2, v_b1, v_b2, pa, pb, decay, trial, sigma, timer) {
    # get IBl predictions for current gamble/participant sampling period
    tau <- sigma * sqrt(2)
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
    ns <- matrix(runif(trial * 6), ncol = 6)
    # ns <- t(apply(ns, 1, function(x) x * (c(pa, (1 - pa), 1, pb, (1 - pb), 1) != 0)))
    out <- rep(NA, trial) # simplex probability of choosing option a
    act_out <- matrix(rep(NA, 4 * trial), ncol = 4)
    bv_out <- matrix(rep(NA, 2 * trial), ncol = 2)
    pr_out <- matrix(rep(NA, 4 * trial), ncol = 4)
    out[1] <- 1 / 2
    for (t in 2:trial) {
      noise <- ns[t, ]
      acts <- sapply(1:6, function(x) log(sum((t - mem[[x]]) ^ -decay, na.rm = TRUE))) + sigma * log((1 -  noise) / noise)
      # acts <- ifelse(acts == Inf, -Inf, acts)
      acts <- ifelse(c(pa, (1 - pa), 1, pb, (1 - pb), 1) == 0, -Inf, acts)
      act_out[t,] <- acts[c(1:2, 4:5)]
      
      ps <- c(exp(acts[1:3] / tau) / sum(exp(acts[1:3] / tau)),  exp(acts[4:6] / tau) / sum(exp(acts[4:6] / tau)))
      pr_out[t,] <- ps[c(1:2,4:5)]
      
      vals <- c(ps[1:3] %*% c(v_a1, v_a2, pp), ps[4:6] %*% c(v_b1, v_b2, pp))
      bv_out[t,] <- vals
      
      out[t] <- exp(vals[1]) / sum(exp(vals))
      tmp <- ifelse(
        sample(1:2, 1, prob = exp(exp(vals) / sum(exp(vals)))) == 1, 
        as.character(sample(paste0("1_", c(paste0("1_", v_a1), paste0("2_", v_a2))), 1, prob = c(pa, 1 - pa))), 
        as.character(sample(paste0("2_", c(paste0("1_", v_b1), paste0("2_", v_b2))), 1, prob = c(pb, 1 - pb)))
      )
      mem[names(mem) == tmp][[1]] <- c(mem[names(mem) == tmp][[1]], t)
    }
    tmp <- data.frame(cbind(idx, out, act_out, bv_out, 1:trial, pr_out))
    names(tmp) <- c("idx", "out", "aa_1", "aa_2", "ab_1", "ab_2", "bv_a", "bv_b", "trial", "pa_1", "pa_2", "pb_1", "pb_2")
    tmp
  }
  
  values <- reactiveValues(hot = init) # table setup based on initial values at top
  
  output$hot = renderRHandsontable({
    DT = NULL
    if (!is.null(input$hot)) {
      DT = setDT(hot_to_r(input$hot))
      values[["hot"]] = DT
    } else if (!is.null(values[["hot"]])) {
      DT = values[["hot"]]
    }
    if (!is.null(DT))
      names(DT) <- name_in
    DT$p1 <- ifelse(DT$p1 > 1 | DT$p1 < 0, .5, DT$p1)
    DT$p2 <- 1 - DT$p1
    DT$ev <- DT$V1 * DT$p1 + DT$V2 * DT$p2
    names(DT) <- name_out
    rownames(DT) <- rows
    rhandsontable(DT) %>%
      hot_col(col = "O1") %>%
      hot_col(col = "O2") %>%
      hot_col(col = "p(O1)") %>%
      hot_col(col = "p(O2)", readOnly = TRUE) %>%
      hot_col(col = "EV", readOnly = TRUE, renderer = "function(instance, td, row, col,          prop, value, cellProperties) {
              Handsontable.TextCell.renderer.apply(this, arguments);
              td.style.background = 'lightblue';
  }")
      
})
  
  p_dat <- reactive({
    input$go
    isolate(port <- values[["hot"]])
    isolate(names(port) <- name_in)
    isolate(rownames(port) <- rows)
    isolate(do.call("rbind", lapply(1:input$subj, iblCalc, port$V1[1], port$V2[1], port$V1[2], port$V2[2], port$p1[1], port$p1[2], input$decay, input$trial, input$sigma, input$subj
    )))
  })
  
  p_dat %>% ggvis(~trial, ~out) %>%
    layer_smooths(stroke := "orange", fill := "orange", se = TRUE) %>%
    group_by(trial) %>%
    summarise(out = mean(out)) %>%
    layer_lines() %>%
    add_axis("y", title = "p(A)") %>%
    add_axis("x", title = "Trial") %>%
    scale_numeric("y", domain = c(0, 1), nice = TRUE) %>%
    summarise(est   = paste("p(A) =", round(mean(out), 2)),
              trial = max(trial) + 1) %>%
    mutate(out = 0.5) %>%
    layer_text(text := ~est) %>%
    bind_shiny("iblPlot")
  
  bv_dat <- reactive({
    p_dat() %>%
      dplyr::select(bv_a, bv_b, idx, trial) %>%
      dplyr::filter(trial > 1) %>%
      gather(opts, bv, -trial, -idx) %>%
      mutate(Option = ifelse(opts == "bv_a", "A", "B"))
  })
  
  a_dat <- reactive({
    p_dat() %>%
      dplyr::select(idx, aa_1:ab_2, trial) %>%
      dplyr::filter(trial > 1) %>%
      gather(opts, acts, -trial, -idx) %>%
      filter(acts > -Inf) %>%
      mutate(opts = ifelse(opts == "aa_1", "Option A, Outcome 1", ifelse(opts == "aa_2", "Option A, Outcome 2", ifelse(opts == "ab_1", "Option B, Outcome 1", "Option B, Outcome 2"))))
  })
  
  pr_dat <- reactive({
    p_dat() %>%
      dplyr::select(idx, pa_1:pb_2, trial) %>%
      dplyr::filter(trial > 1) %>%
      gather(opts, acts, -trial, -idx) %>%
      filter(acts > -Inf) %>%
      mutate(opts = ifelse(opts == "pa_1", "Option A, Outcome 1", ifelse(opts == "pa_2", "Option A, Outcome 2", ifelse(opts == "pb_1", "Option B, Outcome 1", "Option B, Outcome 2"))))
  })
  
  bv_dat %>% ggvis(~trial, ~bv) %>%
    group_by(Option) %>%
    layer_smooths(stroke = ~Option, se = TRUE) %>%
    group_by(Option, trial) %>%
    summarise(bv = mean(bv)) %>%
    layer_lines() %>%
    # layer_points(fill = ~Option, opacity := .5) %>%
    add_axis("y", title = "Blended Value") %>%
    add_axis("x", title = "Trial") %>%
    bind_shiny("bvPlot")
  
  a_dat %>% ggvis(~trial, ~acts) %>%
    group_by(opts) %>%
    layer_smooths(stroke = ~opts, fill = ~opts, se = TRUE) %>%
    # layer_points(fill = ~opts, opacity := .5) %>%
    add_axis("y", title = "Activation") %>%
    add_axis("x", title = "Trial") %>%
    bind_shiny("actsPlot")
  
  pr_dat %>% ggvis(~trial, ~acts) %>%
    group_by(opts) %>%
    layer_smooths(stroke = ~opts, fill = ~opts, se = TRUE) %>%
    # layer_points(fill = ~opts, opacity := .5) %>%
    add_axis("y", title = "Probability of Retrieval") %>%
    add_axis("x", title = "Trial") %>%
    scale_numeric("y", domain = c(0, 1), nice = TRUE) %>%
    bind_shiny("probPlot")
}

# Run the application 
shinyApp(ui = ui, server = server)

