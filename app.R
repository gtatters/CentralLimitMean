# Central Limit Theorem for Means - Base R only

#shinylive::export(appdir = "../CentralLimitMean/", destdir = "docs")
#httpuv::runStaticServer("docs/", port = 8008)

# Define UI --------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Central Limit Theorem for Means", windowTitle = "CLT for means"),
  
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        radioButtons("dist", "Parent distribution (population):",
                     c("Normal"       = "rnorm",
                       "Uniform"      = "runif",
                       "Right skewed" = "rlnorm",
                       "Left skewed"  = "rbeta"),
                     selected = "rnorm"),
        
        uiOutput("mu"),
        uiOutput("sd"),
        uiOutput("minmax"),
        uiOutput("skew"),
        
        sliderInput("n", "Sample size:",   value = 30,  min = 2,   max = 500),
        br(),
        sliderInput("k", "Number of samples:", value = 200, min = 10, max = 1000)
      ),
      
      helpText("Glenn Tattersall, PhD"),
      helpText("For use in BIOL 3P96 - Biostatistics")
    ),
    
    mainPanel(
      tabsetPanel(
        type = "tabs",
        
        tabPanel(
          title = "Population Distribution",
          plotOutput("pop.dist", height = "500px"),
          br()
        ),
        
        tabPanel(
          title = "Samples",
          br(),
          plotOutput("sample.dist", height = "600px"),
          div(h3(textOutput("num.samples")), align = "center"),
          br()
        ),
        
        tabPanel(
          title = "Sampling Distribution",
          fluidRow(
            column(width = 7,
                   br(), br(),
                   div(textOutput("CLT.descr"), align = "justify")),
            column(width = 5,
                   br(),
                   plotOutput("pop.dist.two", width = "85%", height = "200px"))
          ),
          fluidRow(
            column(width = 12,
                   br(),
                   plotOutput("sampling.dist"),
                   div(textOutput("sampling.descr", inline = TRUE), align = "center"),
                   br(),
                   # --- NEW: normality diagnostic panel ---
                   wellPanel(
                     div(strong("Normality Check"), align = "center"),
                     br(),
                     div(textOutput("normality.descr"), align = "justify")
                   ))
          )
        )
      )
    )
  )
)

# Define server ----------------------------------------------------------------

seed <- as.numeric(Sys.time())

server <- function(input, output, session) {
  
  # Dynamic UI: mean slider (Normal) ----
  output$mu <- renderUI({
    if (input$dist == "rnorm")
      sliderInput("mu", "Mean:", value = 0, min = -40, max = 50)
  })
  
  # Dynamic UI: SD slider (Normal) ----
  output$sd <- renderUI({
    if (input$dist == "rnorm")
      sliderInput("sd", "Standard deviation:", value = 20, min = 1, max = 30)
  })
  
  # Dynamic UI: min/max slider (Uniform) ----
  output$minmax <- renderUI({
    if (input$dist == "runif")
      sliderInput("minmax", "Lower and Upper Bounds",
                  value = c(5, 15), min = 0, max = 20)
  })
  
  # Guard: uniform range cannot be zero ----
  observeEvent(input$minmax, {
    req(input$minmax)
    if (input$minmax[1] == input$minmax[2]) {
      if      (input$minmax[1] == 0)  updateSliderInput(session, "minmax", value = c(0, 1))
      else if (input$minmax[2] == 20) updateSliderInput(session, "minmax", value = c(19, 20))
      else    updateSliderInput(session, "minmax", value = c(input$minmax[2], input$minmax[2] + 1))
    }
  })
  
  # Dynamic UI: skew selector (LogNormal / Beta) ----
  output$skew <- renderUI({
    if (input$dist %in% c("rlnorm", "rbeta"))
      selectInput("skew", "Skew:",
                  choices  = c("Low skew" = "low", "Medium skew" = "med", "High skew" = "high"),
                  selected = "low")
  })
  
  # Helper: draw n values from the chosen distribution ----
  rand_draw <- function(dist, n, mu, sd, min, max, skew) {
    switch(dist,
           rbeta = {
             req(skew)
             params <- switch(skew,
                              low  = list(n = n, shape1 = 5, shape2 = 2),
                              med  = list(n = n, shape1 = 5, shape2 = 1.5),
                              high = list(n = n, shape1 = 5, shape2 = 1))
             do.call(rbeta, params)
           },
           rnorm = {
             req(mu, sd)
             rnorm(n, mean = mu, sd = sd)
           },
           rlnorm = {
             req(skew)
             sdlog <- switch(skew, low = 0.25, med = 0.5, high = 1)
             rlnorm(n, meanlog = 0, sdlog = sdlog)
           },
           runif = {
             req(min, max)
             runif(n, min = min, max = max)
           }
    )
  }
  
  rep_rand_draw <- repeatable(rand_draw)
  
  # Reactive: large population vector ----
  parent <- reactive({
    rep_rand_draw(input$dist, 1e5, input$mu, input$sd,
                  input$minmax[1], input$minmax[2], input$skew)
  })
  
  # Reactive: matrix of k samples, each of size n ----
  samples <- reactive({
    pop <- parent()
    replicate(input$k, sample(pop, input$n, replace = TRUE))
  })
  
  # Reactive: vector of sample means ----
  sample_means <- reactive({
    colMeans(samples())
  })
  
  u_min <- reactive({ req(input$minmax); input$minmax[1] })
  u_max <- reactive({ req(input$minmax); input$minmax[2] })
  
  # Helper: base-R histogram + density with annotation ----
  plot_dist <- function(x, title, col_fill, x_pos, y_pos, m, s,
                        x_lim = NULL, cex_base = 1) {
    
    d   <- density(x)
    brk <- pretty(x, n = 45)
    
    if (is.null(x_lim)) x_lim <- range(x)
    
    hist(x, breaks = brk, freq = FALSE,
         col  = col_fill, border = "white",
         main = title, xlab = "x", ylab = "Density",
         xlim = x_lim,
         cex.main = cex_base * 1.1,
         cex.axis = cex_base,
         cex.lab  = cex_base)
    
    lines(d, col = col_fill, lwd = 2)
    
    label <- paste0("mean of x = ", m, "\nSD of x = ", s)
    text(x_pos, y_pos, labels = label, col = "black", cex = cex_base * 0.9)
  }
  
  # Shared logic for pop.dist and pop.dist.two ----
  make_pop_plot <- function(cex_base = 1.4) {
    
    distname <- switch(input$dist,
                       rnorm  = "Population distribution: Normal",
                       rlnorm = "Population distribution: Right skewed",
                       rbeta  = "Population distribution: Left skewed",
                       runif  = "Population distribution: Uniform")
    
    pop    <- parent()
    m_pop  <- round(mean(pop), 2)
    sd_pop <- round(sd(pop), 2)
    d      <- density(pop)
    x_rng  <- diff(range(pop))
    y_pos  <- max(d$y) - 0.2 * max(d$y)
    
    if (input$dist == "rnorm") {
      req(input$mu)
      x_lim <- c(min(-100, min(pop)), max(100, max(pop)))
      x_pos <- ifelse(input$mu > 0,
                      min(-100, min(pop)) + 0.15 * diff(x_lim),
                      max(100,  max(pop)) - 0.15 * diff(x_lim))
      plot_dist(pop, distname, "#195190", x_pos, y_pos,
                m_pop, sd_pop, x_lim = x_lim, cex_base = cex_base)
      
    } else if (input$dist == "runif") {
      if (u_min() == u_max()) return(invisible(NULL))
      x_pos <- max(pop) - 0.1 * x_rng
      y_pos2 <- max(d$y) * 1.3
      plot_dist(pop, distname, "#195190", x_pos, y_pos2,
                m_pop, sd_pop, cex_base = cex_base)
      
    } else if (input$dist == "rlnorm") {
      x_pos <- max(pop) - 0.1 * x_rng
      plot_dist(pop, distname, "#195190", x_pos, y_pos,
                m_pop, sd_pop, cex_base = cex_base)
      
    } else if (input$dist == "rbeta") {
      x_pos <- max(pop) - 0.25 * x_rng
      plot_dist(pop, distname, "#195190", x_pos, y_pos,
                m_pop, sd_pop, cex_base = cex_base)
    }
  }
  
  # Plot 1a: Population (Tab 1) ----
  output$pop.dist <- renderPlot({ make_pop_plot(cex_base = 1.4) })
  
  # Plot 1b: Population (Tab 3 sidebar, smaller) ----
  output$pop.dist.two <- renderPlot({ make_pop_plot(cex_base = 0.85) })
  
  # Helper: stacked dot plot in base R ----
  stacked_dotplot <- function(x, main = "", cex_dot = 0.8, col_dot = "#195190",
                              n_bins = 15) {
    
    brks  <- seq(min(x), max(x), length.out = n_bins + 1)
    binw  <- brks[2] - brks[1]
    bins  <- findInterval(x, brks, rightmost.closed = TRUE)
    bins  <- pmax(1L, pmin(bins, n_bins))
    
    bin_cx <- brks[bins] + binw / 2
    counts <- table(bins)
    
    max_stack <- max(counts)
    
    plot(NULL,
         xlim = c(min(brks) - binw, max(brks) + binw),
         ylim = c(0, max_stack + 1),
         main = main, xlab = "", ylab = "",
         yaxt = "n", bty  = "l",
         cex.main = 1, cex.axis = 1.6)
    
    rank_in_bin <- ave(seq_along(bins), bins, FUN = seq_along)
    
    symbols(bin_cx, rank_in_bin - 0.5,
            circles  = rep(binw / 2 * 0.9, length(x)),
            inches   = FALSE,
            add      = TRUE,
            bg       = col_dot,
            fg       = col_dot)
  }
  
  # Plot 2: 8 sample dot-plots ----
  output$sample.dist <- renderPlot({
    
    y <- samples()
    
    op <- par(mfrow = c(2, 4), mar = c(3, 1, 3, 1), oma = c(0, 0, 0, 0))
    on.exit(par(op))
    
    for (i in 1:8) {
      col_data <- y[, i]
      m   <- round(mean(col_data), 2)
      s   <- round(sd(col_data),   2)
      
      stacked_dotplot(col_data, main = paste("Sample", i))
      
      mtext(paste0("x\u0304 = ", m, "  SD = ", s),
            side = 3, line = -1.5, cex = 1.0, col = "black")
    }
  })
  
  # Text: number of samples ----
  output$num.samples <- renderText({
    paste0("... continuing to Sample ", input$k, ".")
  })
  
  # Plot 3: Sampling distribution ----
  output$sampling.dist <- renderPlot({
    
    n      <- input$n
    k      <- input$k
    pop    <- parent()
    m_pop  <- round(mean(pop), 4)
    sd_pop <- round(sd(pop),   4)
    
    means   <- sample_means()
    m_samp  <- round(mean(means), 4)
    sd_samp <- round(sd(means),   4)
    
    if (input$dist == "runif" && u_min() == u_max()) return(invisible(NULL))
    
    d      <- density(means)
    x_rng  <- diff(range(means))
    y_pos  <- max(d$y) - 0.1 * max(d$y)
    x_pos  <- if (m_samp > 0) min(means) + 0.1 * x_rng
    else             max(means) - 0.1 * x_rng
    
    hist(means, freq = FALSE,
         col = "#009499", border = "white",
         main = "Sampling Distribution*",
         xlab = "Sample means", ylab = "",
         cex.main = 1.4, cex.axis = 1.2, cex.lab = 1.2)
    
    lines(d, col = "#009499", lwd = 2)
    
    # --- NEW: overlay ideal normal curve using CLT-predicted mean and SE ----
    se_clt <- sd_pop / sqrt(n)
    x_seq <- seq(min(means) - 4*se_clt, max(means) + 4*se_clt, length.out = 300)
    lines(x_seq, dnorm(x_seq, mean = m_samp, sd = se_clt),
          col = "#E07B39", lwd = 2.5, lty = 2)
    
    legend("topright",
           legend = c("Observed density", "Ideal normal (CLT prediction)"),
           col    = c("#009499", "#E07B39"),
           lwd    = c(2, 2.5),
           lty    = c(1, 2),
           bty    = "n",
           cex    = 1.0)
    
    label <- paste0("mean of x\u0304 = ", round(m_samp, 2),
                    "\nSE of x\u0304 = ", round(sd_samp, 4))
    text(x_pos, y_pos, labels = label, col = "black", cex = 1.1)
  })
  
  # Text: sampling distribution caption ----
  output$sampling.descr <- renderText({
    distname <- switch(input$dist,
                       rnorm  = "normal population",
                       rlnorm = "right skewed population",
                       rbeta  = "left skewed population",
                       runif  = "uniform population")
    paste0("*Distribution of means of ", input$k, " random samples, each consisting of ",
           input$n, " observations from a ", distname)
  })
  
  # Text: CLT explanation ----
  output$CLT.descr <- renderText({
    pop   <- parent()
    m_pop <- round(mean(pop), 2)
    s_pop <- round(sd(pop),   2)
    n     <- input$n
    se    <- round(s_pop / sqrt(n), 2)
    
    paste0("According to the Central Limit Theorem (CLT), the distribution of sample means ",
           "(the sampling distribution) should be nearly normal. The mean of the sampling ",
           "distribution should be approximately equal to the population mean (", m_pop, ") ",
           "and the standard error (the standard deviation of sample means) should be ",
           "approximately equal to the SD of the population divided by square root of sample size (",
           s_pop, "/sqrt(", n, ") = ", se, "). Below is our sampling distribution graph. ",
           "To help compare, the population distribution plot is also displayed on the right.")
  })
  
  # --- NEW: Text: normality diagnostic ----
  output$normality.descr <- renderText({
    
    means <- sample_means()
    
    if (input$dist == "runif" && u_min() == u_max())
      return("No data to evaluate.")
    
    # Shapiro-Wilk is limited to 5000 observations; subsample if needed
    sw_data <- if (length(means) > 5000) sample(means, 5000) else means
    sw      <- shapiro.test(sw_data)
    W       <- round(sw$statistic, 4)
    p       <- signif(sw$p.value,  3)
    n       <- input$n
    k       <- input$k
    
    # Plain-language verdict based on W (more stable than p for teaching)
    verdict <- if (W >= 0.98) {
      "This is an excellent match — the sampling distribution is very close to a perfect normal curve."
    } else if (W >= 0.95) {
      "This is a good match — the sampling distribution closely follows a normal curve."
    } else if (W >= 0.90) {
      "This is a moderate match — the sampling distribution is roughly normal, but some deviation is visible."
    } else {
      "This is a poor match — the sampling distribution departs noticeably from normal. Try increasing the sample size (n)."
    }
    
    p_note <- if (p < 0.05) {
      paste0("The p-value (", p, ") flags a statistically significant departure from normality. ",
             "With ", k, " sample means, the test is very sensitive and can detect even tiny, ",
             "unimportant deviations — so focus on the W statistic and the curve overlap above.")
    } else {
      paste0("The p-value (", p, ") gives no statistically significant evidence against normality.")
    }
    
    paste0("Shapiro\u2013Wilk normality test on the ", k, " sample means: ",
           "W = ", W, ".  ", verdict, "  ", p_note)
  })
}

# Launch app -------------------------------------------------------------------
shinyApp(ui = ui, server = server)