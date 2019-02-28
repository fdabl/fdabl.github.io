library('shiny')
library('mvtnorm')
library('ggplot2')
library('latex2exp')
theme_set(theme_bw())


plot_all <- function(sd1, sd2, rho, x2val, limits = c(-10, 10), densmult = 5) {
  S <- rbind(
    c(sd1^2, sd1*sd2*rho),
    c(sd1*sd2*rho, sd2^2)
  )
  
  x1 <- seq(limits[1], limits[2], length.out = 100)
  x2 <- x1
  grid <- expand.grid(x1 = x1, x2 = x2)
  
  d <- data.frame(grid, prob = dmvnorm(expand.grid(x1, x2), mean = c(0, 0), S))
  
  cond_mean <- rho*sd1/sd2 * x2val
  cond_sd <- sd1 * (1 - rho^2)
  densx1 <- densmult * dnorm(d$x1, cond_mean, cond_sd)
  
  p <- ggplot(d, aes(x = x1, y = x2, z = prob)) +
    geom_contour(aes(colour = '1')) +
    coord_fixed(xlim = limits, ylim = limits - c(-2, 2)) +
    geom_segment(aes(x = cond_mean, xend = cond_mean, y = -10, yend = max(densx1)), linetype = 'dotted', color = 'grey60') +
    # geom_vline(xintercept = cond_mean, linetype = 'dotted') +
    stat_function(fun = function(x) densmult * dnorm(x, 0, sd1), aes(colour = '2')) +
    stat_function(fun = function(x) densmult * dnorm(x, cond_mean, cond_sd), size = 1.2, aes(colour = '3')) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(), 
      axis.line = element_line(colour = 'black'),
      text = element_text(size = 20),
      legend.position = 'right'
    ) +
    scale_colour_manual(
      name = '',
      labels = c('Joint', 'Marginal', 'Conditional'),
      values = c('blue', 'black', 'purple')
    ) +
    ylab(TeX('$$X_2$$')) +
    xlab(TeX('$$X_1$$'))
  
  p
}

ui <- fluidPage(
  # titlePanel("Joint, Marginal & Conditional Gaussians"),
  titlePanel("Welcome to the Gaussians!"),
  sidebarLayout(
    sidebarPanel(
      tags$head(tags$script('$(document).on("shiny:connected", function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            Shiny.onInputChange("innerHeight", window.innerHeight);
                            });
                            $(window).resize(function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            Shiny.onInputChange("innerHeight", window.innerHeight);
                            });
                            ')),
      withMathJax(),
      sliderInput("rho",
                  # "Correlation:",
                  label = "Correlation \\( \\rho \\)",
                  min = -1,
                  max = 1,
                  value = .4, step = .05
      ),
      
      sliderInput(
        "x2",
        label = "Value of \\( X_2 \\)",
        min = -8,
        max = 8,
        value = 2
      ),
      
      sliderInput("sd1",
                  label = "Standard Deviation of \\( X_1 \\)",
                  min = 1,
                  max = 3,
                  value = 2,
                  step = 0.1
      ),
      
      sliderInput("sd2",
                  label = "Standard Deviation of \\( X_2 \\)",
                  min = 1,
                  max = 3,
                  value = 2,
                  step = 0.1
      )
      ),
    mainPanel(
      plotOutput("distPlot")
    )
    )
  )


server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    sd1 <- input$sd1
    sd2 <- input$sd2
    rho <- input$rho
    plot_all(sd1 = sd1, sd2 = sd2, rho = rho, x2val = input$x2, densmult = 10)
    
  },
  
  height = reactive(ifelse(!is.null(input$innerWidth),input$innerHeight*3/5,0)),
  width = reactive(ifelse(!is.null(input$innerWidth),input$innerWidth*3/5,0))
  )
  
}

shinyApp(ui = ui, server = server)