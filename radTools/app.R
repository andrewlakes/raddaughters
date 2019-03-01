#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(plotly)
library(Scale)
library(DiagrammeR)

library(shiny)


unitList <- list("mCi" = 2.22E9, 
                 "uCi" = 2.22E6, 
                 "nCi" = 2.22E3,
                 "pCi" = 2.22,
                 "dpm" = 1,
                 "kBq" = 6E4,
                 "mBq" = 6E7,
                 "gBq" = 6E10,
                 "tBq" = 6E13,
                 "Ci" = 2.22E12,
                 "Bq" = 60,
                 "dps" = 60,
                 "kCi" = 2.22E15
                 )



# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Radiation Calculator"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(10,
              sliderInput("bins",
                         "Number of bins:",
                         min = 1,
                         max = 50,
                         value = 30)
          )
        ),
        fluidRow(
          column(5,
                 numericInput("input",
                              "Value:",
                              value = 1)
                 ),
          column(5,
                 selectInput("units_in",
                             "Units:",
                             choices = unitList, selected = 1)
                 )
        ),
        fluidRow(
          column(5,
                 h5(textOutput('result'))
          ),
          column(5,
                 selectInput("units_out",
                             "Units:",
                             choices = unitList, selected = 1)
          )
        ),
        fluidRow(
          textInput("isotope",
                    'Isotope:',
                    value = '227AC')
        ),
        
        radioButtons("radio", h5("scheme choice"),
                                 choices = c("Isotopes", 
                                             "Nonsense"),
                                 selected = "Isotopes")
        # selectInput("select",
        #             "Scheme choice:", 
        #             choics)
        
        
        
        
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot"),
         
         textOutput("Iso"),
         
         grVizOutput("decayfigure")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
   
   output$result <- renderText({
     formatC(as.numeric(input$units_in) * 
             input$input / 
             as.numeric(input$units_out), format = "e", digits = 4)
   })
   
   output$Iso <- renderText({
     
     decayChain <- readRDS('decayLists/227AC')
     
     decayChain[[1]]$t12
     
     
   })
  
   output$decayfigure <- renderGrViz({
     plotarg <- switch(input$radio,
                       "Isotopes" = "digraph { graph [overlap = true, fontsize = 10]                      
node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50, fillcolor = '#DDFFEB', fontname = Helvetica]
                       1[label='@@1'] 2[label='@@2'] 3[label='@@3'] 4[label='@@4'] 5[label='@@5'] 6[label='@@6'] 7[label='@@7'] 8[label='@@8'] 9[label='@@9']
                       edge[color=black] 1->2 
   } 
                       [1]: '227AC' \n [2]: '223FR' \n [3]: '227TH' \n [4]: '223RA' \n [5]: '223RA' \n [6]: '219RN' \n [7]: '219RN' \n [8]: '215PO' \n [9]: '211PB' ",
                       "Nonsense" = "
     digraph Isotopes {
                       
                       
                       graph [overlap = true, fontsize = 10]
                       
                       
                       node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50,
                       fillcolor = '#DDFFEB', fontname = Helvetica]
                       A [label = '@@1']; B; C; D; E; F
                       
                       node [shape = circle, fixedsize = true, width = 0.9, penwidth = 2.0,]
                       1; 2; 3; 4; 5; 6; 7; 8
                       
                       edge[color=black]
                       A->1 B->2 B->3 B->4 C->A
                       1->D E->A 2->4 1->5 1->F
                       E->6 4->6 5->7 6->7 3->8
                       C->B
   }
                       
                       [1]: 'top'
                       ")
     grViz(plotarg)
   })
   
   
   
   
}

# Run the application 
shinyApp(ui = ui, server = server)

