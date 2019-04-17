library(shiny)
library(romero.gateway)
library(EBImage)

ui <- fluidPage(
  
  titlePanel("Image segmentation with OMERO and EBImage"),
  
  verticalLayout(
    h4("OMERO Server:"),
    inputPanel(
      textInput('server', "Hostname", placeholder = "demo.openmicroscopy.org", value = "demo.openmicroscopy.org"),
      textInput('user', "Username", value = ""),
      passwordInput('pass', "Password", value = "")
    ),
    h4("Image:"),
    inputPanel(
      numericInput('image', "Image ID", value = 0, min = 0, step = 1),
      numericInput('c', "Channel", value = 1, min = 1, step = 1),
      numericInput('z', "Z", value = 1, min = 1, step = 1),
      numericInput('t', "T", value = 1, min = 1, step = 1)
    ),
    h4("Segementation algorithm settings:"),
    inputPanel(
      numericInput('w', "w", value = 10, min = 1, step = 1),
      numericInput('h', "h", value = 10, min = 1, step = 1),
      numericInput('offset', "offset", value = 0.1, min = 0.01, step = 0.01),
      numericInput('filter', "filter size", value = 3, min = 1, step = 2),
      actionButton('run_button', "Perform Segmentation", icon("refresh"))
    ),
    
    mainPanel(
      width = 12,
      fluidRow(
        column(6, style='padding:10px;', plotOutput('image')),
        column(6, style='padding:10px;', plotOutput('threshold'))
      ),
      fluidRow(
        column(
          12, 
          style='padding:10px;', 
          div(
            style='overflow-x: scroll',
            dataTableOutput('table')
          )
        )
      )
    ),
    
    inputPanel(
      actionButton('save_button', "Save results back to OMERO", icon("save"))
    )
  )
)

server <- function(input, output, session) {
  observeEvent(
    eventExpr = input$run_button,
    handlerExpr = {
      
      s <- isolate(input$server)
      u <- isolate(input$user)
      p <- isolate(input$pass)
      i <- isolate(input$image)
      z <- isolate(input$z)
      t <- isolate(input$t)
      c <- isolate(input$c)
      w <- isolate(input$w)
      h <- isolate(input$h)
      o <- isolate(input$offset)
      f <- isolate(input$filter)
      
      server <- OMEROServer( host = s, username = u, password = p)
      server <- connect(server)
      
      image <- loadObject(server, "ImageData", i)
      
      pixels <- getPixelValues(image, z, t, c)
      
      disconnect(server)
      
      ebimage <- EBImage::Image(data = pixels, colormode = 'Grayscale')
      ebimage <- normalize(ebimage)
      EBImage::display(ebimage, method = "raster")
      
      threshImg <- thresh(ebimage, w=w, h=h, offset=o)
      threshImg <- medianFilter(threshImg, size=f)
      threshImg <- fillHull(threshImg)          
      threshImg <- bwlabel(threshImg)
      
      shapes = computeFeatures.shape(threshImg)
      moments = computeFeatures.moment(threshImg)
      features <<- merge(shapes, moments, by=0, all=TRUE)
      
      output$image <- renderPlot({
        EBImage::display(ebimage, method = "raster")
      })
      
      output$threshold <- renderPlot({
        EBImage::display(colorLabels(threshImg), method = "raster")
      })
      
      output$table <- renderDataTable(
        options = list(
          pageLength = 10),
        features)
    }
  )
  observeEvent(
    eventExpr = input$save_button,
    handlerExpr = {
      rois <- data.frame(x = c(0), y = c(0), rx = c(0), ry = c(0), w = c(0), h = c(0), theta = c(0), text = c('remove'), stringsAsFactors = FALSE)
      for (index in 1:length(features[,1])) {
        x <- features[index,8]
        y <- features[index,9]
        r1 <- features[index,10]
        ecc <- features[index,11]
        r2 <- sqrt(r1^2 * (1 - ecc^2))
        theta <- features[index,12]
        rois <- rbind(rois, c(x, y, r1, r2, NA, NA, -theta, as.character(index)))
      }
      rois <- rois[-1,]
      rownames(rois) <- c()
      
      s <- isolate(input$server)
      u <- isolate(input$user)
      p <- isolate(input$pass)
      i <- isolate(input$image)
      z <- isolate(input$z)
      t <- isolate(input$t)
      server <- OMEROServer( host = s, username = u, password = p)
      server <- connect(server)
      image <- loadObject(server, "ImageData", i)
      
      addROIs(image, z = z, t = t, coords = rois)
      
      csvFile <- "/tmp/ROI_features.csv"
      write.csv(features, file = csvFile)
      fileannotation <- attachFile(image, csvFile)
      
      disconnect(server)
      
      showModal(modalDialog(
        title = "Saved",
        "ROIs and results table have been attached to the image."
      ))
    }
  )
}

shinyApp(ui = ui, server = server)

