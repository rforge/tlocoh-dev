#' Interactive selection of hullsets from a locoh-hullset
#'
#' Visually select one hullset per individual from a \link[tlocoh]{LoCoH-hullset} 
#'
#' @param lhs A LoCoH-hullset object
#' @param selection An object of class \emph{locoh.selection} containing a named list (one per individual) of selected hullset parameter values
#' @param gmap The name of a background image that will be downloaded from Google: \code{"none"}, 
#' \code{"roadmap"}, \code{"satellite"}, \code{"hybrid"}, or \code{"terrain"}. May also be a object of type \code{locoh.gmap}, see Notes.
#'
#' @return An object of class \emph{locoh.selection} containing a list of selected hullset parameter value for each id (see Details)
#'
# '@details This function will open up a Shiny app which allows you to visually inspect hullsets together
#' with plots of isopleth area and isopleth edge:area ratio, and select the value of 'k', 'a', or 'r' that does the
#' best job balancing over- and under-estimation.
#'
#' This function can be useful when you have a locoh-hullset object with multiple hullsets for each individual over a range of \emph{k} or \emph{a} values, and
#' you want to pick one hullset per individual for the rest of the analysis. Note the range of parameter values should be uniform. 
#'
#' This function requires using RStudio and the Shiny package. To return the selection, be sure to click the 'Save and Return' button in the shiny app. The object returned can be passed to the function again as the 'selection' argument.
#'
#' To display an image from Google in the background, set gmap to \code{"roadmap"}, \code{"satellite"}, \code{"hybrid"}, or \code{"terrain"}. 
#' This requires an internet connection. You may also set gmap to an object of type \code{locoh.gmap}, so the image(s) don't have to be 
#' downloaded each time. See \code{lhs.gmap} (in the tlocoh.dev package). 
#'
#' @seealso code{\link{lhs.selection}}, code{\link{lhs.gmap}}
#'
#' @export

lhs.shiny.select <- function(lhs, selection=NULL, gmap="none") {

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    if (!requireNamespace("dismo", quietly=TRUE)) stop("package dismo required to display a background image, please install")
    if (!requireNamespace("rgdal", quietly=TRUE)) stop("package rgdal required to display a background image, please install")
    if (!requireNamespace("raster", quietly=TRUE)) stop("package raster required to display a background image, please install")            

    ## Get the mode for this lhs
    lhs.mode <- unique(sapply(lhs, function(x) x$mode))
    if (length(lhs.mode) > 1) stop("More than one mode in this hullset")
    
    ## Get the ids for this lhs
    lhs.id <- unique(sapply(lhs, function(x) x$id))
    
    ## Prepare gmapData, which is a list of locoh.gmap objects (mostly blank)
    gmap.types <- c("roadmap", "satellite", "hybrid", "terrain", "none")
    gmapData <- lapply(gmap.types, function(x) list())
    #gmapData <- as.list(rep(NA, length(gmap.types)))
    names(gmapData) <- gmap.types
    for (i in 1:length(gmap.types)) class(gmapData[[i]]) <- c("locoh.gmap", "list")
    gmapData[["none"]] <- "none"

    if (inherits(gmap, "locoh.gmap")) {
        gmapType <- gmap[[1]]$type
        gmapData[[gmapType]] <- gmap
    } else if (is.character(gmap)) {
        if (! gmap %in% gmap.types) stop("Unknown value for gmap")
        gmapType <- gmap
    } else {
        gmapType <- "none"
    }
    
    ## Get the parameter k/a/r values
    lhs.akr_vals <- sort(unique(sapply(lhs, function(x) x[[x$mode]])))
    if (length(lhs.akr_vals)==1) stop(paste("This hullset only has a single value for ", lhs.mode, sep=""))
    
    ## Calculate the interval(s) 
    lhs.akr_min <- min(lhs.akr_vals)
    lhs.akr_max <- max(lhs.akr_vals)
    lhs.akr_vals_step <- unique(lhs.akr_vals[2:length(lhs.akr_vals)] - lhs.akr_vals[1:(length(lhs.akr_vals)-1)])
    if (length(lhs.akr_vals_step) > 1) {
        print(lhs.akr_vals_step)
        stop (paste("The increment in", lhs.mode, "is not uniform", sep=" "))
    }
    lhs.akr_step <- lhs.akr_vals_step
    
    lhs.akr_vals_lst <- as.list(as.character(lhs.akr_vals))
    names(lhs.akr_vals_lst) <- as.character(lhs.akr_vals)
    lhs.akr_vals_lst <- c(list("NA"="NA"), lhs.akr_vals_lst)

    if (is.null(selection)) {
        lhs.akr_choice <- as.list(rep("NA", length(lhs.id)))
        names(lhs.akr_choice) <- lhs.id
        attr(lhs.akr_choice, "mode") <- lhs.mode
        class(lhs.akr_choice) <- c("locoh.selection", "list")
    } else {
        if (!is.list(selection)) stop("'selection' must be a list")
        if (length(selection) != length(lhs.id)) stop("incorrect length for 'selection'")
        if (FALSE %in% (names(selection) %in% lhs.id)) stop("element names in selection must match the ids")
        lhs.akr_choice <- selection
    }

    # Grab all isopleth levels
    iso.levels <- sort(unique(as.numeric(sapply(lhs, function(x) sapply(x$isos, function(y) y$polys$iso.level)))))
    iso.levels.lst <- as.list(iso.levels)
    names(iso.levels.lst) <- iso.levels
    
    ## Initialize a status message
    #display_msg <- "---"
    #rvals <- reactiveValues(msg = "---")
        
    shiny_app <- shinyApp(
    
        ui = fluidPage(
 
            # Application title
            title="Select Hullset Parameters",

            fluidRow(
                column(7,
                    radioButtons(inputId="id", label="ID:", choice=lhs.id, selected=lhs.id[1], inline=TRUE)
                ),

                column(3,
                    sliderInput("akr_val",
                            paste(lhs.mode, " value:", sep=""),
                            min=lhs.akr_min,
                            max=lhs.akr_max,
                            step=lhs.akr_step,
                            animate=FALSE,
                            value = lhs.akr_min)
                ),

                column(2,
                    div(actionLink("cmdSaveAKR", icon=icon("arrow-right fa-lg"), label=""), style="display:inline-block; width:10%; vertical-align:middle; position:relative; left:-10px;"),
                    div(selectInput("akr_choice", label="selection", choices=lhs.akr_vals_lst, selected=lhs.akr_choice[[ lhs.id[1] ]]), style="display:inline-block; width:85%; vertical-align:middle;")
                )

            ),

            fluidRow(
              column(8,
                    plotOutput("isoPlot"),
                    br(),

                    fluidRow(
                        column(4,
                            p(strong("Map options")),
                            checkboxInput("chkShowPoints", label="Show points", value=TRUE),
                            selectInputFlex("lstGmap", "Gmap: ", c("none", "roadmap", "satellite", "hybrid", "terrain"), selected = gmapType, selectStyle="display:inline-block; width:110px; vertical-align:middle;"),
                            textOutput("txtMsg"),
                            tags$style(type="text/css", '#txtMsg {display:block; position:relative; width: 80px; color:#F00; top:-20px}')
                        ),
                        
                        column(8,
                            checkboxGroupInput("chkIsoLevels", label="Isopleth levels", choices = iso.levels.lst, selected=iso.levels, inline=TRUE),                            
                            tags$style(type="text/css", '#sngOpacity {width: 80px; display:inline;}'),
                            div(numericInput("sngOpacity", "Opacity (0..1):", value=0.8, min=0, max=1), style="display:block-inline;")
                            
                        )
                    )

                ),

                column(4,
                    plotOutput("isoareaPlot", height="250px"),
                    br(),
                    plotOutput("isoearPlot", height="250px"),
                    div(actionButton("cmdClose", label="Save and Close"), style="text-align:right; padding-top:1em; padding-bottom:1em;"),
                    tags$style(type="text/css", "#cmdClose {background-color:#E8E8E8;}")
                    
                )
            )

        ),

      server = function(input, output, session) {

          # Create a Progress object
          # http://shiny.rstudio.com/articles/progress.html
          #progress <- shiny::Progress$new()
          # Make sure it closes when we exit this reactive, even if there's an error
          #on.exit(progress$close())
          #progress$set(message = "Downloading", value = 0)

          #output$txtMsg <- renderText({rvals$msg})
          
          output$isoearPlot <- renderPlot({
              lhs.plot.isoear(lhs, id=input$id, title="Isopleth Edge:Area Ratio")
              abline(v=input$akr_val)
          })

          output$isoareaPlot <- renderPlot({
              lhs.plot.isoarea(lhs, id=input$id, title="Isopleth Area")
              abline(v=input$akr_val)
          })

          output$isoPlot <- renderPlot({
              ## Create the expression for the plot command
              #exp.str <- paste("plot(lhs, id=input$id, ", lhs.mode, "=input$akr_val, allpts=input$chkShowPoints, gmap=gmap, cex.allpts=0.3, col.allpts=\"grey50\", iso=T, col.iso.opacity=input$sngOpacity, iso.level=as.numeric(input$chkIsoLevels))", sep="")
              exp.str <- paste("plot(lhs, id=input$id, ", lhs.mode, "=input$akr_val, allpts=input$chkShowPoints, gmap=gmapData[[input$lstGmap]], cex.allpts=0.3, col.allpts=\"grey50\", iso=T, col.iso.opacity=input$sngOpacity, iso.level=as.numeric(input$chkIsoLevels))", sep="")
              expr <- parse(text=exp.str)   ## converts the string into an expression object
              eval(expr)
          })

          observe({
              ## Kicks in when id is changed (only)
              ## Changes the akr_choice control and akr_slider control to the saved version for the new id
              input$id
              updateSelectInput(session, "akr_choice", selected=lhs.akr_choice[[input$id]])
              if (input$akr_choice != "NA") {
                  if (as.character(isolate(input$akr_val)) != input$akr_choice) {
                      updateSliderInput(session, "akr_val", value=as.numeric(input$akr_choice))
                  }
              }
          })
          
          observe({
              ## Kicks in when lstGmap changes (only)
              ## Download a background map if needed
              if (!identical(input$lstGmap, "none")) {
                  if (is.null(gmapData[[input$lstGmap]][[input$id]])) {
                      
                      #display_msg <- "pls wait"
                      #output$txtMsg <- renderText({display_msg})   ## DOESN'T UPDATE UNTIL THE END OF THE FUNCTION, NO GOOD
                      
                      ##
                      ## rvals$msg <<- paste("pls wait", input$id) THIS WORKS BUT THE SCREEN ISN'T UPDATED UNTIL *AFTER* THE IMAGE IS DOWNLOADED
                      
                      #Sys.sleep(0.1)
                      #progress$inc(0.3, detail = "Downloading")  ## PROGRESS BAR APPEARS AT THE TOP OF THE SCREEN
                      
                      #invalidateLater(1, session)
                      
                      # gmapData[[input$lstGmap]] <<- lhs.gmap(lhs, gmap=input$lstGmap, status=TRUE) THIS GRABS MAPS FOR ALL THE IDS. DECIDED NOT TO DO 
                      # THIS. IF THE USER WANTS THAT PERFORMANCE, CAN PASS A LOCOH.GMAP OBJECT TO THE FUNCTION
                      
                      gmapData[[input$lstGmap]][[input$id]] <<- lhs.gmap(lhs, id=input$id, gmap=input$lstGmap, status=TRUE)[[ input$id  ]] 
                      
                      #progress$inc(0.3, detail = "Downloading")
                      # display_msg <- "thanks"
                      #rvals$msg <<- paste("thanks ", input$id, sep="")
                      #output$txtMsg <- renderText({display_msg})
                      #output$txtMsg <- renderText({""})
                  }
              }
          })
          
          observe({
              ## Kicks in when input$akr_choice changes (only)
              input$akr_choice
              lhs.akr_choice[[isolate(input$id)]] <<- input$akr_choice
          })
          
          observe({
              ## Stops the app when the user clicks the 'Save and Close' button
              if (input$cmdClose > 0) stopApp(lhs.akr_choice)
          })
          
          observe({
              ## Saves the current value from the slider to akr_choice
              if (input$cmdSaveAKR > 1) updateSelectInput(session, "akr_choice", selected=as.character(isolate(input$akr_val)))
          })


      }

    )
    
    ra <- runApp(shiny_app)
    
    return(ra)
    
}

