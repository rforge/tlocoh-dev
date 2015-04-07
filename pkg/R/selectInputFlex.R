# Customizable selectInput control for Shiny apps
#
# An alternative form of the selectInput control for Shiny apps which allows you to specify the style for the select box

selectInputFlex <- function (inputId, label, choices, selected = NULL, multiple = FALSE, selectize = TRUE, width = NULL, size = NULL, selectStyle=NULL) {
    choices <- shiny:::choicesWithNames(choices)
    if (is.null(selected)) {
        if (!multiple) 
            selected <- shiny:::firstChoice(choices)
    } else selected <- shiny:::validateSelected(selected, choices, inputId)
    
    if (!is.null(size) && selectize) {
        stop("'size' argument is incompatible with 'selectize=TRUE'.")
    }
    
    selectTag <- tags$select(id = inputId, class = if (!selectize) "form-control", size = size, shiny:::selectOptions(choices, selected))
    
    if (multiple) selectTag$attribs$multiple <- "multiple"
    res <- div(class = "form-group shiny-input-container", style = if (!is.null(width)) paste0("width: ", validateCssUnit(width)), shiny:::controlLabel(inputId, label), div(selectTag, style=selectStyle))
    if (!selectize) return(res)
    
    shiny:::selectizeIt(inputId, res, NULL, nonempty = !multiple && !("" %in% choices))
    
}