#' Create HTML Visualization of ClaSSE/GeoSSE Model Parameters
#'
#' @description
#' Creates an interactive HTML visualization of ClaSSE/GeoSSE model parameters, including
#' lambda tensors, mu parameters, and Q matrices for all epochs. The visualization includes
#' toggle-able parameter names and color-coded matrices.
#'
#' @param par.categories.td List containing model parameters and configuration:
#'   \itemize{
#'     \item pars - Parameter values
#'     \item Nstates - Number of states in the model
#'     \item n.epoch - Number of epochs
#'     \item epoch.times - Vector of epoch transition times
#'     \item states - Vector of region/state names
#'   }
#' @param model_name Character string specifying the model name (optional)
#' @param output_file Character string specifying the output HTML file path
#'
#' @return None (invisible). Creates an HTML file at the specified location.
#'
#' @details
#' The function generates an HTML page with:
#' - Model information summary
#' - Lambda tensor matrices for each epoch
#' - Mu parameter tables for each epoch
#' - Q matrices for each epoch
#' - Toggle-able parameter names
#' - Color coding for special matrix elements
#'
#' @examples
#' \dontrun{
#' # Create HTML visualization
#' create_parameter_html(
#'   par.categories.td = model_params,
#'   model_name = "My GeoSSE Model",
#'   output_file = "model_visualization.html"
#' )
#' }
#'
#' @export
create_parameter_html <- function(par.categories.td, model_name = NULL, output_file) {
  # library(knitr)
  # library(kableExtra)

    # Check if output_file is provided
  if (missing(output_file) || is.null(output_file)) {
    stop("output_file must be provided")
  }
  
  # Add .html extension if not present
  if (!grepl("\\.html$", output_file)) {
    output_file <- paste0(output_file, ".html")
  }
  
  Npars <- length(par.categories.td$pars)
  Nstates <- par.categories.td$Nstates
  n.epoch <- par.categories.td$n.epoch
  epoch.times <- par.categories.td$epoch.times
  regions <-par.categories.td$states
  

  #--- preposecess
  args_list <- pars_yaml_to_arrays_td(par.categories.td)

  # genertae parameter names
  cell_names <- vector("list", length = n.epoch)
  for (epoch in seq_along(1:n.epoch)){
    pars <- diversitree:::default.argnames.classe(Nstates)
    pars <- paste0(pars, '.', epoch)
    cell_names.i <- pars_to_arrays(pars, Nstates, regions)
    cell_names[[epoch]] <- cell_names.i
  }


  #--- preposecess

  html_content <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "<title>Model Parameters</title>",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 20px; }",
    "h1 { color: #2c3e50; }",
    "h2 { color: #34495e; margin-top: 30px; }",
    "table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }",
    "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
    "tr:nth-child(even) td:not(:first-child) { background-color: #f2f2f2; }",
    "td:first-child { background-color: #3498db; color: white; }",
    "table.mu-table td:first-child { background-color: inherit; color: inherit; }",
    "table.mu-table tr:nth-child(even) td { background-color: #f2f2f2; }",
    "th { background-color: #3498db; color: white; }",
    "td {",
    "        position: relative;",
    "        padding: 4px 8px 16px 8px;",
    "        line-height: 1.2;",
    "        vertical-align: top;",
    "    }",
    ".cell-wrapper {",
    "        min-height: 20px;",
    "        display: inline-block;",
    "        position: relative;",
    "    }",
    ".coordinate {",
    "        color: #888;",
    "        font-size: 0.8em;",
    "        display: none;",
    "        position: absolute;",
    "        bottom: -12px;",
    "        left: 8px;",
    "        background-color: white;",
    "        padding: 0 2px;",
    "        border-radius: 2px;",
    "        z-index: 1;",
    "    }",
    ".epoch-container {",
    "        display: flex;",
    "        justify-content: space-between;",
    "        gap: 10px;",
    "    }",
    ".epoch {",
    "        flex: 1;",
    "        padding: 0 5px;",
    "    }",
    "    .model-info {",
    "        background-color: #f8f9fa;",
    "        padding: 15px;",
    "        border-radius: 5px;",
    "        margin: 20px 0;",
    "        border-left: 4px solid #3498db;",
    "    }",
    "    .model-info p {",
    "        margin: 5px 0;",
    "        font-family: monospace;",
    "    }",
    "</style>",
    "<script>",
    "function toggleCoordinates() {",
    "  var coords = document.getElementsByClassName('coordinate');",
    "  var display = document.getElementById('showCoords').checked ? 'block' : 'none';",
    "  for(var i = 0; i < coords.length; i++) {",
    "    coords[i].style.display = display;",
    "  }",
    "}",
    "</script>",
    "</head>",
    "<body>",
    "<h1>ClaSSE/GeoSSE model (epoch 1 is the most recent)</h1>",
    "<div class='model-info'>",
    sprintf("        <h3>%s</h3>", if(!is.null(model_name)) model_name else ""),
    sprintf("        <p>Parameter count: %d</p>", Npars),
    sprintf("        <p>Number of states: %d</p>", Nstates),
    sprintf("        <p>Number of epochs: %d</p>", n.epoch),
    sprintf("        <p>Epoch times: %s</p>", paste(epoch.times, collapse=", ")),
    sprintf("        <p>Regions: %s</p>", paste(regions, collapse=", ")),
    "        </div>",
    "<input type='checkbox' id='showCoords' onclick='toggleCoordinates()'> Show parameter names",
    "<br><br>"
  )

  # Mark lower triangle with orange for lambda tensors
  mark_lower_triangle <- function(mat, is_lambda = TRUE) {
    if (!is_lambda) return(mat)
    
    n <- nrow(mat)
    for(i in 1:n) {
      for(j in 1:(i-1)) {  # Below diagonal only
        if(mat[i,j] == "0") {
          mat[i,j] <- cell_spec("0", background = "#FFA726", format = "html")
        }
      }
    }
    return(mat)
  }

  # Function to tag cells based on their position and content
  tag_matrix_cells <- function(mat) {
    n <- nrow(mat)
    m <- ncol(mat)
    tags <- matrix("", n, m)
    
    # Check if it's Q matrix
    is_q <- FALSE
    if ("Q" %in% names(args_list[[1]]) && identical(mat, args_list[[1]]$Q)) {
      is_q <- TRUE
    }
    
    for(i in 1:n) {
      for(j in 1:m) {
        if(mat[i,j] != "0") {
          # Tag all non-zero cells as bold, including Q matrix
          tags[i,j] <- "bold"
        } else if(i > j && !is_q) {
          # Only highlight lower triangle zeros for non-Q matrices
          tags[i,j] <- "lower-zero"
        }
      }
    }
    
    list(values = mat, tags = tags)
  }

  # Function to get cell name from the names list
  get_cell_name <- function(mat, i, j, tensor_name = NULL, epoch_cell_names = NULL) {
    if (is.null(epoch_cell_names)) return("")
    
    cell_value <- mat[i,j]
    if (cell_value == "0") return("")
    
    if (!is.null(tensor_name) && !is.null(epoch_cell_names$lam.tensor[[tensor_name]])) {
        return(epoch_cell_names$lam.tensor[[tensor_name]][i,j])
    } else if (!is.null(epoch_cell_names$Q) && all(dim(mat) == dim(epoch_cell_names$Q))) {
        return(epoch_cell_names$Q[i,j])
    } else if (!is.null(epoch_cell_names$mu) && nrow(mat) == 1) {
        return(epoch_cell_names$mu[j])
    }
    return("")
  }

  # Replace cell_spec calls
  cell_spec <- function(...) {
      kableExtra::cell_spec(...)
  }

  # Update create_table_html function with explicit package references
  create_table_html <- function(mat_list, title = NULL, class = NULL) {
    n_epochs <- length(mat_list)
    
    # Initialize combined matrix with original column headers
    combined_mat <- NULL
    for(epoch in 1:n_epochs) {
      mat <- mat_list[[epoch]]
      if(is.null(combined_mat)) {
        combined_mat <- mat
      } else {
        # Add epoch columns keeping original headers
        combined_mat <- cbind(combined_mat, mat)
      }
    }
    
    # Tag cells and apply styling
    tagged <- tag_matrix_cells(combined_mat)
    
    for(i in 1:nrow(combined_mat)) {
      for(j in 1:ncol(combined_mat)) {
        # Get cell name for the corresponding epoch
        epoch <- ceiling(j/ncol(mat_list[[1]]))
        orig_col <- ((j-1) %% ncol(mat_list[[1]])) + 1
        
        # Get cell name using the correct epoch's cell_names
        if (!is.null(cell_names) && length(cell_names) >= epoch) {
          cell_name <- get_cell_name(mat_list[[epoch]], i, orig_col, title, cell_names[[epoch]])
        } else {
          cell_name <- get_cell_name(mat_list[[epoch]], i, orig_col, title)
        }
        
        # Create wrapper with cell name
        wrapper <- sprintf("%s<span class='coordinate'>%s</span>",
                         combined_mat[i,j], cell_name)
        
        # Update cell styling
        if(tagged$tags[i,j] == "lower-zero") {
          tagged$values[i,j] <- cell_spec(
            wrapper,
            format = "html",
            background = "#343a40",
            color = "white",
            escape = FALSE
          )
        } else if(tagged$tags[i,j] == "bold") {
          tagged$values[i,j] <- cell_spec(
            wrapper,
            format = "html",
            bold = TRUE,
            color = "#dc3545",
            escape = FALSE
          )
        } else {
          tagged$values[i,j] <- cell_spec(
            wrapper,
            format = "html",
            escape = FALSE
          )
        }
      }
    }
    
    # Create table with all epochs
    if (!is.null(class) && class == "mu-table") {
      table_html <- knitr::kable(tagged$values, format = "html", escape = FALSE) %>%
        kableExtra::kable_styling(bootstrap_options = c("hover"),
                                    full_width = FALSE,
                                    position = "left") %>%
        kableExtra::row_spec(1, background = "white", color = "black")
    } else {
      table_html <- knitr::kable(tagged$values, format = "html", escape = FALSE) %>%
        kableExtra::kable_styling(bootstrap_options = c("hover"),
                                    full_width = FALSE,
                                    position = "left")
    }
    
    return(table_html)
  }
  
  # Update main function to handle multiple epochs
  html_content <- c(html_content, 
                   "<div class='epoch-container'>")
  
  # Create content for each epoch
  for(epoch in 1:length(args_list)) {
    html_content <- c(html_content,
                     sprintf("<div class='epoch'><h2>Epoch %d</h2>", epoch))
    
    # Lambda tensors
    html_content <- c(html_content, "<h3>Lambda Tensors</h3>")
    for(name in names(args_list[[epoch]]$lam.tensor)) {
      html_content <- c(html_content,
                      paste("<h4>Lambda Tensor", name, "</h4>"),
                      create_table_html(list(args_list[[epoch]]$lam.tensor[[name]]), 
                                     title = name))
    }
    
    # Mu parameters
    html_content <- c(html_content,
                     "<h3>Mu Parameters</h3>",
                     create_table_html(list(t(args_list[[epoch]]$mu)), 
                                    class = "mu-table"))
    
    # Q matrix
    html_content <- c(html_content,
                     "<h3>Q Matrix</h3>",
                     create_table_html(list(args_list[[epoch]]$Q)))
    
    html_content <- c(html_content, "</div>")
  }
  
  html_content <- c(html_content, "</div>")
  
  # Close HTML
  html_content <- c(html_content,
                    "</body>",
                    "</html>")

  # Write to file
  writeLines(paste(html_content, collapse = "\n"), output_file)

  # Open in browser
  #browseURL(output_file)
}
