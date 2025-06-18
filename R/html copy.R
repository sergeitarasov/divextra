create_parameter_html <- function(args, cell_names = NULL, output_file = "model_parameters.html") {
  library(knitr)
  library(kableExtra)

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
    "        padding: 4px 8px 16px 8px;  /* Increased bottom padding */",
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
    "        background-color: white;  /* White background */",
    "        padding: 0 2px;          /* Small horizontal padding */",
    "        border-radius: 2px;      /* Rounded corners */",
    "        z-index: 1;              /* Ensure coordinates show above background */",
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
    "<h1>Model Parameters</h1>",
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
    if ("Q" %in% names(args) && identical(mat, args$Q)) {
      is_q <- TRUE
    }
    
    for(i in 1:n) {
      for(j in 1:m) {
        if(mat[i,j] != "0") {
          tags[i,j] <- "bold"  # Tag non-zero cells
        } else if(i > j && !is_q) {
          tags[i,j] <- "lower-zero"  # Tag lower triangle zeros
        }
      }
    }
    
    list(values = mat, tags = tags)
  }

  # Function to get cell name from the names list
  get_cell_name <- function(mat, i, j, tensor_name = NULL) {
    if (is.null(cell_names)) return("")
    
    cell_value <- mat[i,j]
    if (cell_value == "0") return("")
    
    # Try to find name in the cell_names list
    if (!is.null(tensor_name) && !is.null(cell_names$lam.tensor[[tensor_name]])) {
      return(cell_names$lam.tensor[[tensor_name]][i,j])
    } else if (identical(mat, args$Q) && !is.null(cell_names$Q)) {
      return(cell_names$Q[i,j])
    } else if (identical(mat, t(args$mu)) && !is.null(cell_names$mu)) {
      return(cell_names$mu[j])  # For mu vector use j as index since matrix is transposed
    }
    return("")
  }

  # Update create_table_html function
  create_table_html <- function(mat, title = NULL, class = NULL) {
    tagged <- tag_matrix_cells(mat)
    tensor_name <- title  # Pass tensor name for lambda tensors
    
    for(i in 1:nrow(mat)) {
      for(j in 1:ncol(mat)) {
        # Get cell name
        cell_name <- get_cell_name(mat, i, j, tensor_name)
        
        # Create wrapper with cell name
        wrapper <- sprintf("%s<span class='coordinate'>%s</span>",
                         mat[i,j], cell_name)
        
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
    
    # Apply different styling for mu table
    if (!is.null(class) && class == "mu-table") {
      table_html <- kable(tagged$values, format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("hover"),  # Removed 'striped'
                      full_width = FALSE,
                      position = "left") %>%
        row_spec(1, background = "white", color = "black")  # Reset first row styling
    } else {
      table_html <- kable(tagged$values, format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("hover"),  # Removed 'striped'
                      full_width = FALSE,
                      position = "left")
    }
    
    return(table_html)
  }

  # Add lambda tensors
  html_content <- c(html_content,
                    "<h2>Lambda Tensors</h2>")

  for(name in names(args$lam.tensor)) {
    html_content <- c(html_content,
                      paste("<h3>Lambda Tensor", name, "</h3>"),
                      create_table_html(args$lam.tensor[[name]], title = name))
  }

  # Add mu parameters
  html_content <- c(html_content,
                    "<h2>Mu Parameters</h2>",
                    create_table_html(t(args$mu), class = "mu-table"))

  # Add Q matrix
  html_content <- c(html_content,
                    "<h2>Q Matrix</h2>",
                    create_table_html(args$Q))

  # Close HTML
  html_content <- c(html_content,
                    "</body>",
                    "</html>")

  # Write to file
  writeLines(paste(html_content, collapse = "\n"), output_file)

  # Open in browser
  browseURL(output_file)
}
