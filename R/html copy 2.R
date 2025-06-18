create_parameter_html <- function(args_list, cell_names = NULL, output_file = "model_parameters.html") {
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
    ".epoch-container {",
    "        display: flex;",
    "        justify-content: flex-start;  # Changed from space-between to flex-start",
    "        gap: 5px;                     # Reduced from 10px to 5px",
    "    }",
    ".epoch {",
    "        flex: 1;",
    "        padding: 0 2px;              # Reduced from 5px to 2px",
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
    
    # Try to find name in the cell_names list
    if (!is.null(tensor_name) && !is.null(epoch_cell_names$lam.tensor[[tensor_name]])) {
      return(epoch_cell_names$lam.tensor[[tensor_name]][i,j])
    } else if (!is.null(epoch_cell_names$Q) && all(dim(mat) == dim(epoch_cell_names$Q))) {
      # Check for Q matrix by dimensions and structure
      return(epoch_cell_names$Q[i,j])
    } else if (!is.null(epoch_cell_names$mu) && nrow(mat) == 1) {
      # Check for mu vector by structure (transposed to 1-row matrix)
      return(epoch_cell_names$mu[j])
    }
    return("")
  }

  # Modify create_table_html function to keep original headers and handle cell names correctly
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
      table_html <- kable(tagged$values, format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("hover"),
                     full_width = FALSE,
                     position = "left") %>%
        row_spec(1, background = "white", color = "black")
    } else {
      table_html <- kable(tagged$values, format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("hover"),
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
  browseURL(output_file)
}
