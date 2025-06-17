create_parameter_html <- function(args, output_file = "model_parameters.html") {
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
    "th, tr td:first-child { background-color: #3498db; color: white; }",
    ".table > tbody > tr > td:first-child { background-color: #3498db !important; color: white !important; }",
    "</style>",
    "</head>",
    "<body>",
    "<h1>Model Parameters</h1>"
  )

  create_table_html <- function(mat, title) {
    table_html <- kable(mat, format = "html", escape = FALSE,
                        caption = title) %>%
      kable_styling(bootstrap_options = c("striped", "hover"),
                    full_width = FALSE,
                    position = "left")
    return(table_html)
  }

  # Add lambda tensors
  html_content <- c(html_content,
                    "<h2>Lambda Tensors</h2>")

  for(name in names(args$lam.tensor)) {
    html_content <- c(html_content,
                      paste("<h3>Lambda Tensor", name, "</h3>"),
                      create_table_html(args$lam.tensor[[name]],
                                        paste("Lambda Tensor", name)))
  }

  # Add mu parameters
  html_content <- c(html_content,
                    "<h2>Mu Parameters</h2>",
                    create_table_html(t(args$mu), "Mu Parameters"))

  # Add Q matrix
  html_content <- c(html_content,
                    "<h2>Q Matrix</h2>",
                    create_table_html(args$Q, "Q Matrix"))

  # Close HTML
  html_content <- c(html_content,
                    "</body>",
                    "</html>")

  # Write to file
  writeLines(paste(html_content, collapse = "\n"), output_file)

  # Open in browser
  browseURL(output_file)
}
