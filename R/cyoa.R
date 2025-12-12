library(jsonlite)
json_file <- Sys.getenv("JSON_OPTIONS");
information <- jsonlite::fromJSON(json_file);
options <- information[["options"]]
paths <- information[["paths"]]
inputs <- information[["inputs"]]
outputs <- information[["outputs"]]
