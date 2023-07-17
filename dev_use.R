#Notes


library(devtools)

#Create a Temporary File
file <-  tempfile()
dir <- tempdir()

#Get something from install example
base <- system.file('sql', package='piana')
sqls <- dir(base, "*sql", f=TRUE)



#UES this
use_r()
use_test()
use_build_ignore("NAMESPACE_OLD")
usethis::use_package("methods")

#Use pipes
use_pipe(export = TRUE)


create_doc <- function(directory = NULL,
                       database_name_name) {
  if (is.null(directory)) directory <- tclvalue(tkchooseDirectory(
    title = "Choose Folder for Input and Output"))
  rmarkdown::render("R/doc_generator.Rmd", output_dir = directory)
}

#Load only exportable functions
load_all(export_all = FALSE)


devtools::check()
