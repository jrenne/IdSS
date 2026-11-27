# Run after installation: Rscript --vanilla inst/validation/check-help-examples.R [output-directory]
args <- commandArgs(trailingOnly = TRUE)
destination <- if (length(args)) args[1] else tempfile("idss-help-examples-")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
destination <- normalizePath(destination)
library(IdSS)
topics <- getNamespaceExports("IdSS")
topics <- sort(topics[vapply(topics, function(name) {
  is.function(getExportedValue("IdSS", name))
}, logical(1))])
libpaths <- file.path(destination, "library-paths.rds")
saveRDS(.libPaths(), libpaths)
runner <- file.path(destination, "run-one-example.R")
writeLines(c(
  "args <- commandArgs(trailingOnly=TRUE)",
  ".libPaths(readRDS(args[2]))",
  "library(IdSS)",
  "set.seed(20261126)",
  'grDevices::pdf(file=tempfile(fileext=".pdf"))',
  'h <- do.call(utils::help, list(topic=args[1],package="IdSS"))',
  'if (!length(h)) stop("Missing help topic: ",args[1])',
  "rd <- utils:::.getHelpFile(h)",
  'if (!any(vapply(rd,function(x)identical(attr(x,"Rd_tag"),paste0(intToUtf8(92),"examples")),logical(1)))) stop("Missing example: ",args[1])',
  'do.call(utils::example,list(topic=args[1],package="IdSS",ask=FALSE,echo=FALSE,run.dontrun=TRUE,run.donttest=TRUE))',
  'cat("EXAMPLE_COMPLETED",intToUtf8(10))',
  "grDevices::dev.off()"
), runner)
results <- lapply(topics, function(topic) {
  log <- file.path(destination, paste0(topic, ".log"))
  started <- proc.time()[["elapsed"]]
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(runner), shQuote(topic), shQuote(libpaths)),
    stdout = log, stderr = log, timeout = 180
  )
  completed <- status == 0L && any(grepl("EXAMPLE_COMPLETED", readLines(log, warn = FALSE), fixed = TRUE))
  cat(if (completed) "PASS" else "FAIL", topic, "\n")
  data.frame(
    function_name = topic, passed = completed, exit_status = status,
    seconds = proc.time()[["elapsed"]] - started, log = log
  )
})
results <- do.call(rbind, results)
utils::write.csv(results, file.path(destination, "results.csv"), row.names = FALSE)
if (!all(results$passed)) stop("Failed help examples; see ", destination)
cat(nrow(results), "public function examples passed in separate fresh R sessions.\n")
