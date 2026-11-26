test_that("every public function has a parsed example using that function", {
  db <- tools::Rd_db("IdSS")
  exports <- getNamespaceExports("IdSS")
  functions <- exports[vapply(exports, function(name) {
    is.function(getExportedValue("IdSS",name))
  }, logical(1))]
  for (name in functions) {
    rd <- db[[paste0(name,".Rd")]]
    expect_false(is.null(rd), info=name)
    if (is.null(rd)) next
    sections <- Filter(function(x) identical(attr(x,"Rd_tag"),paste0(intToUtf8(92),"examples")),rd)
    expect_true(length(sections)>0,info=name)
    if (!length(sections)) next
    path <- tempfile(fileext=".R")
    tools::Rd2ex(rd,out=path,commentDontrun=FALSE,commentDonttest=FALSE)
    expressions <- parse(path)
    unlink(path)
    expect_true(name %in% all.names(expressions,functions=TRUE),info=name)
  }
})
