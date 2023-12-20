getLossFn <- function(criterion){
  expressions <- as.vector(as.character(functionBody(criterion)))
  which <- grepl("result <- mean\\(", expressions)
  if (!(any(which))) stop(deparse(substitute(criterion)),
                          " does not set 'result <- mean(. . .)'")
  expression <- sub("\\)[:space:]*$", "",
                    sub("result <- mean\\(", "result <- ",
                        expressions[which][sum(which)])
  )
  expressions[which] <- expression
  which <- grepl("attr\\(", expressions)
  expressions <- expressions[!which]
  eval(parse(text=paste0("function(y, yhat) ",
                         paste(expressions, collapse="\n"),
                         "}")))
}
