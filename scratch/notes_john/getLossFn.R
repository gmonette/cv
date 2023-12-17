getLossFn <- function(criterion){
  expressions <- as.vector(as.character(functionBody(criterion)))
  which <- grepl("result <- mean\\(", expressions)
  if (!(any(which))) stop(deparse(substitute(criterion)),
                          " does not set 'result <- mean(. . .)'")
  expression <- sub(")$" , "",
                    sub("result <- mean\\(", "",
                        expressions[which][sum(which)])
  )
  eval(parse(text=paste0("function(y, yhat) ", expression)))
}
