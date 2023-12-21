getLossFn <- function(cv){
  eval(parse(text=paste0("function(y, yhat) {",
                         paste(attr(cv, "casewise loss"), collapse="\n"),
                         "}")))
}

mse <- function(y, yhat){
  result <- mean((y - yhat)^2)
  attr(result, "casewise loss") <- c("result <- (y - yhat)^2", "result")
  result
}
