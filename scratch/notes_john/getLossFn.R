getLossFn <- function(cv){
  eval(parse(text=paste0("function(y, yhat) {",
                         attr(cv, "casewise loss"),
                         "}")))
}
