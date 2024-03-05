dataFold <- function(data, indices, sign = -1, ...) {
  # data: an indexable data structure
  # indices: an index vector
  # sign: if -1 (default) select elements not in indices,
  #       if 1 select elements in indices
  UseMethod("dataFold")
}

dataFold.default <- function(data, indices, sign = -1, ...) {
  # default method for rectangular data
  data[sign*indices, ]
}

data(Duncan, package="carData")
Duncan$index <- seq_len(nrow(Duncan))

set.seed(123)
(folds <- cv::folds(n=nrow(Duncan), 5))
(indices <- sort(cv::fold(folds, 2)))

head(dataFold(Duncan, indices))
dim(dataFold(Duncan, indices))

dataFold(Duncan, indices, sign=1)

D <- vector(45, mode="list")
names(D) <- rownames(Duncan)
for (i in 1:45){
  D[[i]] <- Duncan[i, ,drop=FALSE]
}
class(D) <- "rowdata"
head(D)

dataFold.rowdata <- function(data, indices, sign = -1, ...) {
  data[sign*indices]
}

indices
dataFold(D, indices, sign=1)
