library("boot")

bootstrap <- function(input){
# Returns mean and median of the input data set along with 95% confidencre intervals.
  mea <- function(x,y) {return(mean(sample(input, 300, replace=TRUE)))}
  meani <- boot(data=input,statistic=mea, R=1000)
  cimean <- boot.ci(meani, conf=0.95, type="norm")
  med <- function(x,y) {return(median(sample(input, 300, replace=TRUE)))}
  mediani <- boot(data=input,statistic=med, R=1000)
  cimed <- boot.ci(mediani, conf=0.95, type="norm")
  return(list("Mean:", mean(input), cimean$norm, "Median:", median(input), cimed$norm))
}

