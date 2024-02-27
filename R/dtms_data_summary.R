### Method
summary.dmts_data <- function(data,
                              dtms=NULL,
                              fromvar="from",
                              tovar="to",
                              weightvar=NULL) {

    # Weights per transition
    if(is.null(weightvar)) data$COUNT <- 1 else
      data <- dtms_rename(data,weightvar,"COUNT")

    # Aggregate
    formal <- paste0("COUNT~",fromvar,"+",tovar)
    formal <- as.formula(formal)
    result <- aggregate(formal,data,FUN=sum)

    # Order
    ordering <- order(result[,fromvar],result[,tovar])
    result <- result[ordering,]

    # Proportion
    N <- sum(result$COUNT)
    result$PROP <- result$COUNT/N

    # Raw transition probabilities
    probs <- tapply(result$COUNT,
                    result[,fromvar],
                    FUN=function(x) x/sum(x))
    probs <- unlist(probs)
    result$PROB <- probs

    # If dtms is provided match
    if(is.null(dtms)) return(result) else {
      ### CONTINUE HERE
    }
}



