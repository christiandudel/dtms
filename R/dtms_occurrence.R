states <- c("A","A","B","B","A","A",NA,"A","A","C","B","B",NA,NA,"C","C","C","D")

spells <- rle(states)$values

occurences <- table(spells,useNA="always")

expanded_occurences <- sapply(occurences,function(x) 1:x)

names_expanded <- names(expanded_occurences)

occurences <- unlist(expanded_occurences)

ordering <- unlist(lapply(names_expanded,function(x) {
  if(is.na(x)) which(is.na(spells)) else which(x==spells)
    }))

result <- numeric(length(occurences))

result[ordering] <- occurences

result <- inverse.rle(list(values=result,lengths=rle(states)$lengths))

