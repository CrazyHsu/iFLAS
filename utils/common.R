parseArg = function(arg, pattern, msg){
    pattern = paste0('^-', pattern, '=')
    if(grepl(pattern, arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
            stop(paste0('Please specify the value of -', msg))
        }else{
            return(arg.split[2])
        }
    }
}

parseArgNum = function(arg, pattern, msg){
    pattern = paste0('^-', pattern, '=')
    if(grepl(pattern, arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
            stop(paste0('Please specify the value of -', msg))
        }else{
            return(as.numeric(arg.split[2]))
        }
    }
}

parseArgAsNum = function(arg, pattern, msg){
    return(parseArgNum(arg, pattern, msg))
}

parseArgStrs = function(arg, pattern, msg){
    pattern = paste0('^-', pattern, '=')
    if(grepl(pattern, arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
            stop(paste0('Please specify the value of -', msg))
        }else{
            return(strsplit(arg.split[2], ',', fixed = T)[[1]])
        }
    }
}

parseArgNums = function(arg, pattern, msg){
    pattern = paste0('^-', pattern, '=')
    if(grepl(pattern, arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
            stop(paste0('Please specify the value of -', msg))
        }else{
            strs = strsplit(arg.split[2], ',', fixed = T)[[1]]
            return(as.numeric(strs))
        }
    }
}

existNotNull = function(arg){
    if(exists(arg) && !is.null(eval(parse(text=arg)))){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

outline <- function(df, txt = TRUE){

  m <- nrow(df);n <- ncol(df)
  plot( c(1,n), c(min(df), max(df)), type = "n", main = "The outline graph of Data", xlab = "Number" , ylab = "Value")
  for (i in 1:m){
    lines(as.numeric(df[i,]), col=i)
    if (txt == TRUE){
      k <- dimnames(df)[[1]][i]
      text(1+(i-1)%%n, df[i,1+(i-1)%%n], k)
    }
  }
}

# m=matrix(c(99,94,93,100,100,
#          99,88,96,99,97,
#          100,98,81,96,100), byrow=TRUE,nrow=3
#        )
# df=as.data.frame(m)
# 
# outline(df,txt=TRUE)