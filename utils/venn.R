#!/usr/bin/env Rscript
args <- commandArgs()
scriptPath = strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName = basename(scriptPath)
scriptDir = dirname(scriptPath)
args = args[-(1:5)]
source(paste0(scriptDir, '/common.R'))

usage = function(){
    cat(paste0("Usage: ", scriptName) )
    cat(" -p=outputName.pdf input1.lst [inpu2.lst ...]
Note: 1. only support FIVE input
      2. the numbers in result isn't accurate if input is passed in by /dev/* rather than file
Option:
    -p|pdf      FILE    The output figure in pdf[venn.pdf]
    -w|width    INT     The figure width
    -m|main     STR     The main title
    -mainS      DOU     The size of main title[2.5]

    -c|color    STR     The color of each circle's circumference([rainbow], transparent, ...)
    -f|fill     STR     The color of each circle's area[rainbow]
    -lS|lSize   DOU     The size for each area label[2]
    -cN|cNames  STRs    Specify category names
    -cC|cColor  STR     The color for each category name[rainbow]
    -cS|cSize   DOU     The size for each category name[2]
    -cP|cPos    DOU     The position (in degrees) of each category name along the circle
    -lwd        DOU     The width of each circle's circumference[1]
    -lty        INT/STR The dash pattern of each circle's circumference (blank, ...)
    -rD|rDegree DOU     Number of degrees to rotate the entire diagram
    
    -i|inverted         Flip the two-set Venn diagram along its vertical axis
    -margin     DOU     The amount of whitespace around the diagram in grid units [0.1]
    -r|redun            Do not force use only unique elements in each item of the input list
    -h|help             Show help
")
    q(save = 'no')
}

if(length(args) == 0) usage()

myPdf = 'venn.pdf'
mainS = 2.5
color = 'rainbow'
fill = 'rainbow'
lSize = 2
cColor = 'rainbow'
cSize = 2
lwd = 1
margin = 0.1

for(i in 1:length(args)){
    arg = args[i]
    if(arg == '-h' || arg == '-help') usage()
    tmp = parseArg(arg, 'p(df)?', 'p')
    if(!is.null(tmp)){
        myPdf = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'w(idth)?', 'w')
    if(!is.null(tmp)){
        width = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'm(ain)?', 'm')
    if(!is.null(tmp)){
        main = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'mainS', 'mainS')
    if(!is.null(tmp)){
        mainS = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'c(olor)?', 'c')
    if(!is.null(tmp)){
        color = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'f(ill)?', 'f')
    if(!is.null(tmp)){
        fill = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'lS(ize)?', 'lS')
    if(!is.null(tmp)){
        lSize = tmp
        args[i] = NA
        next
    }
    tmp = parseArgStrs(arg, 'cN(ames)?', 'cN')
    if(!is.null(tmp)){
        cNames = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'cC(olor)?', 'cC')
    if(!is.null(tmp)){
        cColor = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'cS(ize)?', 'cS')
    if(!is.null(tmp)){
        cSize = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNums(arg, 'cP(os)?', 'cP')
    if(!is.null(tmp)){
        cPos = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'lwd', 'lwd')
    if(!is.null(tmp)){
        lwd = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'lty', 'lty')
    if(!is.null(tmp)){
        lty = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'rD(egree)?', 'rD')
    if(!is.null(tmp)){
        rDegree = tmp
        args[i] = NA
        next
    }
    tmp = parseArgNum(arg, 'margin', 'margin')
    if(!is.null(tmp)){
        margin = tmp
        args[i] = NA
        next
    }
    if(arg == '-i' || arg == '-inverted'){
        inverted = TRUE
        args[i] = NA
        next
    }
    if(arg == '-r' || arg == '-redun'){
        redun = TRUE
        args[i] = NA
        next
    }
}

args = args[!is.na(args)]
argLen = length(args)

cat(paste0('[DEBUG] ', Sys.time(), ' Check if the following variables are correct as expected:\n'))
cat('\npdf\t'); cat(myPdf)
cat('\nwidth\t'); if(exists('width')) cat(width)
cat('\nmain\t'); if(exists('main')) cat(main)
cat('\nmainS\t'); cat(mainS)
cat('\ncolor\t'); cat(color)
cat('\nfill\t'); cat(fill)
cat('\nlSize\t'); cat(lSize)
cat('\ncNames\t'); if(exists('cNames')) cat(cNames)
cat('\ncColor\t'); cat(cColor)
cat('\ncSize\t'); cat(cSize)
cat('\ncPos\t'); if(exists('cPos')) cat(cPos)
cat('\nlwd\t'); cat(lwd)
cat('\nlty\t'); if(exists('lty')) cat(lty)
cat('\nrDegree\t'); if(exists('rDegree')) cat(rDegree)
cat('\nmargin\t'); cat(margin)
cat('\ninverted\t'); if(exists('inverted')) cat(inverted)
cat('\nredun\t'); if(exists('redun')) cat(redun)
cat('\ninput files\t'); cat(args)
cat('\n')


library(tools)
library(grid)
suppressPackageStartupMessages(library(VennDiagram))

if(exists('width')){
    pdf(myPdf, width = width)
}else{
    pdf(myPdf)
}

data = read.delim(args[1], header = F)
myList = list(data$V1)
if(argLen >= 2){
    for(i in 2:argLen){
        data = read.delim(args[i], header = F)
        myList = c(myList, list(data$V1))
    }
}

if(argLen == 4) color = 'transparent'

names(myList) = basename(file_path_sans_ext(args))

myCmd = 'T = venn.diagram(myList, filename = NULL, cex = lSize, cat.cex = cSize, lwd = lwd'
if(exists('main')) myCmd = paste0(myCmd, ', main = main')
if(exists('mainS')) myCmd = paste0(myCmd, ', main.cex = mainS')

if(color == 'rainbow'){
    myCmd = paste0(myCmd, ', col = rainbow(argLen)')
}else{
    myCmd = paste0(myCmd, ', col = color')
}
if(fill == 'rainbow'){
    myCmd = paste0(myCmd, ', fill = rainbow(argLen)')
}else{
    myCmd = paste0(myCmd, ', fill = fill')
}
if(exists('cNames')) myCmd = paste0(myCmd, ', category.names = cNames')
if(cColor == 'rainbow'){
    myCmd = paste0(myCmd, ', cat.col = rainbow(argLen)')
}else{
    myCmd = paste0(myCmd, ', cat.col = cColor')
}
if(exists('cPos')) myCmd = paste0(myCmd, ', cat.pos = cPos')

if(exists('lty')) myCmd = paste0(myCmd, ', lty = lty')
if(exists('rDegree')) myCmd = paste0(myCmd, ', rotation.degree = rDegree')

if(exists('margin')) myCmd = paste0(myCmd, ', margin = margin')
if(exists('inverted')) myCmd = paste0(myCmd, ', inverted = TRUE')
if(exists('redun')) myCmd = paste0(myCmd, ', force.unique = FALSE')

myCmd = paste0(myCmd, ')')
eval(parse(text = myCmd))

grid.draw(T)

unlink('VennDiagram*.log', force = TRUE)