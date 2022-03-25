#!/bin/env Rscript
args <- commandArgs()
scriptPath = strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName = basename(scriptPath)
scriptDir = dirname(scriptPath)
args = args[-(1:5)]
source(paste0(scriptDir, '/common.R'))

usage = function(){
    cat(paste0("Usage: ", scriptName) )
    cat(" -p=outputName.pdf input1.data inpu2.data [input3.data ...]
Option:
    Common:
    -p|pdf          FILE    The output figure in pdf[figure.pdf]
    -w|width        INT     The figure width
    -height         INT     The figure height
    -m|main         STR     The main title
    -mainS          DOU     The size of main title[22 for ggplot]
    -x|xlab         STR     The xlab
    -y|ylab         STR     The ylab
    -xl|xlog        INT     Transform the X scale to INT base log
    -yl|ylog        INT     Transform the Y scale to INT base log
    -x1             INT     The xlim start
    -x2             INT     The xlim end
    -y1             INT     The ylim start[0 for -f]
    -y2             INT     The ylim end[1 for -f]
    -ng|noGgplot            Draw figure in the style of R base rather than ggplot
    -ho|-horizontal DOU     Draw a horizontal dash line at y=DOU
    -frac                   Draw the Y axis in fraction
    -h|help                 Show help
    
    ggplot specific:

    -legendT        STR     The title shown on legend[Group]

    -a|alpha        DOU     The alpha of bar body
    -alphaV         STR     The column name to apply alpha (V3, V4, ...)
    -alphaT         STR     The title of alpha legend[Alpha]
    -alphaTP        POS     The title position of alpha legend[horizontal: top, vertical:right]
    -alphaLP        POS     The label position of alpha legend[horizontal: top, vertical:right]
    -alphaD         STR     The direction of alpha legend (horizontal, vertical)
    -c|color        STR     The color of line
    -colorTP        POS     The title position of color legend[horizontal: top, vertical:right]
    -colorLP        POS     The label position of color legend[horizontal: top, vertical:right]
    -colorD         STR     The direction of color legend (horizontal, vertical)
    -l|linetype     INT     The line type
    -linetypeV      STR     The column name to apply linetype (V3, V4,...)
    -linetypeT      STR     The title of linetype legend[Line Type]
    -linetypeTP     POS     The title position of linetype legend[horizontal: top, vertical:right]
    -linetypeLP     POS     The label position of linetype legend[horizontal: top, vertical:right]
    -linetypeD      STR     The direction of linetype legend (horizontal, vertical)
    -s|size         DOU     The size of bar boundary
    -sizeV          STR     The column name to apply size (V3, V4,...)
    -sizeT          STR     The title of size legend[Size]
    -sizeTP         POS     The title position of size legend[horizontal: top, vertical:right]
    -sizeLP         POS     The label position of size legend[horizontal: top, vertical:right]
    -sizeD          STR     The direction of size legend (horizontal, vertical)
          
    -fp|flip            Flip the Y axis to horizontal
    -facet      STR     The facet type (facet_wrap, facet_grid)
    -facetM     STR     The facet model (eg: '. ~ V3', 'V3 ~ .', 'V3 ~ V4', '. ~ V3 + V4', ...)
    -facetScl   STR     The axis scale in each facet ([fixed], free, free_x or free_y)

    -noGuide                Don't show the legend guide
    -lgPos          POS     The legend position[horizontal: top, vertical:right]
    -lgPosX         [0,1]   The legend relative postion on X
    -lgPosY         [0,1]   The legend relative postion on Y
    -lgTtlS         INT     The legend title size[22]
    -lgTxtS         INT     The legend text size[20]
    -lgBox          STR     The legend box style (horizontal, vertical)

    -xPer               Show X label in percentage
    -yPer               Show Y label in percentage
    -xComma             Show X label number with comma seperator
    -yComma             Show Y label number with comma seperator
    -axisRatio  DOU     The fixed aspect ratio between y and x units

    -annoTxt    STRs    The comma-seperated texts to be annotated
    -annoTxtX   INTs    The comma-seperated X positions of text
    -annoTxtY   INTs    The comma-seperated Y positions of text
Skill:
    Legend title of alpha, color, etc can be set as the same to merge their guides
")
  q(save='no')
}

if(length(args) == 0) usage()

myPdf = 'figure.pdf'
alphaT = 'Alpha'
linetypeT = 'Line Type'
sizeT = 'Size'
lgTtlS = 22
lgTxtS = 20
showGuide = TRUE
mainS = 22
myLegendTitle = 'Group'

for(i in 1:length(args)){
    arg = args[i]

    if(arg == '-frac'){
        fraction = TRUE
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'ho(rizontal)?', 'ho')
    if(!is.null(tmp)){
        myHorizontal = tmp
        args[i] = NA
        next
    }
    
    tmp = parseArg(arg, 'legendT', 'legendT')
    if(!is.null(tmp)){
        myLegendTitle = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'a(lpha)?', 'a')
    if(!is.null(tmp)){
        myAlpha = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'c(olor)?', 'c')
    if(!is.null(tmp)){
        color = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'l(inetype)?', 'l')
    if(!is.null(tmp)){
        linetype = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 's(ize)?', 's')
    if(!is.null(tmp)){
        size = tmp
        args[i] = NA
        next
    }
    
    if(arg == '-noGuide'){
        showGuide = FALSE
    }
    tmp = parseArg(arg, 'lgPos', 'lgPos')
    if(!is.null(tmp)){
        lgPos = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'lgPosX', 'lgPosX')
    if(!is.null(tmp)){
        lgPosX = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'lgPosY', 'lgPosY')
    if(!is.null(tmp)){
        lgPosY = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'lgTtlS', 'lgTtlS')
    if(!is.null(tmp)){
        lgTtlS = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'lgTxtS', 'lgTxtS')
    if(!is.null(tmp)){
        lgTxtS = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'lgBox', 'lgBox')
    if(!is.null(tmp)){
        lgBox = tmp
        args[i] = NA
        next
    }
        
    if(arg == '-fp' || arg =='-flip'){
        flip = TRUE
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'facet', 'facet')
    if(!is.null(tmp)){
        myFacet = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'facetM', 'facetM')
    if(!is.null(tmp)){
        facetM = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'facetScl', 'facetScl')
    if(!is.null(tmp)){
        facetScl = tmp
        args[i] = NA
        next
    }
    if(arg == '-xPer'){
        xPer = TRUE
        args[i] = NA
        next
    }
    if(arg == '-yPer'){
        yPer = TRUE
        args[i] = NA
        next
    }
    if(arg == '-xComma'){
        xComma = TRUE
        args[i] = NA
        next
    }
    if(arg == '-yComma'){
        yComma = TRUE
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'axisRatio', 'axisRatio')
    if(!is.null(tmp)){
        axisRatio = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'annoTxt', 'annoTxt')
    if(!is.null(tmp)){
        annoTxt = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'annoTxtX', 'annoTxtX')
    if(!is.null(tmp)){
        annoTxtX = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'annoTxtY', 'annoTxtY')
    if(!is.null(tmp)){
        annoTxtY = tmp
        args[i] = NA
        next
    }
        
    if(arg == '-h' || arg == '-help') usage()
    tmp = parseArg(arg, 'p(df)?', 'p')
    if(!is.null(tmp)){
        myPdf = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'w(idth)?', 'w')
    if(!is.null(tmp)){
        width = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'height?', 'height')
    if(!is.null(tmp)){
        height = tmp
        args[i] = NA
        next
    }
    if(arg == '-ng' || arg == '-noGgplot'){
        noGgplot = TRUE
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'x1', 'x1')
    if(!is.null(tmp)){
        x1 = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'x2', 'x2')
    if(!is.null(tmp)){
        x2 = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'y1', 'y1')
    if(!is.null(tmp)){
        y1 = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'y2', 'y2')
    if(!is.null(tmp)){
        y2 = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'xl(og)?', 'xl')
    if(!is.null(tmp)){
        xLog = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'yl(og)?', 'yl')
    if(!is.null(tmp)){
        yLog = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'm(ain)?', 'm')
    if(!is.null(tmp)){
        main = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'mainS', 'mainS')
    if(!is.null(tmp)){
        mainS = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'x(lab)?', 'x')
    if(!is.null(tmp)){
        xLab = tmp
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'y(lab)?', 'y')
    if(!is.null(tmp)){
        yLab = tmp
        args[i] = NA
        next
    }
}

args = args[!is.na(args)]
if(length(args) < 2) stop('Please specify two input files at least')


if(exists('width') && !exists('height')){
    pdf(myPdf, width = width, height = width * 0.6)
}else if(!exists('width') && exists('height')){
    pdf(myPdf, width = height * 1.6, height = height)
}else if(exists('width') && exists('height')){
    pdf(myPdf, width = width, height = height)
}else{
    pdf(myPdf)
}

library(tools)

fileNames = basename(file_path_sans_ext(args))
if(exists('noGgplot')){
    myColors = rainbow(length(args))
    data = read.delim(args[1], header = F)
    if(exists('fraction')) data[2] = data[[2]]/sum(data[2])
    myCmd = 'plot(data, type = "l", col = myColors[1]'
    if(exists('x1') && exists('x2')) myCmd = paste0(myCmd, ', xlim = c(x1, x2)')
    if(exists('myMain')) myCmd = paste0(myCmd, ', main = myMain')
    if(exists('myXlab')) myCmd = paste0(myCmd, ', xlab = myXlab')
    if(exists('myYlab')) myCmd = paste0(myCmd, ', ylab = myYlab')
    myCmd = paste0(myCmd, ')')
    eval(parse(text = myCmd))
    for(i in 2:length(args)){
        file = args[i]
        data = read.delim(file, header = F)
        if(exists('fraction')) data[2] = data[[2]] / sum(data[2])
        lines(data, col = myColors[i])
    }
    legend('topright', legend = fileNames, col = myColors, lty = 1, lwd = 1.5)
}else{
    library(ggplot2)
    data = cbind(read.delim(args[1], header = F), Group = fileNames[1])
    if(exists('fraction')) data[2] = data[[2]] / sum(data[2])
    for(i in 2:length(args)){
        file = args[i]
        fileName = fileNames[i]
        newData = cbind(read.delim(file, header = F), Group = fileName)
        if(exists('fraction')) newData[2] = newData[[2]] / sum(newData[2])
        data = rbind(data, newData)
    }
    p = ggplot(data, aes(x = V1, y = V2, color = Group)) + guides(color = guide_legend(title = myLegendTitle))
    myCmd = 'p = p + geom_line(show_guide = showGuide'
    if(exists('myAlpha')) myCmd = paste0(myCmd, ', alpha = myAlpha')
    if(exists('color')) myCmd = paste0(myCmd, ', color = color')
    if(exists('size')) myCmd = paste0(myCmd, ', size = size')
    myCmd = paste0(myCmd, ')')
    if(exists('myFacet')){
        myCmd = paste0(myCmd, ' + ', myFacet, '("' + facetM + '"')
        if(exists('facetScl')) myCmd = paste0(myCmd, ', scale = facetScl')
        myCmd = paste0(myCmd, ')')
    }
    eval(parse(text = myCmd))

    if(exists('lgPos')) p = p + theme(legend.position = lgPos)
    if(exists('lgPosX') && exists('lgPosY')) p = p + theme(legend.position = c(lgPosX, lgPosY))
    p = p + theme(legend.title = element_text(size = lgTtlS), legend.text = element_text(size = lgTxtS))
    if(exists('lgBox')) p = p + theme(legend.box = lgBox)
    if(exists('xPer')) p = p + scale_x_continuous(labels = percent)
    if(exists('yPer')) p = p + scale_y_continuous(labels = percent)
    if(exists('xComma')) p = p + scale_x_continuous(labels = comma)
    if(exists('yComma')) p = p + scale_y_continuous(labels = comma)
    if(exists('axisRatio')) p = p + coord_fixed(ratio = axisRatio)
    if(exists('annoTxt')) p = p + annotate('text', x = as.numeric(strsplit(annoTxtX, ',', fixed = T)),
                                           y = as.numeric(strsplit(annoTxtY, ',', fixed = T)),
                                           label = strsplit(annoTxt, ',', fixed = T))
    
    if(exists('x1') && exists('x2')) p = p + coord_cartesian(xlim = c(x1, x2))
    if(exists('y1') && exists('y2')) p = p + coord_cartesian(ylim = c(y1, y2))
    if(exists('x1') && exists('x2') && exists('y1') && exists('y2')) p = p + coord_cartesian(xlim = c(x1, x2), ylim = c(y1, y2))
    if(exists('xLog') || exists('yLog')){
        library(scales)
        if(exists('xLog')) p = p + scale_x_continuous(trans = log_trans(xLog)) + annotation_logticks(sides = 'b')
        if(exists('yLog')) p = p + scale_y_continuous(trans = log_trans(yLog)) + annotation_logticks(sides = 'l')
        p = p + theme(panel.grid.minor = element_blank())
    }
    if(exists('main')) p = p + ggtitle(main)
    p = p + theme(plot.title = element_text(size = mainS))
    if(exists('xLab')) p = p + xlab(xLab)
    p = p + theme(axis.title.x = element_text(size = mainS*0.8), axis.text.x = element_text(size = mainS*0.7))
    if(exists('yLab')) p = p + ylab(yLab)
    p = p + theme(axis.title.y = element_text(size = mainS*0.8), axis.text.y = element_text(size = mainS*0.7))
    if(exists('myHorizontal')) p = p + geom_hline(yintercept = myHorizontal, linetype = "longdash", size = 0.3)
    
    p
}
