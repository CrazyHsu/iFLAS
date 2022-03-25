#!/bin/env Rscript
args <- commandArgs()
scriptPath = strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName = basename(scriptPath)
scriptDir = dirname(scriptPath)
args = args[-(1:5)]
source(paste0(scriptDir, '/common.R'))

usage = function(){
    cat(paste0("Usage: ", scriptName) )
    cat(" -p=outputName.pdf <input.tsv
Option:
    -p|pdf      FILE    The output figure in pdf[figure.pdf]
    -w|width    INT     The figure width
    -height     INT     The figure height
    -m|main     STR     The main title
    -mainS      DOU     The size of main title[20]
    -x|xlab     STR     The xlab
    -y|ylab     STR     The ylab
    -xl|xlog    INT     Transform the X scale to INT base log
    -yl|ylog    INT     Transform the Y scale to INT base log
    -x1         INT     The xlim start
    -x2         INT     The xlim end
    -y1         INT     The ylim start
    -y2         INT     The ylim end
    -h|help             Show help

    -b|bin      DOU     The bin number on X and Y
    -binWidthX  DOU     The bin width on X
    -binWidthY  DOU     The bin width on Y
    -fSum       STR     The summary funtion ([mean], identity, median, var, quantile, ...)

    -fL         STR     The low color of color filling
    -fH         STR     The high color of color filling

    -facet      STR     The facet type (facet_wrap, facet_grid)
    -facetM     STR     The facet model (eg: '. ~ V4', 'V4 ~ .', 'V4 ~ V5', '. ~ V4 + V5', ...)
    -facetScl   STR     The axis scale in each facet ([fixed], free, free_x or free_y)

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
	q(save = 'no')
}

lgTtlS = 15
lgTxtS = 15
showGuide = TRUE
myPdf = 'figure.pdf'
mainS = 20
bin = 30
fSum = 'mean'

if(length(args) >= 1){
    for(i in 1:length(args)){
        arg = args[i]
        
        tmp = parseArgAsNum(arg, 'b(in)?', 'bin'); if(!is.null(tmp)) bin = tmp
        tmp = parseArgAsNum(arg, 'binWidthX', 'binWidthX'); if(!is.null(tmp)) binWidthX = tmp
        tmp = parseArgAsNum(arg, 'binWidthY', 'binWidthY'); if(!is.null(tmp)) binWidthY = tmp
        tmp = parseArg(arg, 'fSum', 'fSum'); if(!is.null(tmp)) fSum = tmp
        
        tmp = parseArg(arg, 'fL(ow)?', 'fLow'); if(!is.null(tmp)) fillLow = tmp
        tmp = parseArg(arg, 'fH(igh)?', 'fHigh'); if(!is.null(tmp)) fillHigh = tmp
        
        tmp = parseArg(arg, 'facet', 'facet'); if(!is.null(tmp)) myFacet = tmp
        tmp = parseArg(arg, 'facetM', 'facetM'); if(!is.null(tmp)) facetM = tmp
        tmp = parseArg(arg, 'facetScl', 'facetScl'); if(!is.null(tmp)) facetScl = tmp
        if(arg == '-xPer') xPer = TRUE
        if(arg == '-yPer') yPer = TRUE
        if(arg == '-xComma') xComma = TRUE
        if(arg == '-yComma') yComma = TRUE
        tmp = parseArgAsNum(arg, 'axisRatio', 'axisRatio'); if(!is.null(tmp)) axisRatio = tmp
        tmp = parseArg(arg, 'annoTxt', 'annoTxt'); if(!is.null(tmp)) annoTxt = tmp
        tmp = parseArg(arg, 'annoTxtX', 'annoTxtX'); if(!is.null(tmp)) annoTxtX = tmp
        tmp = parseArg(arg, 'annoTxtY', 'annoTxtY'); if(!is.null(tmp)) annoTxtY = tmp
        
        if(arg == '-h' || arg == '-help') usage()
        tmp = parseArg(arg, 'p(df)?', 'p'); if(!is.null(tmp)) myPdf = tmp
        tmp = parseArgAsNum(arg, 'w(idth)?', 'w'); if(!is.null(tmp)) width = tmp
        tmp = parseArgAsNum(arg, 'height', 'height'); if(!is.null(tmp)) height = tmp
        tmp = parseArgAsNum(arg, 'x1', 'x1'); if(!is.null(tmp)) x1 = tmp
        tmp = parseArgAsNum(arg, 'x2', 'x2'); if(!is.null(tmp)) x2 = tmp
        tmp = parseArgAsNum(arg, 'y1', 'y1'); if(!is.null(tmp)) y1 = tmp
        tmp = parseArgAsNum(arg, 'y2', 'y2'); if(!is.null(tmp)) y2 = tmp
        tmp = parseArgAsNum(arg, 'xl(og)?', 'xl'); if(!is.null(tmp)) xLog = tmp
        tmp = parseArgAsNum(arg, 'yl(og)?', 'yl'); if(!is.null(tmp)) yLog = tmp
        tmp = parseArg(arg, 'm(ain)?', 'm'); if(!is.null(tmp)) main = tmp
        tmp = parseArgAsNum(arg, 'mainS', 'mainS'); if(!is.null(tmp)) mainS = tmp
        tmp = parseArg(arg, 'x(lab)?', 'x'); if(!is.null(tmp)) xLab = tmp
        tmp = parseArg(arg, 'y(lab)?', 'y'); if(!is.null(tmp)) yLab = tmp
    }
}
if(exists('width') && !exists('height')){
    pdf(myPdf, width = width, height = width * 0.6)
}else if(!exists('width') && exists('height')){
    pdf(myPdf, width = height * 1.6, height = height)
}else if(exists('width') && exists('height')){
    pdf(myPdf, width = width, height = height)
}else{
    pdf(myPdf)
}

data = read.delim(file('stdin'), header = F)

if(exists('noGgplot')){

}else{
    library(ggplot2)
    p = ggplot(data, aes(V1, V2, z = V3))
    
    myCmd = paste0('p = p + stat_summary2d(bins = bin')
    if(exists("binWidthX") && exists("binWidthY")) myCmd = paste0(myCmd, ', binwidth = c(binWidthX, binWidthY)')
    if(exists('fSum')) myCmd = paste0(myCmd, ', fun = ', fSum)
    if(exists('fill')) myCmd = paste0(myCmd, ', fill = fill')
    myCmd = paste0(myCmd, ')')
    
    if(exists('myFacet')){
        myCmd = paste0(myCmd, ' + ', myFacet, '("' + facetM + '"')
        if(exists('facetScl')) myCmd = paste0(myCmd, ', scale = facetScl')
        myCmd = paste0(myCmd, ')')
    }
    eval(parse(text = myCmd))
    
    if(exists('fillLow') && exists('fillHigh')) p = p + scale_fill_continuous(low=fillLow, high=fillHigh, name = fSum)
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
        if(exists('xLog')) p = p + scale_x_continuous(trans = log_trans(xLog))
        if(exists('yLog')) p = p + scale_y_continuous(trans = log_trans(yLog))
        p = p + annotation_logticks() + theme(panel.grid.minor = element_blank())
    }
    if(exists('main')) p = p + ggtitle(main)
    p = p + theme(plot.title = element_text(size = mainS, hjust = 0.5))
    if(exists('xLab')) p = p + xlab(xLab) + theme(axis.title.x = element_text(size = mainS*0.8), axis.text.x = element_text(size = mainS*0.7), legend.text = element_text(size = mainS*0.7), legend.title = element_text(size = mainS*0.7))
    if(exists('yLab')) p = p + ylab(yLab) + theme(axis.title.y = element_text(size = mainS*0.8), axis.text.y = element_text(size = mainS*0.7), legend.text = element_text(size = mainS*0.7), legend.title = element_text(size = mainS*0.7))
    p
}