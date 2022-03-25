#!/usr/bin/env Rscript
args <- commandArgs()
scriptPath = strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName = basename(scriptPath)
scriptDir = dirname(scriptPath)
args = args[-(1:5)]
source(paste0(scriptDir, '/common.R'))
library(tools)

usage = function(){
    cat(paste0("Usage: ", scriptName) )
    cat(" -p=outputName.pdf input1.tsv inpu2.tsv [input3.tsv ...]
Option:
    Common:
    -p|pdf          FILE    The output figure in pdf[figure.pdf]
    -w|width        INT     The figure width
    -height         INT     The figure height
    -m|main         STR     The main title
    -mainS          DOU     The size of main title[20 for ggplot]
    -x|xlab         STR     The xlab
    -y|ylab         STR     The ylab
    -xl|xlog        INT     Transform the X scale to INT base log
    -yl|ylog        INT     Transform the Y scale to INT base log
    -x1             INT     The xlim start
    -x2             INT     The xlim end
    -y1             INT     The ylim start
    -y2             INT     The ylim end
    -ng|noGgplot            Draw figure in the style of R base rather than ggplot
    -h|help                 Show help

    R base specific:
    -no|noOutlier           Don't draw outlier
    
    ggplot specific:
    -ho|horizontal  DOU     Draw a horizontal line
    -nJ|nJitter             Do not draw jitter
    -n|notch                Draw notch
    -oC|oColor      STR     Outlier color
    -oS|oSize       DOU     Outlier size

    -a|alpha        DOU     The alpha of box body
    -alphaV         STR     The column name to apply alpha (V3, V4, ...)
    -alphaT         STR     The title of alpha legend[Alpha]
    -alphaTP        POS     The title position of alpha legend[horizontal: top, vertical:right]
    -alphaLP        POS     The label position of alpha legend[horizontal: top, vertical:right]
    -alphaD         STR     The direction of alpha legend (horizontal, vertical)
    -c|color        STR     The color of box boundary
    -colorV         STR     The column name to apply color (V3, V4,...)
    -colorC                 Continuous color mapping
    -colorT         STR     The title of color legend[Color]
    -colorTP        POS     The title position of color legend[horizontal: top, vertical:right]
    -colorLP        POS     The label position of color legend[horizontal: top, vertical:right]
    -colorD         STR     The direction of color legend (horizontal, vertical)
    -f|fill         STR     The color of box body
    -fillV          STR     The column name to apply fill (V3, V4,...)
    -fillT          STR     The title of fill legend[Fill]
    -fillTP         POS     The title position of fill legend[horizontal: top, vertical:right]
    -fillLP         POS     The label position of fill legend[horizontal: top, vertical:right]
    -fillD          STR     The direction of fill legend (horizontal, vertical)
    -l|linetype     INT     The line type
    -linetypeV      STR     The column name to apply linetype (V3, V4,...)
    -linetypeT      STR     The title of linetype legend[Line Type]
    -linetypeTP     POS     The title position of linetype legend[horizontal: top, vertical:right]
    -linetypeLP     POS     The label position of linetype legend[horizontal: top, vertical:right]
    -linetypeD      STR     The direction of linetype legend (horizontal, vertical)
    -shape          STR     The shape of box
    -shapeV         STR     The column name to apply shape (V3, V4,...)
    -shapeT         STR     The title of shape legend[Shape]
    -shapeTP        POS     The title position of shape legend[horizontal: top, vertical:right]
    -shapeLP        POS     The label position of shape legend[horizontal: top, vertical:right]
    -shapeD         STR     The direction of shape legend (horizontal, vertical)
    -s|size         DOU     The size of box body
    -sizeV          STR     The column name to apply size (V3, V4,...)
    -sizeT          STR     The title of size legend[Size]
    -sizeTP         POS     The title position of size legend[horizontal: top, vertical:right]
    -sizeLP         POS     The label position of size legend[horizontal: top, vertical:right]
    -sizeD          STR     The direction of size legend (horizontal, vertical)
    -weight         DOU     The weight of box
    -weightV        STR     The column name to apply weight (V3, V4,...)
    -weightT        STR     The title of weight legend[weight]
    -weightTP       POS     The title position of weight legend[horizontal: top, vertical:right]
    -weightLP       POS     The label position of weight legend[horizontal: top, vertical:right]
    -weightD        STR     The direction of weight legend (horizontal, vertical)
                    
    -noGuide                Don't show the legend guide
    -lgPos          POS     The legend position[horizontal: top, vertical:right]
    -lgPosX         [0,1]   The legend relative postion on X
    -lgPosY         [0,1]   The legend relative postion on Y
    -lgTtlS         INT     The legend title size[15]
    -lgTxtS         INT     The legend text size[15]
    -lgBox          STR     The legend box style (horizontal, vertical)

    -fp|flip                Flip the Y axis to horizontal
    -facet          STR     The facet type (facet_wrap, facet_grid)
    -facetM         STR     The facet model (eg: '. ~ V3', 'V3 ~ .', 'V3 ~ V4', '. ~ V3 + V4', ...)
    -facetScl       STR     The axis scale in each facet ([fixed], free, free_x or free_y)

    -xPer                   Show X label in percentage
    -yPer                   Show Y label in percentage
    -xComma                 Show X label number with comma seperator
    -yComma                 Show Y label number with comma seperator
    -axisRatio      DOU     The fixed aspect ratio between y and x units

    -annoTxt        STRs    The comma-seperated texts to be annotated
    -annoTxtX       INTs    The comma-seperated X positions of text
    -annoTxtY       INTs    The comma-seperated Y positions of text
    
Skill:
    Legend title of alpha, color, etc can be set as the same to merge their guides
")
  q(save='no')
}

if(length(args) == 0) usage()

alphaT = 'Alpha'
colorT = 'Color'
fillT = 'Fill'
linetypeT = 'Line Type'
shape = 'Shape'
sizeT = 'Size'
weightT = 'Weight'
lgTtlS = 15
lgTxtS = 15
showGuide = TRUE
myPdf = 'figure.pdf'
mainS = 20
outlier = TRUE
notch = FALSE
oColor = 'red'

for(i in 1:length(args)){
    arg = args[i]
    if(arg == '-no' || arg == '-noOutlier'){
        outlier = FALSE
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'ho(rizontal)?', 'ho')
    if(!is.null(tmp)){
        horizontal = tmp
        args[i] = NA
        next
    }
    if(arg == '-nJ' || arg == '-nJitter'){
        noJitter = TRUE
        args[i] = NA
        next
    }
    if(arg == '-n' || arg == '-notch'){
        notch = TRUE
        args[i] = NA
        next
    }
    tmp = parseArg(arg, 'oC(olor)?', 'oC')
    if(!is.null(tmp)){
        oColor = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'oS(ize)?', 'oS')
    if(!is.null(tmp)){
        oSize = tmp
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
    tmp = parseArgAsNum(arg, 'f(ill)?', 'f')
    if(!is.null(tmp)){
        fill = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 'shape', 'shape')
    if(!is.null(tmp)){
        shape = tmp
        args[i] = NA
        next
    }
    tmp = parseArgAsNum(arg, 's(ize)?', 's')
    if(!is.null(tmp)){
        size = tmp
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
    tmp = parseArgAsNum(arg, 'height', 'height')
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

fileNames = basename(file_path_sans_ext(args))
data = data.frame()
for(i in 1:length(args)){
    file = args[i]
    newData = read.delim(file, header = F)
    names(newData)[1] = fileNames[i]
    data = rbind(data, stack(newData))
}

if(exists('noGgplot')){
    myCmd = 'boxplot(values~ind, data = data, outline = outlier'
    if(exists('x1') && exists('x2')) myCmd = paste0(myCmd, ', xlim = c(x1, x2)')
    if(exists('main')) myCmd = paste0(myCmd, ', main = main')
    if(exists('xLab')) myCmd = paste0(myCmd, ', xlab = xLab')
    if(exists('yLab')) myCmd = paste0(myCmd, ', ylab = yLab')
    myCmd = paste0(myCmd, ', cex.axis = mainS * 0.1')
    myCmd = paste0(myCmd, ', cex.lab = mainS * 0.1')
    myCmd = paste0(myCmd, ', pars=list(par(mar=c(5,5,2,2)))')
    myCmd = paste0(myCmd, ')')
    eval(parse(text = myCmd))
}else{
    library(ggplot2)
    p = ggplot(data, aes(x = factor(ind), y = values))
    
    myCmd = paste0('p = p + geom_boxplot(notch = notch, outlier.colour = oColor')
    if(exists('oSize')) myCmd = paste0(myCmd, ', outlier.size = oSize')
    
    if(exists('myAlpha')) myCmd = paste0(myCmd, ', alpha = myAlpha')
    if(exists('color')) myCmd = paste0(myCmd, ', color = color')
    if(exists('fill')) myCmd = paste0(myCmd, ', fill = fill')
    if(exists('shape')) myCmd = paste0(myCmd, ', shape = shape')
    if(exists('size')) myCmd = paste0(myCmd, ', size = size')
    myCmd = paste0(myCmd, ')')
    if(exists('myFacet')){
        myCmd = paste0(myCmd, ' + ', myFacet, '("' + facetM + '"')
        if(exists('facetScl')) myCmd = paste0(myCmd, ', scale = facetScl')
        myCmd = paste0(myCmd, ')')
    }
    eval(parse(text = myCmd))
    
    if(!exists('noJitter')) p = p + geom_jitter()
    
    if(exists('flip')) p = p + coord_flip()
    
    if(exists('xPer')) p = p + scale_x_continuous(labels = percent)
    if(exists('yPer')) p = p + scale_y_continuous(labels = percent)
    if(exists('xComma')) p = p + scale_x_continuous(labels = comma)
    if(exists('yComma')) p = p + scale_y_continuous(labels = comma)
    if(exists('axisRatio')) p = p + coord_fixed(ratio = axisRatio)
    if(exists('annoTxt')) p = p + annotate('text', x = as.numeric(strsplit(annoTxtX, ',', fixed = T)),
                                           y = as.numeric(strsplit(annoTxtY, ',', fixed = T)),
                                           label = strsplit(annoTxt, ',', fixed = T))
    
    if(exists('y1') && exists('y2')) p = p + coord_cartesian(ylim = c(y1, y2))
    if(exists('xLog') || exists('yLog')){
        library(scales)
        if(exists('xLog')) p = p + scale_x_continuous(trans = log_trans(xLog)) + annotation_logticks(sides = 'b')
        if(exists('yLog')) p = p + scale_y_continuous(trans = log_trans(yLog)) + annotation_logticks(sides = 'l')
        p = p + theme(panel.grid.minor = element_blank())
    }
    if(exists('main')) p = p + ggtitle(main)
    p = p + theme(plot.title = element_text(size = mainS, hjust = 0.5))
    if(exists('xLab')) p = p + xlab(xLab) + theme(axis.title.x = element_text(size = mainS*0.8), axis.text.x = element_text(size = mainS*0.7))
    if(exists('yLab')) p = p + ylab(yLab)
    p = p + theme(axis.title.y = element_text(size = mainS*0.8), axis.text.y = element_text(size = mainS*0.7))
    
    if(exists('horizontal')) p = p + geom_hline(yintercept = horizontal, linetype = "longdash", size = 0.3)
    p
}