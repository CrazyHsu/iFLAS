#!/usr/bin/env Rscript
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
    Common:
    -p|pdf          FILE    The output figure in pdf[figure.pdf]
    -w|width        INT     The figure width
       height       INT     The figure height
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
    -xLblS          DOU     The X-axis label size[20 for ggplot]
    -xTxtS          DOU     The X-axis text size[18 for ggplot]
    -yTxtS          DOU     The Y-axis text size[18 for ggplot]
    -xAngle         [0,360] The angle of tick labels[0]
    -ng|noGgplot            Draw figure in the style of R base rather than ggplot
    -ho|horizontal  DOU     Draw a horizontal line
    -t|tp                   Transpose the input table
    -dL|drawLine            Use the values in last line to draw a line across the boxes
    -h|help                 Show help

    R base specific:
    -no|noOutlier           Don't draw outlier

    ggplot specific:
    -stack                  Stack the data
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
    q(save = 'no')
}

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
xLblS = 20
xTxtS = 18
yTxtS = 18
xAngle = 0
vJust = 0.5
outlier = TRUE
notch = FALSE
oColor = 'red'

if(length(args) >= 1){
    for(i in 1:length(args)){
        arg = args[i]
        
        if(arg == '-dL' || arg == '-drawLine') drawLine = TRUE
        if(arg == '-t' || arg == '-tp') tp = TRUE
        if(arg == '-stack') toStack = TRUE
        if(arg == '-no' || arg == '-noOutlier') outlier = FALSE
        tmp = parseArgAsNum(arg, 'ho(rizontal)?', 'ho'); if(!is.null(tmp)) horizontal = tmp
        if(arg == '-nJ' || arg == '-nJitter') noJitter = TRUE
        if(arg == '-n' || arg == '-notch') notch = TRUE
        tmp = parseArg(arg, 'oC(olor)?', 'oC'); if(!is.null(tmp)) oColor = tmp
        tmp = parseArgAsNum(arg, 'oS(ize)?', 'oS'); if(!is.null(tmp)) oSize = tmp

        tmp = parseArgAsNum(arg, 'a(lpha)?', 'a'); if(!is.null(tmp)) myAlpha = tmp
        tmp = parseArg(arg, 'alphaV', 'alphaV'); if(!is.null(tmp)) alphaV = tmp
        tmp = parseArg(arg, 'alphaT', 'alphaT'); if(!is.null(tmp)) alphaT = tmp
        tmp = parseArg(arg, 'alphaTP', 'alphaTP'); if(!is.null(tmp)) alphaTP = tmp
        tmp = parseArg(arg, 'alphaLP', 'alphaLP'); if(!is.null(tmp)) alphaLP = tmp
        tmp = parseArg(arg, 'alphaD', 'alphaD'); if(!is.null(tmp)) alphaD = tmp
        tmp = parseArg(arg, 'c(olor)?', 'c'); if(!is.null(tmp)) color = tmp
        tmp = parseArg(arg, 'colorV', 'colorV'); if(!is.null(tmp)) colorV = tmp
        if(arg == '-colorC') colorC = TRUE
        tmp = parseArg(arg, 'colorT', 'colorT'); if(!is.null(tmp)) colorT = tmp
        tmp = parseArg(arg, 'colorTP', 'colorTP'); if(!is.null(tmp)) colorTP = tmp
        tmp = parseArg(arg, 'colorLP', 'colorLP'); if(!is.null(tmp)) colorLP = tmp
        tmp = parseArg(arg, 'colorD', 'colorD'); if(!is.null(tmp)) colorD = tmp
        tmp = parseArgAsNum(arg, 'f(ill)?', 'l'); if(!is.null(tmp)) fill = tmp
        tmp = parseArg(arg, 'f(ill)?', 'f'); if(!is.null(tmp)) fill = tmp
        tmp = parseArg(arg, 'fillV', 'fillV'); if(!is.null(tmp)) fillV = tmp
        tmp = parseArg(arg, 'fillT', 'fillT'); if(!is.null(tmp)) fillT = tmp
        tmp = parseArg(arg, 'fillTP', 'fillTP'); if(!is.null(tmp)) fillTP = tmp
        tmp = parseArg(arg, 'fillLP', 'fillLP'); if(!is.null(tmp)) fillLP = tmp
        tmp = parseArg(arg, 'fillD', 'fillD'); if(!is.null(tmp)) fillD = tmp
        tmp = parseArgAsNum(arg, 'l(inetype)?', 'l'); if(!is.null(tmp)) linetype = tmp
        tmp = parseArg(arg, 'linetypeV', 'linetypeV'); if(!is.null(tmp)) linetypeV = tmp
        tmp = parseArg(arg, 'linetypeT', 'linetypeT'); if(!is.null(tmp)) linetypeT = tmp
        tmp = parseArg(arg, 'linetypeTP', 'linetypeTP'); if(!is.null(tmp)) linetypeTP = tmp
        tmp = parseArg(arg, 'linetypeLP', 'linetypeLP'); if(!is.null(tmp)) linetypeLP = tmp
        tmp = parseArg(arg, 'linetypeD', 'linetypeD'); if(!is.null(tmp)) linetypeD = tmp
        tmp = parseArgAsNum(arg, 'shape', 'shape'); if(!is.null(tmp)) shape = tmp
        tmp = parseArg(arg, 'shapeV', 'shapeV'); if(!is.null(tmp)) shapeV = tmp
        tmp = parseArg(arg, 'shapeT', 'shapeT'); if(!is.null(tmp)) shapeT = tmp
        tmp = parseArg(arg, 'shapeTP', 'shapeTP'); if(!is.null(tmp)) shapeTP = tmp
        tmp = parseArg(arg, 'shapeLP', 'shapeLP'); if(!is.null(tmp)) shapeLP = tmp
        tmp = parseArg(arg, 'shapeD', 'shapeD'); if(!is.null(tmp)) shapeD = tmp
        tmp = parseArgAsNum(arg, 's(ize)?', 's'); if(!is.null(tmp)) size = tmp
        tmp = parseArg(arg, 'sizeV', 'sizeV'); if(!is.null(tmp)) sizeV = tmp
        tmp = parseArg(arg, 'sizeT', 'sizeT'); if(!is.null(tmp)) sizeT = tmp
        tmp = parseArg(arg, 'sizeTP', 'sizeTP'); if(!is.null(tmp)) sizeTP = tmp
        tmp = parseArg(arg, 'sizeLP', 'sizeLP'); if(!is.null(tmp)) sizeLP = tmp
        tmp = parseArg(arg, 'sizeD', 'sizeD'); if(!is.null(tmp)) sizeD = tmp
        tmp = parseArgAsNum(arg, 'weight', 'weight'); if(!is.null(tmp)) weight = tmp
        tmp = parseArg(arg, 'weightV', 'weightV'); if(!is.null(tmp)) weightV = tmp
        tmp = parseArg(arg, 'weightT', 'weightT'); if(!is.null(tmp)) weightT = tmp
        tmp = parseArg(arg, 'weightTP', 'weightTP'); if(!is.null(tmp)) weightTP = tmp
        tmp = parseArg(arg, 'weightLP', 'weightLP'); if(!is.null(tmp)) weightLP = tmp
        tmp = parseArg(arg, 'weightD', 'weightD'); if(!is.null(tmp)) weightD = tmp
        
        if(arg == '-noGuide') showGuide = FALSE
        tmp = parseArg(arg, 'lgPos', 'lgPos'); if(!is.null(tmp)) lgPos = tmp
        tmp = parseArgAsNum(arg, 'lgPosX', 'lgPosX'); if(!is.null(tmp)) lgPosX = tmp
        tmp = parseArgAsNum(arg, 'lgPosY', 'lgPosY'); if(!is.null(tmp)) lgPosY = tmp
        tmp = parseArgAsNum(arg, 'lgTtlS', 'lgTtlS'); if(!is.null(tmp)) lgTtlS = tmp
        tmp = parseArgAsNum(arg, 'lgTxtS', 'lgTxtS'); if(!is.null(tmp)) lgTxtS = tmp
        tmp = parseArg(arg, 'lgBox', 'lgBox'); if(!is.null(tmp)) lgBox = tmp
        
        if(arg == '-fp' || arg =='-flip') flip = TRUE
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
        tmp = parseArgAsNum(arg, 'h(eight)?', 'w'); if(!is.null(tmp)) height = tmp
        if(arg == '-ng' || arg == '-noGgplot') noGgplot = TRUE
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
        tmp = parseArgAsNum(arg, 'xAngle', 'xAngle'); if(!is.null(tmp)) xAngle = tmp
        tmp = parseArgAsNum(arg, 'vJust', 'vJust'); if(!is.null(tmp)) vJust = tmp
        tmp = parseArgAsNum(arg, 'xLblS', 'xLblS'); if(!is.null(tmp)) xLblS = tmp
        tmp = parseArgAsNum(arg, 'xTxtS', 'xTxtS'); if(!is.null(tmp)) xTxtS = tmp
        tmp = parseArgAsNum(arg, 'yTxtS', 'yTxtS'); if(!is.null(tmp)) yTxtS = tmp
    }
}

myCmd = 'pdf(myPdf'
if(exists('width')) myCmd = paste0(myCmd, ', width = width')
if(exists('height')) myCmd = paste0(myCmd, ', height = height')
myCmd = paste0(myCmd, ')')
eval(parse(text = myCmd))

data = read.delim(file('stdin'), header = F)
if(exists('tp')) data = t(data)
if(exists('drawLine')){
    line.df = as.data.frame(data[nrow(data),])
    colnames(line.df) = c('V1')
    data = as.data.frame(data[-nrow(data), ])
}
if(exists('noGgplot')){
    myCmd = 'boxplot(data, outline = outlier'
    if(exists('main')) myCmd = paste0(myCmd, ', main = main')
    if(exists('xLab')) myCmd = paste0(myCmd, ', xlab = xLab')
    if(exists('yLab')) myCmd = paste0(myCmd, ', ylab = yLab')
    myCmd = paste0(myCmd, ')')
    eval(parse(text = myCmd))
    
    if(exists('drawLine')) lines(line.df)
    if(exists('horizontal')) abline(h = horizontal, lty = 2)
}else{
    library(ggplot2)
    library(scales) # for xPer, yPer
    if(exists('toStack')){
        data = stack(data)
        colnames(data) = c('V2', 'V1')
    }
    p = ggplot(data, aes(factor(V1, levels = unique(V1)), V2))
    if(exists('alphaV')){
        p = p + aes_string(alpha = alphaV)
        myCmd = 'p = p + guides(alpha = guide_legend(alphaT'
        if(exists('alphaTP')) myCmd = paste0(myCmd, ', title.position = alphaTP')
        if(exists('alphaLP')) myCmd = paste0(myCmd, ', label.position = alphaLP')
        if(exists('alphaD')) myCmd = paste0(myCmd, ', direction = alphaD')
        myCmd = paste0(myCmd, '))')
        eval(parse(text = myCmd))
    }
    if(exists('colorV')){
        if(exists('colorC')){
            p = p + aes_string(color = colorV)
        }else{
            myCmd = paste0('p = p + aes(color = factor(', colorV, '))'); eval(parse(text = myCmd))
        }
        myCmd = 'p = p + guides(color = guide_legend(colorT'
        if(exists('colorTP')) myCmd = paste0(myCmd, ', title.position = colorTP')
        if(exists('colorLP')) myCmd = paste0(myCmd, ', label.position = colorLP')
        if(exists('colorD')) myCmd = paste0(myCmd, ', direction = colorD')
        myCmd = paste0(myCmd, '))')
        eval(parse(text = myCmd))
    }
    if(exists('fillV')){
        myCmd = paste0('p = p + aes(fill = factor(', fillV, '))'); eval(parse(text = myCmd))
        myCmd = 'p = p + guides(fill = guide_legend(fillT'
        if(exists('fillTP')) myCmd = paste0(myCmd, ', title.position = fillTP')
        if(exists('fillLP')) myCmd = paste0(myCmd, ', label.position = fillLP')
        if(exists('fillD')) myCmd = paste0(myCmd, ', direction = fillD')
        myCmd = paste0(myCmd, '))')
        eval(parse(text = myCmd))
    }
    if(exists('linetypeV')){
        myCmd = paste0('p = p + aes(linetype = factor(', linetypeV, '))'); eval(parse(text = myCmd))
        myCmd = 'p = p + guides(linetype = guide_legend(linetypeT'
        if(exists('linetypeTP')) myCmd = paste0(myCmd, ', title.position = linetypeTP')
        if(exists('linetypeLP')) myCmd = paste0(myCmd, ', label.position = linetypeLP')
        if(exists('linetypeD')) myCmd = paste0(myCmd, ', direction = linetypeD')
        myCmd = paste0(myCmd, '))')
        eval(parse(text = myCmd))
    }
    if(exists('shapeV')){
        myCmd = paste0('p = p + aes(shape = factor(', shapeV, '))'); eval(parse(text = myCmd))
        myCmd = 'p = p + guides(shape = guide_legend(shapeT'
        if(exists('shapeTP')) myCmd = paste0(myCmd, ', title.position = shapeTP')
        if(exists('shapeLP')) myCmd = paste0(myCmd, ', label.position = shapeLP')
        if(exists('shapeD')) myCmd = paste0(myCmd, ', direction = shapeD')
        myCmd = paste0(myCmd, '))')
        eval(parse(text = myCmd))
    }
    if(exists('sizeV')){
        p = p + aes_string(size = sizeV)
        myCmd = 'p = p + guides(size = guide_legend(sizeT'
        if(exists('sizeTP')) myCmd = paste0(myCmd, ', title.position = sizeTP')
        if(exists('sizeLP')) myCmd = paste0(myCmd, ', label.position = sizeLP')
        if(exists('sizeD')) myCmd = paste0(myCmd, ', direction = sizeD')
        myCmd = paste0(myCmd, '))')
        eval(parse(text = myCmd))
    }
    if(exists('weightV')){
        p = p + aes_string(weight = weightV)
        myCmd = 'p = p + guides(weight = guide_legend(weightT'
        if(exists('weightTP')) myCmd = paste0(myCmd, ', title.position = weightTP')
        if(exists('weightLP')) myCmd = paste0(myCmd, ', label.position = weightLP')
        if(exists('weightD')) myCmd = paste0(myCmd, ', direction = weightD')
        myCmd = paste0(myCmd, '))')
        eval(parse(text = myCmd))
    }
    
    myCmd = paste0('p = p + geom_boxplot(show.legend = showGuide, notch = notch, outlier.colour = oColor')
    if(exists('oSize')) myCmd = paste0(myCmd, ', outlier.size = oSize')
    
    if(exists('myAlpha')) myCmd = paste0(myCmd, ', alpha = myAlpha')
    if(exists('color')) myCmd = paste0(myCmd, ', color = color')
    if(exists('fill')) myCmd = paste0(myCmd, ', fill = fill')
    if(exists('shape')) myCmd = paste0(myCmd, ', shape = shape')
    if(exists('size')) myCmd = paste0(myCmd, ', size = size')
    myCmd = paste0(myCmd, ')')
    if(exists('myFacet')){
        myCmd = paste0(myCmd, ' + ', myFacet, '("', facetM, '"')
        if(exists('facetScl')) myCmd = paste0(myCmd, ', scale = facetScl')
        myCmd = paste0(myCmd, ')')
    }
    eval(parse(text = myCmd))
    
    if(!exists('noJitter')) p = p + geom_jitter(alpha = 0.2)
    
    if(exists('flip')) p = p + coord_flip()
    
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
    
    if(exists('y1') && exists('y2')) p = p + coord_cartesian(ylim = c(y1, y2))
    if(exists('xLog') || exists('yLog')){
        library(scales)
        if(exists('xLog')) p = p + scale_x_continuous(trans = log_trans(xLog)) + annotation_logticks(sides = 'b')
        if(exists('yLog')) p = p + scale_y_continuous(trans = log_trans(yLog)) + annotation_logticks(sides = 'l')
        p = p + theme(panel.grid.minor = element_blank())
    }
    if(exists('main')) p = p + ggtitle(main)
    p = p + theme(plot.title = element_text(size = mainS, hjust = 0.5))
    if(exists('xLab')) p = p + xlab(xLab) + theme(axis.title.x = element_text(size = xLblS), 
                                                  axis.text.x = element_text(size = xTxtS, angle = xAngle, vjust = vJust))
    if(exists('yLab')) p = p + ylab(yLab)
    p = p + theme(axis.title.y = element_text(size = mainS*0.8), axis.text.y = element_text(size = yTxtS))
    
    if(exists('horizontal')) p = p + geom_hline(yintercept = horizontal, linetype = "longdash", size = 0.3)
    if(exists('drawLine')) p = p + geom_line(data = line.df, aes(x = 1:nrow(line.df), y = V1))
    p
}