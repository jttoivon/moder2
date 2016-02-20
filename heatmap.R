#!/usr/bin/Rscript

# FIX
# - format of the cell text should be possible to give on command line, FIXED
# - Make the cell font size adjust by the number of decimals in the cell,
#   and possibly by the number of cells per row



#library(RSVGDevice)
#require(graphics)
#require(gplots)
require(pheatmap)
require(RColorBrewer)
#require(Heatplus)
require(grid)

#is.integer.matrix <- function(x)
#{
#  all(x%%1 == 0)
#}

# detect odd/even integers
odd <- function(x) x!=as.integer(x/2)*2
even <- function(x) x==as.integer(x/2)*2 

# Generat a set of n colors which smoothly transition from 'low' to 'mid' to 'high'.
colorpanel <- function(n,low='green',mid='black',high='red')
  {
    if(even(n)) warning("n is even: colors panel will not be symmetric")

    # convert to rgb
    low <- col2rgb(low)
    mid <- col2rgb(mid)
    high <- col2rgb(high)

    # determine length of each component
    lower <- ceiling(n/2)
    upper <- ceiling(n/2)
    
    red  <- c(
              seq(low[1,1], mid [1,1], length=lower),
              seq(mid[1,1], high[1,1], length=upper)
              )/255

    red = red[-lower]

    green <- c(
               seq(low[3,1], mid [3,1], length=lower),
               seq(mid[3,1], high[3,1], length=upper)
               )/255
    green = green[-lower]

    blue <- c(
              seq(low[2,1], mid [2,1], length=lower),
              seq(mid[2,1], high[2,1], length=upper)
              )/255
    blue = blue[-lower]

             
    rgb(red,blue,green)
  }

colorpanel2 <- function(n,low='black',high='red')
  {

    # convert to rgb
    low <- col2rgb(low)
    high <- col2rgb(high)

    red  <- c(seq(low[1,1], high [1,1], length=n))/255
    green <- c(seq(low[3,1], high [3,1], length=n))/255
    blue <- c(seq(low[2,1], high [2,1], length=n))/255
             
    rgb(red,blue,green)
  }



# Generate red-black-green colorscale
redgreen <- function(n) colorpanel(n, 'red', 'black', 'green')

# Generate green-black-red colorscale
greenred <- function(n) colorpanel(n, 'green', 'black', 'red' )

# Generate blue white red colorscale
bluered  <- function(n) colorpanel(n, 'blue','white','red')

red <- function(n) colorpanel2(n, 'black', 'red')

#pie(rep(1,13), col=redgreen(13))

options=c("-p", "-c", "-s", "-w")

output.format="png"
usecounts=FALSE
usewidth=FALSE
title.given=FALSE
myformat="%.1f"
myfontsize=16
Args = commandArgs(TRUE)
print(Args)

while (length(Args) > 0 && substr(Args[1], 1, 1) == "-") {
    option=Args[1]
    if (option == "-p") {
        output.format="pdf"
        print("Using pdf")
        Args=Args[-1]
    }
    else if (option == "-c") {
        usecounts=TRUE
        print("Using counts instead of logodds")
        Args=Args[-1]
    }      
    else if (option == "-s") {
        output.format="svg"
        print("Using svg")
        Args=Args[-1]
    }
    else if (option == "-w") {  
        if (length(Args) < 2)
            stop("Option -w requires integer parameter")
        width=as.numeric(Args[2])
        if (is.na(width) || width <= 0) {
            stop(help.string)
        }
        usewidth=TRUE
        Args=Args[c(-1,-2)]
    }
    else if (option == "-z") {  
        if (length(Args) < 2)
            stop("Option -z requires a numeric parameter")
        myfontsize=as.numeric(Args[2])
        if (is.na(myfontsize) || myfontsize <= 0) {
            stop(help.string)
        }
        Args=Args[c(-1,-2)]
    }
    else if (option == "-f") {  
        if (length(Args) < 2)
            stop("Option -f requires a format string parameter")
        myformat=Args[2]
        Args=Args[c(-1,-2)]
    }
    else if (option == "-t") {  
        if (length(Args) < 2)
            stop("Option -t requires a string parameter")
        title=Args[2]
        title.given=TRUE
        Args=Args[c(-1,-2)]
    }
    else {
        stop("Unknown option ", option, "\n", paste(usage))
    }      
}



help.string=paste("Usage: ./heatmap.R [-p] [-c] [-s] [-w pwm_width] [ -f formatstring ] [ -t title ] inputfile [outputfile]",
	          "-p\tCreate pdf",
		  "-c\tUse counts instead of logodds", 
                  "-s\tCreate svg",
                  "-f formatstring\tFormat string used for number in cell", 
		  "-w width\tUse given width to number columns, no header or row names in inputfile",
		  "-z size\tUse given font size",
		  "-t title\tUse the given title",
		  sep='\n')

		    
if (length(Args) < 1) {
  stop(help.string)
}


input.file=Args[1]

if (length(Args) == 1) {
   dir=dirname(input.file)
   base=sub("\\.[^.]*$", "", basename(input.file))  # remove the extension
   ending = switch(output.format, pdf=".pdf", png=".png", svg=".svg")
   output.file=paste(dir, "/", base, ending, sep="")
} else {
   output.file = Args[2]
}

if (usewidth) {
   cob <- read.csv(input.file, sep="\t", check.names=FALSE, header=FALSE)  # No header and row name in inputfile
} else {
   cob <- read.csv(input.file, sep="\t", check.names=FALSE, header=TRUE, row.names=1)
}

cob
d=dim(cob)
d
rows=d[1]
cols=d[2]

cob <- cob[,1:length(cob)]

if (usewidth) {
   if (rows == 3) {
      row.names(cob) <- c("HT", "HH", "TT")
   } else {
      row.names(cob) <- c("HT", "HH", "TT", "TH")
   }

   min.value = -width+1
   names(cob) <- c(min.value:(cols+min.value-1))
}

cob_matrix <- data.matrix(cob)

#base_size <- 9

switch(output.format, 
 	pdf=pdf(output.file, 16.0, 4.0), 
	png=png(output.file, 800, 200), 
#	svg=svg(output.file, 16.0, 4.0, 
	svg=svg(output.file, 8.0, 2.0, 
	onefile=TRUE))

#devSVG(file='SVG Output.svg', height=6, width=6, onefile=TRUE)

#title=paste("Heatmappi", input.file)
if (!title.given) {
   title=input.file
}

heatmap.image <- function(cob_matrix, rows, cols)
{
  cob_heatmap = image(1:cols, 1:rows, t(cob_matrix)[,rows:1], axes=F, 
	    main=title, xlab="Gap length", ylab="Orientation")
  axis(2, seq(rows, 1), labels=row.names(cob), xpd=TRUE, las=1)
  axis(1, seq(1, cols), labels=names(cob), xpd=TRUE, las=1)
  legend("right", legend=c("4","6","8"))
}

my_draw_colnames = function(coln, ...){
        m = length(coln)
        x = (1:m)/m - 1/2/m
        grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 1.5, hjust = 0.5, rot = 0, gp = gpar(...))
}


heatmap.pretty <- function(cob_matrix, rows, cols)
{
  m = max(abs(cob_matrix))
#  m = 1.0
  assignInNamespace(x="draw_colnames", value="my_draw_colnames",
  ns=asNamespace("pheatmap"))


  my.color = if (usecounts) 
                 colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(200)[101:200]
             else
                 colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

  my.breaks = if (usecounts) seq(0, m, length.out=101) else seq(-m, m, length.out=101)
  mygray="#808080"
  myneutral=my.color[1]
  constant.color = rep(myneutral, 102)
  if (usecounts) {
      #my.color = c(my.color[1], my.color[1], my.color)
      my.color = c(myneutral, my.color[1], my.color)
      my.breaks = c(-0.0003, -0.0001, my.breaks)
  }
  #print(my.color)
  cob_heatmap = pheatmap(cob_matrix,
    color = my.color,
    breaks=my.breaks,
    legend = TRUE,
    cluster_rows = FALSE, cluster_cols = FALSE, 
#			 cellwidth=20, cellheight=20,
#                         display_numbers=TRUE,
      display_numbers=TRUE, #number_format = "%.5f",  # Change this to command line parameter
      number_format = myformat,

#      fontsize_number = 8,
#      fontsize_number = 32,

      fontsize = myfontsize,
	                 main=title, xlab="Gap length", ylab="Orientation", asp=1.0, las=1,
			 width=800, height=200)
#  title(xlab="Gap length", ylab="Orientation")
}

heatmap.heatplus <- function(cob_matrix, rows, cols)
{
#doLegend(1:9, g2r.colors(8), 4)
  cob_heatmap = regHeatmap(cob_matrix, legend=2, dendrogram = list(status = "no"))
  plot(cob_heatmap, main=title)

#  cob_heatmap = regHeatmap(cob_matrix, cluster_rows = FALSE, cluster_cols = FALSE, 
#	                   main=title, xlab="Gap length", ylab="Orientation", asp=1.0, las=1)
}

heatmap.two <- function(cob_matrix, rows, cols)
{
  number_of_colors=27
  cob_heatmap <- heatmap.2(cob_matrix, breaks=seq(-5,5, length=number_of_colors+1), 
	    key=FALSE, keysize=2.0, 
	    trace="none", symkey=FALSE, 
	    dendrogram="none", Rowv=NA, Colv=NA, 
	    main=title, xlab="Gap length", ylab="Orientation",
	    colsep=c(1:cols),
            rowsep=c(1:rows-1),
            sepcolor="white",
            sepwidth=c(0.005,0.005),
#	    lmat=rbind( c(3,2), 
#	    		c(1,4),
#			c(2,2)), 
#			lhei=c(0.25, 1, 1), lwid=c(2,0.5),
	    lmat=rbind( c(3,2), 
	    		c(1,4)),

	    lhei=c(0.5, 2),    # row height
	    lwid=c(2, 0.15),     # column width
#	    margin=c(6,5),
	    margin=c(6,5),               # margins for column and row names
	    col = greenred(number_of_colors), cexRow=1.0, scale="none", asp=1.0)


# In lmat layout
# 1 denotes color key
# 2 denotes column dendrogram
# 3 denotes row dendrogram
# 4 denotes image plot

}

#heatmap.image(cob_matrix, rows, cols)
heatmap.pretty(cob_matrix, rows, cols)
#heatmap.two(cob_matrix, rows, cols)
#heatmap.heatplus(cob_matrix, rows, cols)


