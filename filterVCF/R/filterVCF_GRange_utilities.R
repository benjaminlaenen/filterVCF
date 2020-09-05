#' Plot a Grange object
#' 
#' Simple representation of GRange interval using rectangle. If values are
#' associated with each range multiple track can be visualized at the same time.
#' 
#' @param x A GRange
#' @param xlim=x region to plot, default spans the min and max of the GRange
#' @param main=deparse(substitute(x)) Plot title
#' @param col="black" single value or list of color the length of the GRanges
#' @param colorRange = FALSE Should the value associated to each GRange be represented by a gradient
#' @param colorRamp = c("cadetblue1","brown4") Start and end of the color gradient
#' @param exclude.column = NULL Index of the column to exclude from plot
#' @param sep=0.5 Not in use
#' @param ... paremeter pass to rect() except col and border
#' @return NULL
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~plot
#' @family plot
#' @examples
#' 
#' 
#' @export
plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5, colorRange = FALSE, colorRamp = c("cadetblue1", "brown4"), exclude.column = NULL, ...){
	if(!is.null(exclude.column)){
		values(x) <- values(x)[,-exclude.column]
	}
	height <- ncol(values(x))
   if(height == 0){
      height  <- 1
      values(x) <- 1
   }
   if (is(xlim, "GenomicRanges")){
     xlim <- c(min(start(xlim)), max(end(xlim)))
   }

   plot.new()
   # bins <- disjointBins(IRanges(start(x), end(x) + 1))
   #plot.window(xlim, c(0, max(bins)*(height + sep)))
   #ybottom <- bins * (sep + height) - height
   plot.window(xlim, c(0, (height)))
   ybottom <- 0
   
   if(colorRange){
   	COL = cbind(sort(unique(unlist(values(x)))), c("white", colorRampPalette(colorRamp, bias = 1, interpolate = "linear", alpha = 1)(length(unique(unlist(values(x))))-1)))
   	Mat_x <- as.matrix(values(x))

   	for(i in seq_len(nrow(COL))){
   		Mat_x[Mat_x == as.integer(COL[i,1])] <- COL[i,2]
   	}

   	## Matrix to list
   	col <- split(Mat_x, rep(1:ncol(Mat_x), each = nrow(Mat_x)))
		LegendRamp <- legend(xlim[1], 0, legend = COL[,1], pch =20, col =COL[,2], bty = "n",cex =0.3, horiz = TRUE)

   }

   if(!is.list(col)) col <- as.list(rainbow(ncol(values(x))))

	for(i in seq_len(ncol(values(x))) ){
		rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + i, col=col[[i]], border = NA)
		ybottom = ybottom + i
	}
   title(main)
 axis(1)
 axis(2,at = 0.5:height, labels = colnames(values(x)), las = 2, cex = 0.5)
}
