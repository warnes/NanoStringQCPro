
##' Create a preview of a PNG
##'
##' Generates a resized, vertically-cropped preview version of the input PNG.
##'
##' @param inputFile
##' Input PNG filename
##'
##' @param outputFile
##' Output PNG filename
##'
##' @param width
##' Width (in pixels) for the preview image
##'
##' @param cropHeight
##' Height (in pixels) for the preview image (if the rescaled input is larger
##' than this, it will be cropped)
##'
##' @param res
##' Output PNG resolution (passed to the `res' argument of png())
##'
##' @return
##' A resized, vertically-cropped preview version of the input PNG.
##'
##' @importFrom png readPNG
##'
previewPNG <- function(inputFile, outputFile, width, cropHeight, res=72)
{
  inputRaster <- readPNG(inputFile)

  rescale_factor <- width / ncol(inputRaster)
  height <- min( nrow(inputRaster)*rescale_factor, cropHeight )

  croppedInputRaster <- inputRaster[ 1:(height/rescale_factor), 1:(width/rescale_factor), 1:3 ]
  
  png(filename = outputFile,
      width = width,
      height = height,
      units = "px",
      res = res)

  op <- par(mar = rep(0,4))

  y_max <- nrow(croppedInputRaster)
  x_max <- ncol(croppedInputRaster)
  if (y_max < x_max) {
    x <- 1:x_max
    y <- c(1:y_max, rep(y_max, x_max - y_max))
  } else {
    x <- c(1:x_max, rep(x_max, y_max - x_max))
    y <- 1:y_max
  }  
  plot( y, x, type = "n", xlab = "", ylab = "", axes = FALSE )

  rasterImage( croppedInputRaster, 1, 1, y_max, x_max )

  par(op)
  
  dev.off()
}

