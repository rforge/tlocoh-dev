
isodata.plot <- function(isoc.info, plotme=c("isoarea","isoear")[1], figs.per.page=1, legend=TRUE,
                             title=NULL, title.show=TRUE, subtitle=TRUE, 
                             mar=c(3.3, 3.2, if (title.show) (if (subtitle) 3.2 else 2.3) else 0.5, 0.5), mgp=c(2, 0.7, 0), 
                             png.fn=NULL, png.dir=NULL, png.dir.make=TRUE, png.width=800, png.height=png.width, png.pointsize=12+(png.width-480)/80,
                             png.fn.pre=NULL, png.fn.suf=NULL, png.overwrite=TRUE,
                             panel.num=NULL, panel.num.inside.plot=!title.show, bg="white", legend.space=if (legend) 0.05 else 0, ...) {

    ## isoc.info is a data frame created by lhs.exp.isodata
    ## This is used when you can't have a single object containing all the hullsets due to memory limits
    ## Instead you use lxy.lhs.batch to create single hullset objects and save them to disk

    if (!plotme %in% c("isoarea","isoear")) stop("Unknown value for plotme")
    if (!is(isoc.info, "data.frame")) stop("Unknown object type")
    
    ## See if the output directory exists
    if (is.null(png.fn) && !is.null(png.dir) && !file.exists(png.dir)) {
        if (png.dir.make) {
            dir.made <- dir.create(png.dir)
            if (!dir.made) stop(paste("Couldn't make output folder ", png.dir, sep=""))            
        } else {
            stop(paste(png.dir, " doesn't exist and png.dir.make is False, so can not continue."))
        }
    }
    res <- NULL
    
    kar.vals.all <- sort(unique(isoc.info[["kar"]]))
    iso.levels.all <- sort(unique(isoc.info[["iso.level"]]))
    col.overlay <- rainbow(length(iso.levels.all), end=5/6)
    
    x.mat <- matrix(kar.vals.all, ncol=length(iso.levels.all), nrow=length(kar.vals.all))
    y.mat <- matrix(NA, ncol=length(iso.levels.all), nrow=length(kar.vals.all))
    
    isoc.info <- transform(isoc.info, kar=as.factor(isoc.info$kar), iso.level=as.factor(isoc.info$iso.level))
    iso2colidx <- match(levels(isoc.info[["iso.level"]]), as.character(iso.levels.all))
    paramval2rowidx <- match(levels(isoc.info[["kar"]]), as.character(kar.vals.all))
    
    id.str <- paste(unique(as.character(isoc.info[,"id"])), collapse=".", sep="")
    param.str <- paste(levels(isoc.info[["mode"]]), collapse=".", sep="")
    s.str <- paste(unique(as.character(isoc.info[,"s"])), collapse=".", sep="")
    sort.metric.str <- paste(unique(as.character(isoc.info[,"sort.metric"])), collapse=".", sep="")
    fn.series <- paste(".", param.str, ".vs.", plotme, sep="")
    
    ## Open a PNG device if needed
    if (!is.null(png.dir) || !is.null(png.fn)) {
        png.fn.use <- NULL
        if (is.null(png.fn)) {
            png.fn.use <- file.path(png.dir, paste(png.fn.pre, id.str, ".s", s.str, ".srt-", sort.metric.str, fn.series, png.fn.suf, ".png", sep=""))
        } else {
            png.fn.use <- png.fn
        }
        if (file.exists(png.fn.use) && !png.overwrite) stop(paste(png.fn.use, "exists"))
        png(filename=png.fn.use, height=png.height, width=png.width, pointsize=png.pointsize, bg=bg)
        res <- c(res, png.fn.use)
    }
    opar <- par(mfrow = n2mfrow(figs.per.page), mar=mar, mgp=mgp)
    
    
    ## Construct the plot title
    if (title.show) {
        if (is.null(title)) {
            #hs.sub <- if (subtitle) paste("\n", paste(id.str, collapse=".", sep=""), ", s=", paste(s.str, collapse=".", sep=""), sep="") else ""
            title.use <- paste(param.str, " vs. ", if (plotme=="isoarea") "isopleth area" else "edge:area ratio ", sep="")
        } else {
            title.use <- title
        }
    } else {
        title.use <- NULL
    }
    
    #print("ok populate y.mat");browser()
    for (i in 1:nrow(isoc.info)) {
       y.mat[ paramval2rowidx[as.numeric(isoc.info[i,"kar"])], iso2colidx[as.numeric(isoc.info[i, "iso.level"])]] <- if (plotme=="isoarea") isoc.info[i, "area"] else isoc.info[i, "edge.len"] / isoc.info[i, "area"]
    }
    
    rx <- range(x.mat[,1])
    xlim <- rx - c(legend.space * (rx[2]-rx[1]), 0)
    
    matplot(x.mat, y.mat, xlim=xlim, type="b", xlab=param.str, ylab="area", col=col.overlay, main=title.use, pch=20, lty=3, ...)
    if (legend) legend("topleft", as.character(iso.levels.all), col=col.overlay, lty=1, title="iso level", cex=0.8)
    
    ## Add panel.num
    if (!is.null(panel.num)) {
      if (panel.num.inside.plot) {
          text(x=par("usr")[1], y=par("usr")[4], labels=panel.num, cex=2, adj=c(-0.3,1.2), font=2)
      } else {
          mar.old <- par("mar")
          par(mar=c(0, 0.3, 0.2, 0))
          title(main=panel.num, adj=0, cex.main=2, line=-1.5)
          par(mar=mar.old)
      }
    }
    
    ## Close PNG device
    if (!is.null(png.dir) || !is.null(png.fn)) invisible(dev.off())
        
    
    return(res)

}
