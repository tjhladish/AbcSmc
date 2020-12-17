# pairs2() is almost the same as the default pairs plot, but axes are placed on the bottom and the left, instead of alternating.
# This is much more readable if only the lower triangle is being used to plot data
pairs2 <- 
function (x, labels, panel = points, ..., 
          lower.panel = panel, 
          upper.panel = panel, 
          diag.panel = NULL, 
          text.panel = textPanel, 
          label.pos = 0.5 + has.diag/3, 
          cex.labels = NULL, 
          font.labels = 1, 
          row1attop = TRUE, 
          gap = 1,
          box.col = 'black',
          box.lty = 'solid',
          box.lwd = 1) 
{
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) {
        text(x, y, txt, cex = cex, font = font)
    }
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, oma, ...) {
        if (side%%2 == 1) 
            Axis(x, side = side, xpd = NA, ...)
        else Axis(y, side = side, xpd = NA, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) { plot(...) }
    localLowerPanel <- function(..., main, oma, font.main, cex.main) { lower.panel(...) }
    localUpperPanel <- function(..., main, oma, font.main, cex.main) { upper.panel(...) }
    localDiagPanel <- function(..., main, oma, font.main, cex.main) { diag.panel(...) }
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]])) {
                x[[i]] <- as.numeric(x[[i]])
            }
            if (!is.numeric(unclass(x[[i]]))) {
                stop("non-numeric argument to 'pairs'")
            }
        }
    } else if (!is.numeric(x)) {
        stop("non-numeric argument to 'pairs'")
    }
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) {
        lower.panel <- match.fun(lower.panel)
    }
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) {
        upper.panel <- match.fun(upper.panel)
    }
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) {
        diag.panel <- match.fun(diag.panel)
    }
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) {
        stop("only one column in the argument to 'pairs'")
    }
    has.labs <- TRUE
    if (missing(labels)) {
        labels <- colnames(x)
        if (is.null(labels)) labels <- paste("var", 1L:nc)
    } else if (is.null(labels)) {
        has.labs <- FALSE
    }
    oma <- if ("oma" %in% nmdots) dots$oma else NULL
    main <- if ("main" %in% nmdots) dots$main else NULL
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main)) oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    for (i in if (row1attop) 1L:nc else nc:1L) {
        for (j in 1L:nc) {
            localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, type = "n", ...)
            if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
                #box(col=box.col, lty=box.lty, lwd=box.lwd, ...)
# edited here...
#           if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
#           localAxis(1 + 2 * row1attop, x[, j], x[, i], 
#                       ...)
# draw x-axis
                if (i == nc & j != nc) localAxis(1, x[, j], x[, i], las=2, ...)
# draw y-axis
                if (j == 1 & i != 1) localAxis(2, x[, j], x[, i], las=2, ...)
#           if (j == nc && (i%%2 || !has.upper || !has.lower)) 
#             localAxis(4, x[, j], x[, i], ...)
                mfg <- par("mfg")
                if (i == j) {
                    if (has.diag) localDiagPanel(x, i, ...)
                    #if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
                    if (has.labs) {
                        par(usr = c(0, 1, 0, 1))
                        if (is.null(cex.labels)) {
                            l.wid <- strwidth(labels, "user")
                            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                        }
                        text.panel(0.5, label.pos, labels[i], cex = cex.labels, font = font.labels)
                    }
                } else if (i < j) {
                    localLowerPanel(x, j, i, ...)
                } else {
                    localUpperPanel(x, j, i, ...)
                    #localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
                }
                if (any(par("mfg") != mfg)) stop("the 'panel' function made a new plot")
                box(col=box.col, lty=box.lty, lwd=box.lwd, ...)
            } else {
                par(new = FALSE)
            }
        }
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) dots$font.main else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}

"pairs.panels" <-
function (x,
        sims,
        npar,
        nmet,
        smooth = TRUE,
        density=TRUE,
        digits = 2,
        method="pearson",
        cor=TRUE,
        points.col="#00000033",
        points.pch=20,
        cor.cex=2,
        points.cex=2,
        line_wt=1,
        ...)   #combines a splom, histograms, and correlations
{

#first define all the "little functions"

    "colfunc" <- colorRampPalette(c("#ff0000", "#dddddd", "#0000ff"), space='Lab')
    #"colfunc" <- colorRampPalette(c("#dd0000", "#dddddd", "#0000dd"), space='Lab')
    colscale = colfunc(21) # covers seq(-1, 1, 0.1)
    
#    "panel.hist.density" <-
#        function(x,...) {
#            usr <- par("usr"); on.exit(par(usr))
#            par(usr = c(usr[1:2], 0, 1.5) )
#            h <- hist(x, plot = FALSE)
#            breaks <- h$breaks; nB <- length(breaks)
#            y <- h$counts; y <- y/max(y)
#            rect(breaks[-nB], 0, breaks[-1], y,col=hist.col)
#            tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
#            if(class(tryd) != "try-error") {
#                d$y <- d$y/max(d$y)
#                lines(d)
#            }
#        }

    "panel.hist" <-
        function(d, col1, ...)
        {
            x = as.vector(d[sims==T,col1])
            hist.col = if (col1 <= npar) '#00ff0088' else '#ffa50088';

            usr <- par("usr"); on.exit(par(usr))
            par(usr = c(usr[1:2], 0, 1.5) )
            #h <- hist(x, add=T, col=hist.col)
            inset = (usr[2] - usr[1])*0.01
            h <- hist(x, breaks=seq(usr[1]+inset, usr[2]-inset, length.out=11), plot = FALSE)
            breaks <- h$breaks; nB <- length(breaks)
            y <- h$counts; y <- y/max(y)
            rect(breaks[-nB], 0, breaks[-1], y, col=hist.col, border=NA)
        }

    "panel.smoother.noellipse" <- 
        function (d, col1, col2, pch = par("pch"), col.smooth = "red", span = 2/3, iter = 3, ...) 
        {
            x = as.vector(d[sims==T,col1])
            y = as.vector(d[sims==T,col2])
            r = cor(x, y,use="pairwise",method=method)
            points(x, y, col=points.col, pch = points.pch, cex=points.cex, ...)

            x1 = d[sims==F, col1]
            y1 = d[sims==F, col2]

            # version 3
            color1 = if (col1 <= npar) 'green' else 'orange';
            color2 = if (col2 <= npar) 'green' else 'orange';

            abline(v=x1, col=color1, lty=3, lwd=line_wt)
            abline(h=y1, col=color2, lty=3, lwd=line_wt)
            # version 2
            # color1 = if (col1 <= npar) '#AB82FF' else 'orange';
            # color2 = if (col2 <= npar) '#AB82FF' else 'orange';
            # points(x1, y1, col=color1, pch = 20, cex=points.cex+0.5, ...)
            # if (color1 != color2) {
            #     points(x1, y1, col=color2, pch = 20, cex=(points.cex-1), ...)
            # }

            # version 1
            # points(x1, y1, col=color1, pch = '|', cex=points.cex, ...)
            # points(x1, y1, col=color2, pch = '—', cex=points.cex, ...)

            ok <- is.finite(x) & is.finite(y)
            if (any(ok)) {
                r  <- cor(x, y, use="pairwise", method=method)
                col_idx = round(r*10) + 11
                lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = colscale[col_idx], lwd = line_wt, ...)
                #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
            }
        }

    "panel.cor" <-
        function(d, col1, col2, digits=2, prefix="", ...)
        {
            x = as.vector(d[sims==T,col1])
            y = as.vector(d[sims==T,col2])
            usr <- par("usr"); on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            r  <- cor(x, y,use="pairwise",method=method)
            txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
            txt <- paste(prefix, txt, sep="")
            #if(missing(cor.cex)) {
            #    cor.cex <- min(0.8/strwidth('-0.00'))
            #}
            col_idx = round(r*10) + 11
            text(0.5, 0.5, txt, cex=cor.cex, col=colscale[col_idx], font=2)
        }

#Beginning of the main function
#the organization gets around the problem of passing parameters to pairs by calling different routines
#this is truly clunky, but I don't know how else to do it
#It is possible that the trick I use for rug would work for other options.
#rug is a variable that is local to the entire function but is used only in the hist and hist.density functions.
#having done the same trick for passing method, it is now time to change all of this to pass some of these other options

    old.par <- par(no.readonly = TRUE) # save default, for resetting...
    on.exit(par(old.par))     #and when we quit the function, restore to original values

   # pairs2(x, diag.panel = panel.hist, lower.panel = panel.smoother.noellipse, ...)
    if (cor) {
        upper.panel = panel.cor
    } else {
        upper.panel = NULL
    }
pairs2(x, diag.panel = panel.hist, upper.panel = upper.panel, lower.panel = panel.smoother.noellipse, ...)
#pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.smoother.noellipse, ...)
}
