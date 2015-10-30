require("RSQLite")
require(beanplot)
drv = dbDriver("SQLite")
db = dbConnect(drv, "refit.sqlite")
abc = dbGetQuery(db, 'select J.*, P.*, M.* from jobs J, parameters P, metrics M where J.serial = P.serial and J.serial = M.serial')
extra_serials = which(names(abc)== 'serial')[-1] 
abc = abc[,-c(extra_serials)]
abc = subset(abc, select=-c(startTime, duration, attempts, seed))
abc$post_bool = abc$posterior >= 0

par_cols = 6:11
met_cols = c(12:24)
obs_met = c( 121.943, 0, 2.68004, 33.9603, 147.642, 646.692, 176.722, 1.38442, 0.235294, 0.6, 0.129899, -5.06394, 0.116576 )

incomplete_sets = unique(abc$smcSet[abc$status!='D']) # 6
all_sets = unique(abc$smcSet)
complete_sets = setdiff(all_sets, incomplete_sets)
last_complete_set = max(complete_sets)

pdf('marginal_pars.pdf', width=8.5, height=11)
par(mfrow=c(6,1))
par(mar=c(2.1, 4.1, 1.1, 0.5))
for (col in par_cols) {
    colname = names(abc)[col];
    cat(paste(colname, '\n'))
    beanplot( abc[,colname] ~ abc$post_bool*abc$smcSet, what=c(1,1,1,0), col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1,1, grey(0.7))), main='', ylab=colname, side='both')
}
dev.off()

#ylims = list(c(0,800), c(0,5), c(0,200), c(0,350), c(0,600), c(10,30000), c(10e-1, 10e4), c(-2, 6), c(0, 1), c(0, 1), c(0, 1), c(-6, 4), c(-0.1, 0.1))

# 100 for these values indicates an inability to do a logistic fit.  This is rare, but could be avoided if we made the regression approach more robust.
abc$beta0[abc$beta0 == 100] = NA
abc$beta1[abc$beta1 == 100] = NA

# Plot metrics - tends to be trickier, because distributions can be weird (e.g, highly skewed or long-tailed)
pdf('marginal_mets.pdf', width=8.5, height=22)
par(mfrow=c(13,1))
par(mar=c(2.1, 4.1, 1.1, 0.5))
alpha = 0.025 # plot middle 95% of distributions
for (col in met_cols) {
    met_idx = col - max(par_cols)
    colname = names(abc)[col]

    # calculate reasonable plot limits
    obs_val = obs_met[met_idx]
    val_lims = quantile(abc[,colname], na.rm=T, probs=c(alpha, 1-alpha))
    val_min = min(val_lims[1], obs_val) 
    val_max = max(val_lims[2], obs_val) 

    # filter out NAs
    complete = complete.cases(abc[,colname]) & abc[,colname] > val_lims[1] & abc[,colname] < val_lims[2]
    cat(paste0(colname, ' ', val_lims[1], ', ', val_lims[2], '\n'))

    # plot the stuff
    beanplot( abc[complete, colname] ~ abc$post_bool[complete] * abc$smcSet[complete], what=c(0,1,1,0), 
              col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1, 1, grey(0.7))), 
              #main='', ylab=colname, side='both', ylim=unlist(ylims[met_idx]) )
              main='', ylab=colname, side='both', ylim=c(val_min, val_max), na.rm=T, log='' )
    abline(h=mean(abc[abc$smcSet==last_complete_set, colname], na.rm=T), lty=3)
    abline(h=obs_val, col=2, lwd=1)
}

dev.off()
