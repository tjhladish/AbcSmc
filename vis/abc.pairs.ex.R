source("pairs.panels.R")
d6 = read.table("dengue_predictive_prior-varEIP.00", header=T)

############## from script labeling empirical metrics
m = read.table("dengue.metrics", header=T)
proto = d6[1,]
proto[1,] <- NA
proto[1,names(m)] <- m
dm = rbind(d6, proto)
dm$type = factor(ifelse(is.na(dm$iteration), "met", "sim"))

colors <- c(sim = '#00000012', met = '#00FF00FF')
############## end



#names(dm)[4:12] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
names(dm)[4:13] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Beta', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
#names(dm)[13] = "Autocorr"
names(dm)[14] = "Autocorr"

png("pairs-a-varEIP00.png", width=1800, height=1340, res=150)
#pairs.panels(dm[,4:14], points.col=colors[as.character(dm$type)])#, gap=0.5)
pairs.panels(dm[,4:14], points.col=colors[as.character(dm$type)], box.col='black', box.lwd=0.5)#, gap=0.5)
dev.off()

d6p = dm

names(dm)[c(6,7,10,11)] = c('log10(Daily intros)', 'Mosq/loc', 'Median', 'Std deviation')

#png("pairs-b-varEIP00.png", width=1800, height=1340, res=200)
#pairs.panels(dm[,c(6,7,10,11)], points.col=colors[as.character(dm$type)], box.lwd=1)#, gap=0.5)
#dev.off()
