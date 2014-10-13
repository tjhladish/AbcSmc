source("pairs.panels.R")
d6 = read.table("dengue_predictive_prior-full_ts.06", header=T)

names(d6)[4:12] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
names(d6)[13] = "Autocorr"

png("pairs-full.png", width=1800, height=1340, res=150)
pairs.panels(d6[,4:13])#, gap=0.5)
dev.off()

d6p = d6

names(d6p)[c(6,7,9,10)] = c('log10(Daily intros)', 'Mosq/loc', 'Median', 'Std deviation')

png("pairs-part.png", width=1800, height=1340, res=200)
pairs.panels(d6p[,c(6,7,9,10)])#, gap=0.5)
dev.off()
