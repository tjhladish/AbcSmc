source("pairs.panels.R")
d6 = read.table("dengue_predictive_prior-varEIP.00", header=T)

require('rjson')
abc=fromJSON(file="abc_config.json"); 
nmet=length(abc$metrics); 
npar=length(abc$parameters); 
#factors=c(rep(1,nmet), rep(2,npar)); factors %o% factors
abc_metrics_values=rep(0,length(abc$metrics)); for(i in 1:length(abc$metrics)) { abc_metrics_values[i]=abc$metrics[[i]]$value; } 
abc_metrics_names=rep(0,length(abc$metrics)); for(i in 1:length(abc$metrics)) { abc_metrics_names[i]=abc$metrics[[i]]$name; } 
# we want short names if they're defined
for(i in 1:length(abc$metrics)) { if('short_name' %in% names(abc$metrics[[i]])) abc_metrics_names[i]=abc$metrics[[i]]$short_name; } 


############## from script labeling empirical metrics
#m = read.table("dengue.metrics", header=T)
proto = d6[1,]
proto[1,] <- NA
proto[1,abc_metrics_names] <- abc_metrics_values
dm = rbind(d6, proto)
dm$sim = factor(ifelse(is.na(dm$iteration), F, T))
dm[dm$sim==F,4:(3+npar)] = colMeans(dm[dm$sim==T,4:(3+npar)])

colors <- c(sim = '#00000012', par = 'purple', met = 'orange')
chars <- c(sim = 20, par = '|', met = '—')
############## end



#names(dm)[4:12] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
names(dm)[4:13] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Beta', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
#names(dm)[13] = "Autocorr"
names(dm)[14] = "Autocorr"

png("pairs-a-varEIP00.png", width=1800, height=1340, res=150)
#pairs.panels(dm[,4:14], points.col=colors[as.character(dm$sim)])#, gap=0.5)
pairs.panels(dm[,4:14], dm[,15], npar, nmet, points.col='#00000012', box.col='black', box.lwd=0.5)#, gap=0.5)
dev.off()

d6p = dm

names(dm)[c(6,7,10,11)] = c('log10(Daily intros)', 'Mosq/loc', 'Median', 'Std deviation')

png("pairs-b-varEIP00.png", width=1800, height=1340, res=200)
pairs.panels(dm[,c(6,7,10,11)], dm[,15], npar=2, nmet=2, points.col='#00000012', box.lwd=1)#, gap=0.5)
dev.off()
