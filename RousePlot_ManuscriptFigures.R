## This code will plot PrimingOct2010 data with model predictions
## Stephen Denton, Feb 2010
rm(list=ls())
graphics.off()
require("Hmisc")

create.files = F

load("PrimingObsData_1Beta.Rdata")
pred.acc.res = pred.acc

load("PrimingObsData_Final.Rdata")
print(pred.acc)
nSubjs = length(subjNums)  # 91
deg.fr = 90 #nSubjs-1
colnames(mean.acc) = c("Neither (NF)", "Target (TF)", "Foil (FF)")
colnames(mean.acc.pr) = c('Short', 'Long')
colnames(mean.acc.fl) = c('Short', 'Long')

# Generate barplot of overall data
plot.ylim = c(0.3,1.05)
leg.x = 12.25
leg.y = 1.075

errorbars = qt(.975,deg.fr) * sem
dotSize = 1.3 
dotlwd = 3
errBarWidth = 1.5
clwd = 1.5

fig.width = 6.5
fig.height = 6.5
if (create.files) {
    png(filename = "Data+Model_Fig.png", width=fig.width, height=fig.height,
        units = "in", pointsize = 12, res = 600)
} else { windows(fig.width,fig.height) }
layout(matrix(1:4,2,2,byrow=T))
par(cex.lab=1)
par(cex.axis=.8)
par(cex.main=1)
leg.cex = 0.7
leg.labs = c("Neither (NP)", "Target (TP)", "Foil (FP)")
par(adj=0.5)
par(tck=0.025)
par(mar=c(2.7,2.7,3,1), mgp=c(1.35,.2,0))
meanData = mean.acc
barHts = errorbars 
x = barplot((meanData[,,1,1]), beside=T, legend.text=T, xpd=F, 
    args.legend=list(x=leg.x, y=leg.y, title='Priming Condition', legend=leg.labs, cex=leg.cex, xpd=NA),
    ylim=plot.ylim,ylab='Accuracy', xlab='Flanking Condition' )
title(main='A1: Short Prime, Short Flanker', adj=0.0 ) 
errbar(x, meanData[,,1,1], meanData[,,1,1]+barHts[,,1,1], 
                           meanData[,,1,1]-barHts[,,1,1],add=T,pch=' ',lwd=errBarWidth,cap=0.04)
points(x, pred.acc[,,1,1], pch=20, cex = dotSize, lwd=dotlwd)
points(x, pred.acc.res[,,1,1], pch=4, cex = dotSize, lwd=clwd)
abline(.75,0,lty=2)

x = barplot((meanData[,,2,1]), beside=T, legend.text=T, xpd=F, 
    args.legend=list(x=leg.x, y=leg.y, title='Priming Condition', legend=leg.labs, cex=leg.cex, xpd=NA), 
    ylim=plot.ylim,ylab='Accuracy', xlab='Flanking Condition' )
title(main='A2: Long Prime, Short Flanker', adj=0.0 ) 
errbar(x, meanData[,,2,1], meanData[,,2,1]+barHts[,,2,1], 
                           meanData[,,2,1]-barHts[,,2,1],add=T,pch=' ',lwd=errBarWidth,cap=0.04)
points(x, pred.acc[,,2,1], pch=20, cex=dotSize, lwd=dotlwd)
points(x, pred.acc.res[,,2,1], pch=4, cex = dotSize, lwd=clwd)
abline(.75,0,lty=2)

x = barplot((meanData[,,1,2]), beside=T, legend.text=T, xpd=F, 
    args.legend=list(x=leg.x, y=leg.y, title='Priming Condition', legend=leg.labs, cex=leg.cex, xpd=NA),
    ylim=plot.ylim,ylab="Accuracy", xlab='Flanking Condition' )
title(main='B1: Short Prime, Long Flanker', adj=0.0 ) 
errbar(x, meanData[,,1,2], meanData[,,1,2]+barHts[,,1,2], 
                           meanData[,,1,2]-barHts[,,1,2],add=T,pch=' ',lwd=errBarWidth,cap=0.04)
points(x, pred.acc[,,1,2], pch=20, cex = dotSize, lwd=dotlwd)
points(x, pred.acc.res[,,1,2], pch=4, cex = dotSize, lwd=clwd)
abline(.75,0,lty=2)

x = barplot((meanData[,,2,2]), beside=T, legend.text=T, xpd=F, 
    args.legend=list(x=leg.x, y=leg.y, title='Priming Condition', legend=leg.labs, cex=leg.cex, xpd=NA),
    ylim=plot.ylim,ylab="Accuracy", xlab='Flanking Condition' )
title(main='B2: Long Prime, Long Flanker', adj=0.0 )
errbar(x, meanData[,,2,2], meanData[,,2,2]+barHts[,,2,2], 
                           meanData[,,2,2]-barHts[,,2,2],add=T,pch=' ',lwd=errBarWidth,cap=0.04)
points(x, pred.acc[,,2,2], pch=20, cex = dotSize, lwd=dotlwd)
points(x, pred.acc.res[,,2,2], pch=4, cex = dotSize, lwd=clwd)
abline(.75,0,lty=2)

if (create.files) { dev.off() }


fig.width = 6.5
fig.height = 6.5
if (create.files) {
    png(filename = "Data+Model_Mar_Fig.png", width=fig.width, height=fig.height,
        units = "in", pointsize = 12, res = 600)
} else { windows(fig.width,fig.height) }
layout(matrix(1:4,2,2,byrow=T))
par(cex.lab=1)
par(cex.axis=.8)
par(cex.main=1)
leg.cex = 0.7
par(adj=0.5)
par(tck=0.025)
par(mar=c(2.7,2.7,3,1), mgp=c(1.35,.2,0))

dotSize = 2
leg.x = 8.1
leg.y = 1.02

# Calculate priming effects and flanker effects separately
meanData = mean.acc.pr
barHts = qt(.975,deg.fr) * sem.pr
model.pred = apply(pred.acc, c(1,3,4), mean)
model.pred.res = apply(pred.acc.res, c(1,3,4), mean)
x = barplot(meanData[,,1], beside=T, legend.text=T, xpd=F,
        args.legend=list(x=leg.x, y=leg.y, horiz=T, title='Priming Condition', legend=leg.labs, cex=leg.cex, xpd=NA), 
        ylim=plot.ylim, ylab='Accuracy', xlab='Prime Duration' )
title(main='Short Flanker', adj=0, line=-.2, font.main=1)
text(-.4,1.15,'A: Marginal Priming Effects', xpd=NA, cex=1.3, font=2, adj=0)
errbar(x, meanData[,,1], meanData[,,1]+barHts[,,1], 
          meanData[,,1]-barHts[,,1], add=T, pch=' ',lwd=1.5,cap=0.04)
points(x, model.pred[,,1], pch=20, cex=dotSize, lwd=dotlwd)
points(x, model.pred.res[,,1], pch=4, cex=dotSize, lwd=clwd)
abline(.75,0,lty=2)
x = barplot(meanData[,,2], beside=T, legend.text=T, xpd=F,
         args.legend=list(x=leg.x, y=leg.y, horiz=T, title='Priming Condition', legend=leg.labs, cex=leg.cex, xpd=NA), 
        ylim=plot.ylim, ylab='Accuracy',, xlab='Prime Duration' )
title(main='Long Flanker', adj=0, line=-.2, font.main=1)
errbar(x, meanData[,,2], meanData[,,2]+barHts[,,2], 
          meanData[,,2]-barHts[,,2], add=T, pch=' ',lwd=1.5,cap=0.04)
points(x, model.pred[,,2], pch=20, cex=dotSize, lwd=dotlwd)
points(x, model.pred.res[,,2], pch=4, cex=dotSize, lwd=clwd)
abline(.75,0,lty=2)

leg.labs = c("Neither (NF)", "Target (TF)", "Foil (FF)")
meanData = mean.acc.fl
barHts = qt(.975,deg.fr) * sem.fl
model.pred = apply(pred.acc, c(2,3,4), mean)
model.pred.res = apply(pred.acc.res, c(2,3,4), mean)
x = barplot(meanData[,,1], beside=T, legend.text=T, xpd=F,
      args.legend=list(x=leg.x, y=leg.y, horiz=T, title='Flanker Condition', legend=leg.labs, cex=leg.cex, xpd=NA),
        ylim=plot.ylim, ylab='Accuracy',, xlab='Prime Duration' )
title(main='Short Flanker', adj=0, line=-.2, font.main=1)
text(-.4,1.15,'B: Marginal Flanker Effects', xpd=NA, cex=1.3, font=2, adj=0)
errbar(x, meanData[,,1], meanData[,,1]+barHts[,,1], 
          meanData[,,1]-barHts[,,1], add=T, pch=' ',lwd=1.5,cap=0.04)
points(x, model.pred[,,1], pch=20, cex=dotSize, lwd=dotlwd)
points(x, model.pred.res[,,1], pch=4, cex=dotSize, lwd=clwd)
abline(.75,0,lty=2)
x = barplot(meanData[,,2], beside=T, legend.text=T, xpd=F,
      args.legend=list(x=leg.x, y=leg.y, horiz=T, title='Flanker Condition', legend=leg.labs, cex=leg.cex, xpd=NA),
        ylim=plot.ylim, ylab='Accuracy',, xlab='Prime Duration' )
title(main='Long Flanker', adj=0, line=-.2, font.main=1)
errbar(x, meanData[,,2], meanData[,,2]+barHts[,,2], 
          meanData[,,2]-barHts[,,2], add=T, pch=' ',lwd=1.5,cap=0.04)
points(x, model.pred[,,2], pch=20, cex=dotSize, lwd=dotlwd)
points(x, model.pred.res[,,2], pch=4, cex=dotSize, lwd=clwd)
abline(.75,0,lty=2)
abline(1.2,0,lty=1, xpd=NA, lwd=2)

if (create.files) { dev.off() }
