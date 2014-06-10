#!/usr/bin/env Rscript

argv = commandArgs(trailingOnly = TRUE)

setwd('~/Desktop')

depths1 = read.table('XUC025_HS.bam.depth', header = TRUE, sep = "\t")[,3]
depths2 = read.table('XUC025_HS_WGA.bam.depth', header = TRUE, sep = "\t")[,3]
depths3 = read.table('XUA175_HS_1.bam.depth', header = TRUE, sep = "\t")[,3]
depths4 = read.table('XUA175_HS.bam.depth', header = TRUE, sep = "\t")[,3]



maxmean = max(c(mean(depths1), mean(depths2), mean(depths3), mean(depths4)))
correction_factor1 = maxmean/mean(depths1) #divide 1 by this
correction_factor2 = maxmean/mean(depths2) #divide 2 by this
correction_factor3 = maxmean/mean(depths3) #divide 3 by this
correction_factor4 = maxmean/mean(depths4) #divide 4 by this


#genome-wide depth plot
#jpeg('Genomewide_coverage.jpg')
d = depths1
sites_per_window = 1000 
num_windows = length(d)/sites_per_window
groups = c(rep(seq(1,num_windows), each = floor(sites_per_window)), rep(sites_per_window,length(d) - floor(sites_per_window)*floor(num_windows)))
means1 = rowsum(d, groups)/sites_per_window
means1 = means1*correction_factor1
means1 = sort(means1)

d = depths2
sites_per_window = 1000
num_windows = length(d)/sites_per_window
groups = c(rep(seq(1,num_windows), each = floor(sites_per_window)), rep(sites_per_window,length(d) - floor(sites_per_window)*floor(num_windows)))
means2 = rowsum(d, groups)/sites_per_window
means2 = means2*correction_factor2
means2 = sort(means2)

d = depths3
sites_per_window = 1000
num_windows = length(d)/sites_per_window
groups = c(rep(seq(1,num_windows), each = floor(sites_per_window)), rep(sites_per_window,length(d) - floor(sites_per_window)*floor(num_windows)))
means3 = rowsum(d, groups)/sites_per_window
means3 = means3*correction_factor3
means3 = sort(means3)

d = depths4
sites_per_window = 1000
num_windows = length(d)/sites_per_window
groups = c(rep(seq(1,num_windows), each = floor(sites_per_window)), rep(sites_per_window,length(d) - floor(sites_per_window)*floor(num_windows)))
means4 = rowsum(d, groups)/sites_per_window
means4 = means4*correction_factor4
means4 = sort(means4)


#draw
jpeg('WGA_eff.jpg', width = 800, height = 600)

plot(means1[means1<300], type = 'l', xlab = '1kb Windows, Sorted By Coverage Depth', ylab = 'Mean-Equalized Coverage Depth', main = "Effect of WGA on Windowed Sequencing Coverage Depth", cex.lab=1.5, cex.main = 2, cex.axis = 1.5, col = 'green', lwd = 3)
#points(which(means1<15), means1[means1<15], type = 'l', col = 'red')
abline(h = mean(means1), col = 1, lty=2)
abline(h = 15, col = 'red', lty=2)

points(means2[means2<300], type = 'l', col = 'blue', lwd = 3)
#points(which(means2<15), means2[means2<15], type = 'l', col = 'red')
abline(h = 15, col = 'red', lty=2)
legend('topleft', 'Strains', c('XUC025_HS', 'XUC025_HS_WGA'), col = c('green', 'blue'), cex = 1.5, lwd = 3)


text(1500,100, 'Scaled Mean', cex = 1.75)
text(700, 35, "Callable", cex = 1.75)
dev.off()



#plot(density(depths))
#print("depths 1:")
#print(summary(depths1))
#print("depths 2:")
#print(summary(depths2))

jpeg('density1_2_trimmed.jpg')
plot(density(depths1[depths1 < 500]), main = "Genomewide Per-nucleotide Coverage Depth", ylim = c(0, 0.015))
#dev.off()
#jpeg('density2_trimmed.jpg')
lines(density(depths2[depths2 < 500]), col = 'green')
legend('topright', 'Strains', c('XUC025_HS', 'XUC025_HS_WGA'), col = c('blue', 'green'), pch = c(1, 1))
dev.off()
