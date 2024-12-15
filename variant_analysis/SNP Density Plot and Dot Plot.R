####################################
# SNP Density Plot and SNP dot plot
# 
# Author: Jessica Luc, 2022
####################################

dprop = read.csv("run2propdensity.csv",header=T)

#SNP Density Plot
plot(dprop$POS,dprop$efk3b.i12, type = 'b', col = 1, xlab='Position',ylab='Fraction', 
     xlim=c(0,30000), ylim = c(0,1)) +
  lines(dprop$POS,dprop$efk3b.i24, type = 'b', col = 2) +
  legend(0,0.8, legend=c('Efk3B-hACE2-12', 'Efk3B-hACE2-24'), col = c(1,2), lty = 1:1, cex=0.8)
text(x = 23616, y = 0.5, labels = 23616, pos = 3)

#Dot Plot
plot(dprop$efk3b.i12, dprop$efk3b.i24, type = 'p', xlab = '12 hr post infection', ylab='24 hr post infection',
     main = 'Efk3B-hACE2', ylim = c(0,1), xlim = c(0,1))
lines(x=c(0,1), y=c(0,1), type = 'l')
#label all points
text(dprop$efk3b.i12, dprop$efk3b.i24, labels = dprop$POS)
#label POS 23616 point
text(x = 0.289, y= 0.66, labels = 23616, pos = 2)
text(x=0, y = 1, labels= "23616 (G => C)", pos = 4)
text(x=0, y=0.96, labels= "Arg => Pro Furin site mutation", pos = 4)

