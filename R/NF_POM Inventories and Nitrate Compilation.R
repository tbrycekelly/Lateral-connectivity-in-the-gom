library(TheSource)
library(openxlsx)
library(geospt)



load('_rdata/NF_CTD Data (final).rdata')
nitrate = read.xlsx('Data/NF_Nutrient Data Final.xlsx', startRow = 2)
nitrate$Pressure = NA
nitrate$Temperature = NA
nitrate$Salinity = NA

for (i in 1:nrow(nitrate)) {
  if (is.na(nitrate$Bottle[i])) {
    l = which(bottle$Cruise == nitrate$Cruise[i] & bottle$Cast == nitrate$Cast[i])
    nitrate$Bottle[i] = bottle$Bottle[l[which.min((bottle$Depth[l] - nitrate$Depth[i])^2)]]
  }
  
  k = which(bottle$Cruise == nitrate$Cruise[i] & bottle$Cast == nitrate$Cast[i] & bottle$Bottle == nitrate$Bottle[i])
  nitrate$SigmaTheta[i] = bottle$Rho[k] - 1000
  nitrate$Pressure[i] = bottle$Pressure[k]
  nitrate$Temperature[i] = bottle$Temperature[k]
  nitrate$Salinity[i] = bottle$Salinity[k]
}

write.xlsx(nitrate, file = 'Data/temp.xlsx')



## Plots of Nitrate -Sigma Theta Relationship
par(mfrow=c(1,2))
l = nitrate$Depth > 50 & nitrate$Depth < 300
plot(nitrate$SigmaTheta[l], nitrate$Nitrate[l], col = make.pal(nitrate$Pressure[l], min = 30, max = 200), pch = 16)
plot(nitrate$SigmaTheta[!l], nitrate$Nitrate[!l], col = make.pal(nitrate$Pressure[!l], min = 30, max = 200), pch = 16)


l = which(l & !is.na(nitrate$Nitrate))
rel = smooth.spline(nitrate$SigmaTheta[l], nitrate$Nitrate[l], spar = 0.7)
lines(rel$x, rel$y, col = 'red', lwd = 2)

nitrate$PotTemp = calc.ptemp(S = nitrate$Salinity, T = nitrate$Temperature, P = nitrate$Pressure, P.ref = 0)
plot.TS(nitrate$Salinity, nitrate$PotTemp, xlim = c(32, 38), ylim = c(3, 30),
        col.pch = make.pal(nitrate$Pressure, min = 20, max = 200),
        freezing.line = F)

par(mfrow = c(2,3))
for (cycle in 1:5) {
  l = which(nitrate$Cycle == cycle)
  plot(nitrate$Nitrate[l], nitrate$Pressure[l], ylim = c(200,0), xlim = c(0,20), xaxs = 'i', yaxs = 'i', cex = 2)
  grid(); box()
  
  k = which(downcast$Cruise == nitrate$Cruise[l[1]] & downcast$Cast %in% nitrate$Cast[l])
  points(approx(rel$x, rel$y, rule = 2, xout = downcast$Rho[k]-1000)$y, downcast$Pressure[k], col = 'red', pch = 20)
  points(nitrate$Nitrate[l], nitrate$Pressure[l], cex = 2)
}

####### POM Inventories

load('_rdata/NF_CTD Data (final).rdata')
pom = read.xlsx('Data/NF_POM Data.xlsx', startRow = 2)
pom = pom[!is.na(pom$Cast),]
pom$mld = NA

inventory = unique(pom[,c('Cruise', 'Cast')])
inventory$mld = NA
inventory$pon = NA
inventory$nit = NA
inventory$ammon = NA
inventory$Cycle = NA
inventory$Date = make.time()

for (i in 1:nrow(inventory)) {
  l = which(downcast$Cruise == inventory$Cruise[i] & downcast$Cast == inventory$Cast[i])
  
  if (length(l) > 2) {
    inventory$mld[i] = calc.mld(downcast$Depth[l], downcast$Rho[l])
    
    k = which(inventory$Cruise[i] == pom$Cruise & inventory$Cast[i] == pom$Cast)
    profile = approx(x = pom$Depth[k], y = pom$N[k], rule = 2, xout = c(1:inventory$mld[i]))$y
    profile = approx(x = pom$Depth[k], y = pom$N[k], rule = 2, xout = c(1:60))$y
    inventory$pon[i] = mean(profile)
    inventory$Cycle[i] = pom$Cycle[k[1]]
    inventory$Date[i] = conv.time.excel(as.numeric(pom$Date[k[1]]))
  }
  
  k = which(nut$Cruise == inventory$Cruise[i] & nut$Cast == inventory$Cast[i])
  
  if (length(k) > 1) {
    profile = approx(x = nut$Depth[k], y = nut$Nitrate[k], rule = 2, xout = c(1:60))$y
    inventory$nit[i] = mean(profile)
    profile = approx(x = nut$Depth[k], y = nut$Ammonium[k], rule = 2, xout = c(1:60))$y
    inventory$ammon[i] = mean(profile)
  }
}

inventory = inventory[inventory$Cycle %in% c(1:5),]
inventory$Cycle[inventory$Cruise == 'NF18' & inventory$Cycle == '1'] = 4
inventory$Cycle[inventory$Cruise == 'NF18' & inventory$Cycle == '2'] = 5


par(mfrow = c(2,1))
plot(get.julian(inventory$Date), inventory$pon, col = cols$cycles[as.numeric(inventory$Cycle)], xlab = 'Julian Day', ylab = 'PON Inventory',
     pch = 16, ylim = c(0, 2), yaxs = 'i')
plot(get.julian(inventory$Date), inventory$nit, col = cols$cycles[as.numeric(inventory$Cycle)], xlab = 'Julian Day', ylab = 'PON Inventory',
     pch = 16, yaxs = 'i')

nutri = data.frame(Cast = unique(nut$Cast), Cycle = NA, nit = NA, ammon = NA, Cruise = NA)

for (i in 1:nrow(nutri)) {
  k = which(nut$Cast == nutri$Cast[i])
  
  if (length(k) > 1 & max(nut$Depth[k]) < 200) {
    nutri$Cruise[i] = nut$Cruise[k[1]]
    nutri$Cycle[i] = nut$Cycle[k[1]]
    
    profile = approx(x = nut$Depth[k], y = nut$Nitrate[k], rule = 2, xout = c(1:60))$y
    nutri$nit[i] = mean(profile)
    #profile = approx(x = nut$Depth[k], y = nut$Ammonium[k], rule = 2, xout = c(1:60))$y
    #inventory$ammon[i] = mean(profile)
  }
}

### Chl

chl = read.xlsx('Data/NF_Chl Samples.xlsx', startRow = 2)
chl$Datetime = conv.time.excel(chl$Datetime)

chl$ZChl = NA

for (time in unique(chl$Datetime)) {
  l = which(chl$Datetime == time)
  profile = approx(chl$Depth[l], chl$Chl[l], xout = c(1:max(chl$Depth[l])), rule = 2)
  profile = approx(chl$Depth[l], chl$Chl[l], xout = c(1:60), rule = 2)
  chl$ZChl[l] = sum(profile$y[profile$x <= max(chl$Depth[l])])
}

plot(get.julian(chl$Datetime), chl$ZChl, col = cols$cycles[chl$Cycle], pch = 16, ylim = c(0,30), yaxs = 'i',
     ylab = 'Chl (mg Chl m-3)', xlab = 'Julian Day')


R24 = 1
R4 = 2
N15 = 0.01

R0 = (R24 - R4) / (24-4)*-4 + R4
for (i in 1:5) {
  nit = 1 / 4 * N15 * (R4 - R0)
  
  R0 = nit * 4 / N15 + R4
  
  message('Nit = ', nit, '\t\tR0 = ', R0)
}


ss = read.xlsx('Data/NF_In Situ NO3 and CO3 2019.10.24.xlsx', sheet = 'HCO3 Summary', startRow = 2)
ss$Date = conv.time.excel(ss$Date)
inventory = unique(ss[,c('Cruise', 'Cast', 'Cycle')])
inventory$npp = NA
inventory$nit = NA
inventory$ammon = NA
inventory$Date = make.time()
inventory$Day = NA
for (i in 1:nrow(inventory)) {
  l = which(inventory$Cruise[i] == ss$Cruise & inventory$Cast[i] == ss$Cast)
  
  if (length(l) > 1) {
    inventory$Date[i] = ss$Date[l[1]]
    inventory$Day[i] = ss$Day[l[1]]
    
    profile = approx(x = ss$Depth[l], y = ss$Nitrate[l], rule = 2, xout = c(1:60))$y
    inventory$nit[i] = mean(profile)
    
    profile = approx(x = ss$Depth[l], y = ss$Ammonium[l], rule = 2, xout = c(1:60))$y
    inventory$ammon[i] = mean(profile)
    
    profile = approx(x = ss$Depth[l], y = ss$Prim.Prod[l], rule = 2, xout = c(1:60))$y
    inventory$npp[i] = mean(profile)
  }
}

pdf('_figures/Supplemental Figure 1 (part2).pdf')
{
  par(plt = c(0.2, 0.7, 0.7, 0.9))
  plot(inventory$Day, inventory$npp*60, pch = 16, ylim = c(0, 50),
       yaxs = 'i', xaxt = 'n', xlab = '', ylab = 'NPP\n(mmol C m-2)', cex = 1.5)
  points(inventory$Day, inventory$npp*60, col = cols$cycles[inventory$Cycle], pch = 16, cex = 1.2)
  mtext('c', adj = 0.05, line = -1)
  
  par(plt = c(0.2, 0.7, 0.5, 0.695), new = T)
  plot(inventory$Day, inventory$nit*60, pch = 16, ylim = c(0, 10),
       yaxs = 'i', xaxt = 'n', xlab = '', ylab = 'Nitrate\n(mmol N m-2)', cex = 1.5)
  points(inventory$Day, inventory$nit*60, col = cols$cycles[inventory$Cycle], pch = 16, cex = 1.2)
  mtext('d', adj = 0.05, line = -1)
  
  par(plt = c(0.2, 0.7, 0.3, 0.495), new = T)
  plot(inventory$Day, inventory$ammon*60, pch = 16, ylim = c(0, 10),
       yaxs = 'i', ylab = 'Ammonium\n(mmol N m-2)', xlab = 'Experimental Day', xaxt = 'n', cex = 1.5)
  points(inventory$Day, inventory$ammon*60, col = cols$cycles[inventory$Cycle], pch = 16, cex = 1.2)
  axis(1, at = c(1:4))
  mtext('e', adj = 0.05, line = -1)
  
  par(plt = c(0.7, 0.9, 0.3, 0.8), new = T)
  plot.new()
  legend('top', legend = c(1:5), col = cols$cycles, pch = 16, cex = 1.3)

}
dev.off()
