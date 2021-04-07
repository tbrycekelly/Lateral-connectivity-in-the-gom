library(TheSource)
library(openxlsx)

cols = list(obs = '#4A9586',
            obs2 = '#A5D3CA',
            cycles = c('#4A9586', '#8DC7BB', '#DCEDEA', '#D54FD5', '#E697E6'),
            satellite = '#25A0C5',
            satellite2 = '#8ED6EA',
            nemuro = '#D73E68',
            nemuro2 = '#F0B9C8')

data = read.xlsx('Data/NF_ST and Z Inventories.xlsx')
st = read.xlsx('Data/NF_SedTrap Data.xlsx', sheet = 'Flux', startRow = 2)
nitrate = read.xlsx('Data/NF_In Situ NO3 and CO3 2019.10.24.xlsx', sheet = 'Summary Data', startRow = 2)


{ ## Plot N flux profiles
  par(plt = c(0.2, 0.35, 0.2, 0.9), new = F)
  plot(NULL, NULL, xlim = c(0,1.6e3), ylim = c(210,0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)', xlab = 'Nitrogen Flux\n(umol N m-2 d-1)')
  grid(); box(); mtext('d', line = -1.5, side = 3, adj = 0.02, cex = 1.5)
  for (i in unique(st$Cycle)) {
    l = which(st$Cycle == i)
    lines(st$Norg[l]/15*1e3, st$Depth[l], lwd = 3, col = cols$cycles[i])
    points(st$Norg[l]/15*1e3, st$Depth[l], pch = 16)
    add.error.bars(x = st$Norg[l]/15*1e3, s.x = st$sNorg[l]/15, y = st$Depth[l], s.y = 0)
    text(x = st$Norg[min(l)]/15*1e3, st$Depth[min(l)], i, pos = 3, col = 'darkgrey')
  }
}


{ ## Plot N flux profiles
  par(plt = c(0.2, 0.4, 0.2, 0.9), new = F)
  plot(NULL, NULL, xlim = c(0,2e3), ylim = c(210,0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)', xlab = 'Nitrogen Flux (umol N m-2 d-1)')
  grid(); box(); mtext('A', line = -1.5, side = 3, adj = 0.02, cex = 1.5)
  for (i in unique(st$Cycle)) {
    l = which(st$Cycle == i)
    lines(st$Norg[l]/15*1e3, st$Depth[l], lwd = 2, col = cols$cycles[i])
    points(st$Norg[l]/15*1e3, st$Depth[l], pch = 16)
    add.error.bars(x = st$Norg[l]/15*1e3, s.x = st$sNorg[l]/15, y = st$Depth[l], s.y = 0)
    text(x = st$Norg[min(l)]/15*1e3, st$Depth[min(l)], i, pos = 3, col = 'darkgrey')
  }
  
  par(plt = c(0.403, 0.6, 0.2, 0.9), new = T)
  plot(NULL, NULL, xlim = c(0, 10), ylim = c(210, 0), yaxs = 'i', xaxs = 'i', ylab = '', yaxt = 'n', xlab = 'Nitrate Uptake (mmol N m-2 d-1)')
  grid(); box(); mtext('B', line = -1.5, side = 3, adj = 0.02, cex = 1.5)
  for (i in 1:5) {
    l = which(st$Cycle == i) ## St Data for pairing
    l = l[order(st$Depth[l])] ## order from top to bottom
    ll = which(nitrate$Cycle == i)
    nit = rep(NA, length(l))
    
    for (j in 1:length(l)) {
      nit[j] = sum(approx(x = nitrate$Depth[ll], y = nitrate$Rho[ll], rule = 2, xout = c(1:st$Depth[l[j]]))$y) * 1e-3 # mmol N m-2 d-1
    }
    lines(nit, st$Depth[l]); points(nit, st$Depth[l], pch = 20)
  }
  
  par(plt = c(0.603, 0.95, 0.2, 0.9), new = T)
  plot(NULL, NULL, xlim = c(0, 10), ylim = c(210, 0), yaxs = 'i', xaxs = 'i', ylab = '', yaxt = 'n', xlab = 'Nitrate Uptake/ST Export')
  grid(); box(); mtext('C', line = -1.5, side = 3, adj = 0.02, cex = 1.5)
  for (i in 1:5) {
    l = which(st$Cycle == i) ## St Data for pairing
    l = l[order(st$Depth[l])] ## order from top to bottom
    ll = which(nitrate$Cycle == i)
    nit = rep(NA, length(l))
    
    for (j in 1:length(l)) {
      nit[j] = sum(approx(x = nitrate$Depth[ll], y = nitrate$Rho[ll], rule = 2, xout = c(1:st$Depth[l[j]]))$y) * 1e-3 # mmol N m-2 d-1
    }
    lines(nit/st$Norg[l] * 15, st$Depth[l]); points(nit/st$Norg[l] * 15, st$Depth[l], pch = 20)
  }
  
}



st$CN = (st$Corg / 12) / (st$Norg / 14)
{ ## Plot C:N profiles
  par(plt = c(0.51, 0.65, 0.2, 0.9), new = T)
  plot(NULL, NULL, xlim = c(5,14), ylim = c(210,0), yaxs = 'i', xaxs = 'i', ylab = '', xlab = 'C:N Ratio', yaxt = 'n')
  grid(); box(); mtext('B', line = -1.5, side = 3, adj = 0.02, cex = 1.5)
  for (i in unique(st$Cycle)) {
    l = which(st$Cycle == i)
    lines(st$CN[l], st$Depth[l], lwd = 2)
    points(st$CN[l], st$Depth[l], pch = 20)
    text(x = st$CN[max(l)], st$Depth[max(l)], i, pos = 1, col = 'darkgrey')
  }
}

{ ## e-ratio
  par(plt = c(0.2, 0.45, 0.2, 0.9))
  plot(data$STC / data$ZNPP, data$Depth, ylim = c(210,0), pch = 15, yaxs = 'i', xlim = c(0, 0.65), xaxs = 'i',
       col = make.pal(data$Cycle, pal = kovesi.rainbow, min = 0, max = 6), ylab = 'Depth (m)', xlab = 'e-ratio')
  grid(); box()
  
  par(plt = c(0.45, 0.7, 0.2, 0.9), new = TRUE)
  plot(data$STC, data$Depth, ylim = c(210,0), pch = 15, yaxs = 'i', xlim = c(0, 140), xaxs = 'i', yaxt = 'n',
       col = make.pal(data$Cycle, pal = kovesi.rainbow, min = 0, max = 6), ylab = '', xlab = 'ST Export')
  grid(); box()
  
  par(plt = c(0.7, 0.95, 0.2, 0.9), new = TRUE)
  plot(data$ZNPP, data$Depth, ylim = c(210,0), pch = 15, yaxs = 'i', xlim = c(0, 450), xaxs = 'i', yaxt = 'n',
       col = make.pal(data$Cycle, pal = kovesi.rainbow, min = 0, max = 6), ylab = '', xlab = 'ZNPP')
  grid(); box()
}
dev.off()

#### ISOTOPES

plot(NULL, NULL, xlim = c(-2,2), ylim = c(220,0), yaxs = 'i', xaxs = 'i',
     ylab = 'Depth', xlab = '13C Deviation')

cycles = unique(data$Cycle)
col = make.qual.pal(cycles)
for (i in 1:length(cycles)) {
  l = which(data$Cycle == cycles[i])
  l = l[order(data$Depth[l])]
  
  lines(data$del13Cp[l], data$Depth[l] + i*2, lwd = 2, col = col[i])
  add.error.bars(data$del13Cp[l], data$del13Cp.sd[l], data$Depth[l] + i*2, 0, col = col[i])
}


plot(NULL, NULL, xlim = c(-3,3), ylim = c(220,0), yaxs = 'i', xaxs = 'i',
     ylab = 'Depth', xlab = '15N Deviation')

for (i in 1:length(cycles)) {
  l = which(data$Cycle == cycles[i])
  l = l[order(data$Depth[l])]
  
  lines(data$del15Np[l], data$Depth[l] + i*2, lwd = 2, col = col[i])
  add.error.bars(data$del15Np[l], data$del15Np.sd[l], data$Depth[l] + i*2, 0, col = col[i])
}


#### Euphotic Zone Only ####

data = data[data$Depth < 160,]

plot(data$del13C, data$del15N, col = make.pal(data$Depth, 255, 30, 150), pch = 16,
     xlim = c(-24, -21), xaxs = 'i', ylim = c(-1, 6), yaxs = 'i', xlab = expression(delta^13~C),
     ylab = expression(delta^15~N))
add.error.bars(data$del13C, data$del13C.sd, data$del15N, data$del15N.sd, col = 'grey')

abline(h = 3.2, lty = 2, col = 'dark blue')
abline(h = -0.5, lty = 2, col = 'dark grey')
abline(h = 2, lty = 2, col = 'grey')

mtext('Subsurface NO3', 4, cex = 0.7, at = 3.2, col = 'dark blue')
mtext('POM', 4, cex = 0.7, at = 2, col = 'grey')
mtext('N2 Fixation', 4, cex = 0.7, at = -0.5, col = 'dark grey')

text.default(data$del13C, data$del15N, labels = data$Cycle, adj = 2, pos = 3)

