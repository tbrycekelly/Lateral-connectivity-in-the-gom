library(TheSource)
library(openxlsx)

cols = list(obs = '#4A9586',
            obs2 = '#A5D3CA',
            cycles = c('#4A9586', '#8DC7BB', '#DCEDEA', '#D54FD5', '#E697E6'),
            satellite = '#25A0C5',
            satellite2 = '#8ED6EA',
            nemuro = '#D73E68',
            nemuro2 = '#F0B9C8')


## Load data
load('_rdata/NF_Kz Estimates 2020.04.27.rdata')
load('_rdata/NF_CTD Data (final).rdata')
#nut = read.xlsx('Data/NF_In Situ NO3 and CO3 2019.10.24.xlsx', startRow = 2, sheet = 1)
#nut$Nitrate[nut$Nitrate < 0.01] = 0.01
nut = read.xlsx('data/NF_Nutrient Data Final.xlsx', startRow = 2, sheet = 2)
load('_rdata/NF_Kz Estimates 2020.04.27.rdata')
st = read.xlsx('Data/NF_SedTrap Data.xlsx', sheet = 'Flux', startRow = 2)

nitrate = nut
## Create interpolated profiles of nitrate
interp.nitrate = data.frame(Cruise = NA, Cycle = NA, Depth = NA, Sigma = NA, Nitrate = NA)

for (cycle in c(1:5)) {
  l = which(cycle == nitrate$Cycle & !is.na(log10(nitrate$Nitrate)) & is.finite(log10(nitrate$Nitrate)))
  d = which(nitrate$Cast[l[1]] == downcast$Cast & nitrate$Cruise[l[1]] == downcast$Cruise)
  
  if (length(l) > 1) {
    nitrate.new = rep(mean(nitrate$Nitrate[l]), length(d))
    
    for (i in 1:length(d)) {
      dist = 0*abs(downcast$Rho[d[i]] - 1000 - nitrate$SigmaTheta[l])^3 + 5^3 + abs(downcast$Depth[d[i]] - nitrate$Depth[l])^3
      w = 1 / dist
      nitrate.new[i] = sum((nitrate$Nitrate[l]) * w) / sum(w)
    }
    
    nitrate.new = approx(nitrate$Depth[l], nitrate$Nitrate[l], xout = downcast$Depth[d], rule = 2)$y
    
    interp.nitrate = rbind(interp.nitrate,
                           data.frame(Cruise = nitrate$Cruise[l[1]],
                                      Cycle = cycle,
                                      Depth = downcast$Depth[d],
                                      Sigma = downcast$Rho[d]-1000,
                                      Nitrate = nitrate.new))
  }
}
interp.nitrate = interp.nitrate[-1,]


#### Supplemental Figure 1.
#
if (F) { # Three panel figure with Nitrate, Kz and flux estimates
  for (cycle in 1:5) {
    par(mfrow=c(1,3), plt = c(0.3, 1, 0.2, 0.8))
    l = which(interp.nitrate$Cycle == cycle)
    l.obs = which(nut$Cycle == cycle)
    depth = seq(1, 200, by = 0.1)
    
    plot(NULL, NULL, xlim = c(0, 20), ylim = c(200, 0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)', xlab = 'Nitrate (uM)')
    grid(); box()
    lines((interp.nitrate$Nitrate[l]), interp.nitrate$Depth[l], lwd = 2)
    points((nut$Nitrate[l.obs]), nut$Depth[l.obs], pch = 16, cex = 1.5, col = 'black')
    points((nut$Nitrate[l.obs]), nut$Depth[l.obs], pch = 16, col = cols$cycles[cycle])
    
    
    par(plt = c(0.01, 1, 0.2, 0.8))
    ## Plot Physics
    plot(NULL, NULL, xlim = c(-8, -4), ylim = c(200,0), yaxs = 'i', xaxs = 'i', xaxt = 'n', ylab = '', xlab = 'Kz (m2 s-1)', yaxt = 'n', main = paste0('Cycle ', cycle))
    add.log.axis(1, grid.major = T); grid(nx = NA, ny = 4); box()
    ll = which(kz$cycle == cycle)
    lines(log10(kz$kz[ll]), kz$depth[ll], col = cols$cycles[cycle], lwd = 2)
    
    par(plt = c(0.01, 0.9, 0.2, 0.8))
    plot(NULL, NULL, xlim = c(-2, 3), ylim = c(200,0), yaxs = 'i', xaxs = 'i', ylab = '', xlab = 'Nitrate Flux (umol N m-2 d-1)', yaxt = 'n', xaxt = 'n')
    grid(); abline(v = 0, col = 'grey'); box()
    J = 86400 * diff(interp.nitrate$Nitrate[l]) / diff(interp.nitrate$Depth[l]) * approx(kz$depth[ll], kz$kz[ll], rule = 2, xout = interp.nitrate$Depth[l[-1]])$y
    lines(log10(J)+3, interp.nitrate$Depth[l[-1]], col = '#008844', lwd = 2)
    add.log.axis(1)
    
  }
}

## Supplemental Figure 2
#pdf('_figures/Nitrate and Kz Profiles.pdf')
{ # Three panel figure with Nitrate, Kz and flux estimates
  par(mfrow=c(1,3), plt = c(0.3, 1, 0.2, 0.8))
  plot(NULL, NULL, xlim = c(0, 20), ylim = c(200, 0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)', xlab = 'Nitrate (uM)')
  grid(); box()
  
  for (cycle in 1:5) {
    l = which(interp.nitrate$Cycle == cycle)
    l.obs = which(nut$Cycle == cycle)
    depth = seq(1, 200, by = 0.1)
    
    lines((interp.nitrate$Nitrate[l]), interp.nitrate$Depth[l], lwd = 2, col = cols$cycles[cycle])
    points((nut$Nitrate[l.obs]), nut$Depth[l.obs], pch = 16, cex = 1.5, col = 'black')
    points((nut$Nitrate[l.obs]), nut$Depth[l.obs], pch = 16, col = cols$cycles[cycle])
  }
  legend('topright', legend = c(1:5), col = cols$cycles, lwd = 2, title = 'Cycle')
    
  par(plt = c(0.01, 1, 0.2, 0.8))
  ## Plot Physics
  plot(NULL, NULL, xlim = c(-8, -4), ylim = c(200,0), yaxs = 'i', xaxs = 'i', xaxt = 'n', ylab = '', xlab = 'Kz (m2 s-1)', yaxt = 'n')
  add.log.axis(1, grid.major = T); grid(nx = NA, ny = 4); box()
  for (cycle in 1:5) {
    ll = which(kz$cycle == cycle)
    lines(log10(kz$kz[ll]), kz$depth[ll], col = cols$cycles[cycle], lwd = 2)
  }
  
  par(plt = c(0.01, 0.9, 0.2, 0.8))
  plot(NULL, NULL, xlim = c(-2, 3), ylim = c(200,0), yaxs = 'i', xaxs = 'i', ylab = '', xlab = 'Nitrate Flux (umol N m-2 d-1)', yaxt = 'n', xaxt = 'n')
  grid(); abline(v = 0, col = 'grey'); box()
  for (cycle in 1:5) {
    l = which(interp.nitrate$Cycle == cycle)
    l.obs = which(nut$Cycle == cycle)
    ll = which(kz$cycle == cycle)
    J = 86400 * diff(interp.nitrate$Nitrate[l]) / diff(interp.nitrate$Depth[l]) * approx(kz$depth[ll], kz$kz[ll], rule = 2, xout = interp.nitrate$Depth[l[-1]])$y
    lines(log10(J)+3, interp.nitrate$Depth[l[-1]], col = cols$cycles[cycle], lwd = 2)
  }
  add.log.axis(1)
}
#dev.off()


budget = read.xlsx('Data/NF_Nitrogen Budget.xlsx', startRow = 2)

{ # Calculate nitrate flux for budgets
  for (i in 1:nrow(budget)) {
    l = which(nut$Cycle == budget$Cycle[i] & !is.na(nut$Nitrate))
    depth = seq(1, 200, by = 0.1)
    
    int = gridIDW(depth, rep(10, length(depth)), nut$Depth[l], rep(0, length(l)), nut$Nitrate[l], p  = 4, 1, 1, 1)
    int[int < 0.02] = 0.02
    ll = which(kz$cycle == budget$Cycle[i])
    J = 86400 * diff(int) * approx(kz$depth[ll], kz$kz[ll], rule = 2, xout = depth[-1])$y ## mmol m-4 m2 d-1
    
    k = which(abs(depth-budget$Depth[i]) < 5)
    budget$J.Nitrate[i] = 1 / mean(1 / J[k]) * 1e3  # umol N m-2 d-1
  }
}

st$J = NA
{ # TCalculate nitrate flux for ST data
  for (i in 1:nrow(st)) {
    l = which(nut$Cycle == st$Cycle[i] & !is.na(nut$Nitrate))
    depth = seq(1, 300, by = 0.1)
    
    int = gridIDW(depth, rep(10, length(depth)), nut$Depth[l], rep(0, length(l)), nut$Nitrate[l], p  = 4, 1, 1, 1)
    int[int < 0.02] = 0.02
    ll = which(kz$cycle == st$Cycle[i])
    J = 86400 * diff(int) * approx(kz$depth[ll], kz$kz[ll], rule = 2, xout = depth[-1])$y
    
    k = which(abs(depth-st$Depth[i]) < 5)
    st$J[i] = 1 / mean(1 / J[k]) * 1e3
  }
}

budget$J.Nitrate[budget$J.Nitrate < 0.1] = 0.1

#mort = read.xlsx('Data/NF_Mesozooplankton Data from Landry.xlsx', sheet = 'Mortality Summary', startRow = 2)

{
  par(mfrow = c(1,1), plt = c(0.3, 0.6, 0.2, 0.9))
  
  plot(NULL, NULL, ylim = c(15.5, 0.5), xlim = c(-1, 3.5), yaxt = 'n', yaxs = 'i', xlab = 'Flux (umol N m-2 d-1)', xaxt = 'n', ylab = '')
  #rect(-2, 0,4,20, col = 'darkgrey')
  axis(2, at = c(3, 8, 13), labels = c('Sed Trap\nExport', 'Zooplankton\nExcretion', 'Vertical\nNitate Supply'), tick = F, line = 1.5, las = 1)
  axis(2, at = c(1:15), label  = c(1:5, 1:5, 1:5), las = 1)
  #add.log.axis(1, grid.major = T)
  abline(v = c(-1:5), lty = 3, col = 'grey')
  axis(1, c(0:5), 10^c(0:5))
  axis(1, -1, '<0.1')
  abline(h = 5.5, lty = 2)
  abline(h = 10.5, lty = 2)
  
  points(log10(budget$SML.Flux), c(1:5), col = cols$cycles, pch = 22, cex = 1.5)
  #points(log10(budget$SML.Flux), c(1:5), pch = 1, cex = 1.5)
  add.error.bars.logx(budget$SML.Flux, budget$sSML.Flux, budget$Cycle, 0, 10)
  
  points(log10(budget$ST.N), c(1:5), col = cols$cycles, cex = 1.5, pch = 16)
  add.error.bars.logx(budget$ST.N, budget$sST.N, budget$Cycle, 0, 10)
  
  points(log10(budget$Ex), c(1:5)+5, col = cols$cycles, cex = 1.5, pch = 16)
  add.error.bars.logx(budget$Ex, budget$sEx, budget$Cycle+5, 0, 10)
  
  #for (i in 1:5) {
  #  lines(x = log10(mort$Mortality[mort$Cycle == i]), rep(i,2)+5, col = cols$obs2)
  #}
  #points(log10(mort$Mortality), c(1:5, 1:5)+5, col= cols$obs2, cex = 1, pch = 20)
  
  points(log10(budget$J.Nitrate), c(1:5)+10, col = cols$cycles, cex = 1.5, pch = c(16))
  points(log10(budget$J.Nitrate[1:3]), c(1:3)+10, col = 'white', cex = 1, pch = 16)
  add.error.bars.logx(budget$J.Nitrate, budget$J.Nitrate/2, budget$Cycle+10, 0, 10)
  #legend('bottomright', legend = c(1:5), col = col[1:5], title = 'Cycle', pch = 1, cex = 1.5)
  
  par(new = T, plt = c(0.1, 0.9, 0.91, 0.96), mfrow = c(1,1))
  plot.new()
  legend('top', legend = c(1:5, ''), col = c(cols$cycles, 'white'), pch = 16, cex = 1.2, horiz = T, bg = NA, box.col = 'white')
  #mtext('Cycle', adj = 0.24, line = -1.6)
}



{
  pdf('_figures/NF_Figure 2.pdf')
  par(mfrow = c(1,1), plt = c(0.3, 0.6, 0.2, 0.9))
  
  plot(NULL, NULL, ylim = c(20.5, 0.5), xlim = c(-1, 3.5), yaxt = 'n', yaxs = 'i', xlab = 'Flux (umol N m-2 d-1)', xaxt = 'n', ylab = '')
  #rect(-2, 0,4,20, col = 'darkgrey')
  axis(2, at = c(3, 8, 13, 18), labels = c('Sed Trap\nExport', 'Zooplankton\nExcretion', 'Vertical\nNitate Supply', 'N2 Fixation'), tick = F, line = 1.5, las = 1)
  #axis(2, at = c(1:20), label  = c(1:5, 1:5, 1:5, 1:5), las = 1)
  axis(2, at = c(0.5,5.5, 10.5, 15.5,20.5), labels = NA)
  #add.log.axis(1, grid.major = T)
  abline(v = c(-1:5), lty = 3, col = 'grey')
  axis(1, c(0:5), 10^c(0:5))
  axis(1, -1, '<0.1')
  abline(h = c(5.5, 10.5, 15.5), lty = 2)
  
  points(log10(budget$SML.Flux), c(1:5), col = cols$cycles, pch = 15, cex = 1.5)
  points(log10(budget$SML.Flux), c(1:5), col = 'white', pch = 15, cex = 1)
  #points(log10(budget$SML.Flux), c(1:5), col = cols$cycles, pch = 22, cex = 1.5)
  #points(log10(budget$SML.Flux), c(1:5), pch = 1, cex = 1.5)
  add.error.bars.logx(budget$SML.Flux, budget$sSML.Flux, budget$Cycle, 0, 10)
  
  points(log10(budget$ST.N), c(1:5), col = cols$cycles, cex = 1.5, pch = 16)
  add.error.bars.logx(budget$ST.N, budget$sST.N, budget$Cycle, 0, 10)
  
  points(log10(budget$Ex), c(1:5)+5, col = cols$cycles, cex = 1.5, pch = 16)
  add.error.bars.logx(budget$Ex, budget$sEx, budget$Cycle+5, 0, 10)
  
  #for (i in 1:5) {
  #  lines(x = log10(mort$Mortality[mort$Cycle == i]), rep(i,2)+5, col = cols$obs2)
  #}
  #points(log10(mort$Mortality), c(1:5, 1:5)+5, col= cols$obs2, cex = 1, pch = 20)
  
  points(log10(budget$J.Nitrate), c(1:5)+10, col = cols$cycles, cex = 1.5, pch = c(16))
  points(log10(budget$J.Nitrate[1:3]), c(1:3)+10, col = 'white', cex = 1, pch = 16)
  add.error.bars.logx(budget$J.Nitrate, budget$J.Nitrate/2, budget$Cycle+10, 0, 10)
  #legend('bottomright', legend = c(1:5), col = col[1:5], title = 'Cycle', pch = 1, cex = 1.5)
  
  ## N2 fix
  lines(x = log10(c(0.71,4.7)), y = c(16,16), col = cols$cycles[1])
  lines(x = log10(c(0.9,0.97)), y = c(17,17), col = cols$cycles[2])
  lines(x = log10(c(1.79,4.54)), y = c(18,18), col = cols$cycles[3])
  lines(x = log10(c(0.03,1.05)), y = c(19,19), col = cols$cycles[4])
  lines(x = log10(c(0.08,0.75)), y = c(20,20), col = cols$cycles[5])
  
  points(x = log10(c(0.71,4.7)), y = c(16,16), col = cols$cycles[1], pch = 16, cex = 1.5)
  points(x = log10(c(0.9,0.97)), y = c(17,17), col = cols$cycles[2], pch = 16, cex = 1.5)
  points(x = log10(c(1.79,4.54)), y = c(18,18), col = cols$cycles[3], pch = 16, cex = 1.5)
  points(x = log10(c(0.03,1.05)), y = c(19,19), col = cols$cycles[4], pch = 16, cex = 1.5)
  points(x = log10(c(0.08,0.75)), y = c(20,20), col = cols$cycles[5], pch = 16, cex = 1.5)
  
  ## Angies N2 Fixation
  points(x = log10(63), y = 17, col = 'black', pch = 15, cex = 1.5)
  points(x = log10(63), y = 17, col = cols$cycles[2], pch = 15)
  points(x = -1, y = 18, col = cols$cycles[3], pch = 15, cex = 1.5)
  points(x = -1, y = 18, col = 'white', pch = 15, cex = 1)
  points(x = -1, y = 16, col = cols$cycles[1], pch = 15, cex = 1.5)
  points(x = -1, y = 16, col = 'white', pch = 15, cex = 1)
  points(x = -1, y = 19, col = cols$cycles[4], pch = 15, cex = 1.5)
  points(x = -1, y = 19, col = 'white', pch = 15, cex = 1)
  points(x = -1, y = 20, col = cols$cycles[5], pch = 15, cex = 1.5)
  points(x = -1, y = 20, col = 'white', pch = 15, cex = 1)
  
  
  mtext(c('a','b','c','d'), side =2 , las = 1, at = c(1,6,11,16), line = 0.5)
  
  par(new = T, plt = c(0.1, 0.9, 0.91, 0.96), mfrow = c(1,1))
  plot.new()
  legend('top', legend = c(1:5, ''), col = c(cols$cycles, 'white'), pch = 16, cex = 1.1, horiz = T, bg = NA, box.col = 'white')
  #mtext('Cycle', adj = 0.2, line = -1.6)
  
  dev.off()
}





R = c(0.7573443293866458, 1.2009578823463967, 1.7438894855397251, 2.3781837949642965, 3.7317354877583604, 2.8196641371619413, 4.474906045914127)
par = c(10.38525963149209, 49.24623115578049, 180.5695142378584, 299.83249581239863, 598.6599664991677, 896.1474036850993, 1097.152428810729)

plot(par, R)

fit = parameter.search(n = 6, progression = 6, splits = 20, cost = function(alpha, Ek, E, x) {(TheSource::model.Webb1974(alpha, Ek, E) - x)^2}, E = par, x = R, bounds = data.frame(min = c(0,1), max = c(0.01,1000)))
fit$min
lines(c(1:2e3), model.Webb1974(fit$min$alpha, fit$min$Ek, c(1:2e3)))






