library(TheSource)
library(openxlsx)
library(pals)


cols = list(obs = '#4A9586',
            obs2 = '#A5D3CA',
            cycles = c('#4A9586', '#8DC7BB', '#DCEDEA', '#D54FD5', '#E697E6'),
            satellite = '#25A0C5',
            satellite2 = '#8ED6EA',
            nemuro = '#D73E68',
            nemuro2 = '#F0B9C8')

## Load Drifter positions
#drifter = get_gyre(start_date = NF17()[1], end_date = NF17()[2])
#drifter = rbind(drifter, get_gyre(start_date = NF18()[1], end_date = NF18()[2]))

#drifter$Datetime = conv.time.unix(drifter$Datetime)
#drifter = data.frame(Device = drifter$DeviceName, Datetime = drifter$Datetime, Lat = drifter$Latitude, Lon = drifter$Longitude, Wet = drifter$SubmergedBoolean)

#save(drifter, file = 'Positions/NF_Drifter Positions.rdata')
#write.xlsx(drifter, file = 'Positions/NF_Drifter Positions.xlsx')

load('_rdata/NF_Drifter Positions.rdata')


## Ship Positions
#pos = read.table('Positions/NF17 Ship Pos.csv', header = TRUE)
#pos18 = read.table('Positions/NF18 Ship Pos.txt', header = TRUE)


{## Plot map FIGURE 1
  pdf('_figures/NF_Figure1 Map.pdf')
  par(plt = c(0.2,0.8,0.1,0.8), mfrow = c(1,1))
  map = make.map('coastlineWorldFine', lon.min = -93, lon.max = -82, lat.min = 20, lat.max = 30.5, land.col = '#414141',
                 p = '+proj=aea +lat_1=24 +lat_2=29 +lon_0=-90', dlon = 3, dlat = 3)
  #bathy = get.bathy(map, 2)
  add.map.bathy.shade(map, bathy, zlim = c(-4e3, 200), filled = F)
  #add.map.bathy(map, bathy.gom, levels = c(-20, -55, -200, -1e3, -2e3), bathy.col = greyscale(7)[-1])
  add.map.bathy(map, bathy, levels = c(-200, -2e3), bathy.col = greyscale(5)[-1], drawlabels = F, bathy.lty = c(3,1))
  
  ## Cycle 1
  l = which(drifter$Device == 'OHM-MS-0001' & drifter$Datetime >= NF17.1()[1] & drifter$Datetime < NF17.1()[2])
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = 'black', pch = 1, cex = 1.2)
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = cols$cycles[1], pch = 20)
  ## Cycle 2
  l = which(drifter$Device == 'OHM-MS-0001' & drifter$Datetime >= NF17.2()[1] & drifter$Datetime < NF17.2()[2])
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = 'black', pch = 1, cex = 1.2)
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = cols$cycles[2], pch = 20)
  ## Cycle 3
  l = which(drifter$Device == 'OHM-MS-0001' & drifter$Datetime >= NF17.3()[1] & drifter$Datetime < NF17.3()[2])
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = 'black', pch = 1, cex = 1.2)
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = cols$cycles[3], pch = 20)
  ## Cycle 4
  l = which(drifter$Device == 'OHM-MS-0002' & drifter$Datetime >= NF18.1()[1] & drifter$Datetime < NF18.1()[2])
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = 'black', pch = 1, cex = 1.2)
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = cols$cycles[4], pch = 20)
  ## Cycle 5
  l = which(drifter$Device == 'OHM-MS-0001' & drifter$Datetime >= NF18.2()[1] & drifter$Datetime < NF18.2()[2])
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = 'black', pch = 1, cex = 1.2)
  add.map.points(drifter$Lon[l], drifter$Lat[l], col = cols$cycles[5], pch = 20)
  
  add.map.line(lon = c(-90.5,-90.5,-86.5,-86.5,-85.5,-85.5,-86.5,-86.5,-88.5,-88.5,-90.5),
               lat = c(27.5,25.5,25.5,26.5,26.5,27.5,27.5,28.5,28.5,27.5,27.5), lwd = 2, col = '#FAFEFB')
  
  mapScalebar('topleft')
  #add.colorbar(-200, 4e3, pal = ocean.ice, rev = T, x.pos = 0.75, width = 0.025, height = 0.9)
  #redraw.map(map)
  
  dev.off()
}






data = read.xlsx('Data/NF_In Situ NO3 and CO3 2019.10.24.xlsx', sheet = 'Cycle Summary', startRow = 2)
data2 = read.xlsx('Data/NF_In Situ NO3 and CO3 2019.10.24.xlsx', sheet = 'Summary Data', startRow = 2)
chl = read.xlsx('Data/NF_Chl Samples.xlsx', sheet = 'Summary', startRow = 2)
load('_rdata/NF_CTD Data (final).rdata')
nut = read.xlsx('data/NF_Nutrient Data Final.xlsx', startRow = 2, sheet = 2)



add.error.bars.logx = function(x, s.x, y, s.y, base, col = 'black') {
  
  if (length(s.x) == 1) {
    s.x = rep(s.x, length(x))
  }
  if (length(s.y) == 1) {
    s.y = rep(s.y, length(y))
  }
  
  ## Remove NAs
  l = !is.na(x) & !is.na(y)
  x = x[l]
  s.x = s.x[l]
  y = y[l]
  s.y = s.y[l]
  
  for (i in 1:length(x)) {
    if(!is.na(s.y[i])) {
      lines(x = rep(log(x[i], base), 2),
            y = c(y[i] + s.y[i], y[i] - s.y[i]),
            col = col)
    }
    
    if(!is.na(s.x[i])) {
      lines(x = log(c(x[i] + s.x[i], max(x[i] - s.x[i], 1e-23)), base),
            y = rep(y[i], 2),
            col = col)
    }
  }
}

{
  #pdf('_figures/NF_Figure1 Other.pdf', width = 6, height = 7)
  ## Nitrate
  par(plt = c(0.2, 0.36, 0.2, 0.9))
  #plot(NULL, NULL, xlim = c(-1.5,1.5), ylim = c(210,0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)',
  #     xlab = 'Nitrate (uM N)', xaxt = 'n')
  #grid(ny = 4, nx = NA); add.log.axis(1, grid.major = T, grid = F); box()
  #mtext('a', side = 3, line = -1.6, adj = 0.95, cex = 1.3)
  #rect(-2, -10, -1, 500, col = '#00000020', border = F)
  #for (i in 1:5) {
  #  l = which(nut$Cycle == i)
  #  ap = approx(nut$Depth[l], log10(nut$Nitrate[l]), xout = c(1:max(nut$Depth[l])), rule = 2)
  #  lines(ap$y, ap$x, col = cols$cycles[i], lwd = 2)
  #}
  
  plot(NULL, NULL, xlim = c(0,16), ylim = c(250,0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)',
       xlab = 'Nitrate (uM)')
  grid();# add.log.axis(1, grid.major = T, grid = F); box()
  mtext('a', side = 3, line = -1.6, adj = 0.05, cex = 1.3)
  #rect(-2, -10, -1, 500, col = '#00000020', border = F)
  for (i in 1:5) {
    l = which(nut$Cycle == i)
    ap = approx(nut$Depth[l], (nut$Nitrate[l]), xout = c(1:max(nut$Depth[l])), rule = 2)
    lines(ap$y, ap$x, col = cols$cycles[i], lwd = 2)
  }
  box()
  
  ## Inset
  par(plt = c(0.28, 0.36, 0.62, 0.9), new = T)
  plot(NULL, NULL, xlim = c(0,0.8), ylim = c(100,0), yaxs = 'i', xaxs = 'i', ylab = '', xlab = '',
       yaxt = 'n')
  rect(-20, -100, 20, 100, col = 'white'); axis(side = 2, labels = F, at = c(0, 25, 50, 75, 100))
  grid(nx = 4, ny = 4); box()
  for (i in 1:5) {
    l = which(nut$Cycle == i)
    ap = approx(nut$Depth[l], (nut$Nitrate[l]), xout = c(1:100), rule = 2)
    lines(ap$y, ap$x, col = cols$cycles[i], lwd = 2)
  }
  box()
  
  #### Chloropyll
  par(plt = c(0.361, 0.54, 0.2, 0.9), new = T)
  plot(NULL, NULL, xlim = c(0,0.8), ylim = c(250, 0), yaxs = 'i', yaxt = 'n', xaxs = 'i', ylab = '', xlab = 'Chl\n(ug Chl L-1)')
  grid(); box()
  
  for (i in 1:5) {
    cast = unique(data2$Cast[data2$Cycle == i])
    l = which(downcast$Cast %in% cast)
    message(i, '  ', length(cast))
    profile = approx(round(downcast$Depth[l]), downcast$Chl[l], xout = seq(1, 340, by = 1), ties = 'mean')
    lines(profile$y, profile$x, col = cols$cycles[i], lwd = 2)
    #text(profile$y[which.max(profile$x)], max(profile$x), i, pos = 1, col = 'darkgrey', cex = 1.3)
  }
  mtext('b', side = 3, line = -1.6, adj = 0.95, cex = 1.3)
  box()
  
  
  ## NPP
  par(plt = c(0.541, 0.71, 0.2, 0.9), new = T)
  plot(NULL, NULL, xlim = c(0,600), ylim = c(250,0), yaxs = 'i', xaxs = 'i', ylab = '', yaxt = 'n',
       xlab = 'NPP\n(umol C m-3 d-1)', xaxt = 'n')
  axis(1, at = c(0,200,400))
  abline(v = c(200,400), lty = 3, col = 'grey')
  abline(h = c(50,100,150,200), lty = 3, col = 'grey')
  box(); mtext('c', side = 3, line = -1.6, adj = 0.05, cex = 1.3)
  
  for (i in 5:1) {
    dd = runif(1,-1,1)*0
    l = which(data$Cycle == i)
    lines(data$NPP[l]*1000, data$Depth[l]+dd, lwd = 3, col = cols$cycles[i])
    text(data$NPP[max(l)]*1000, data$Depth[max(l)]+2, i, pos = 1, col = 'darkgrey', cex = 1.3)
    
    points(data$NPP*1000, data$Depth+dd, pch = 20)
    add.error.bars(data$NPP[l]*1000, data$sNPP[l]*1000, data$Depth[l]+dd, 0, col = cols$obs)
  }
  box()
  
  par(plt = c(0.711, 0.88, 0.2, 0.9), new = T)
  plot(NULL, NULL, xlim = c(0,1.6e3), ylim = c(250,0), yaxt = 'n', yaxs = 'i', xaxs = 'i', ylab = '', xlab = 'Nitrogen Flux\n(umol N m-2 d-1)')
  grid(); box(); mtext('d', side = 3, line = -1.6, adj = 0.05, cex = 1.3)
  
  for (i in unique(st$Cycle)) {
    l = which(st$Cycle == i)
    lines(st$Norg[l]/15*1e3, st$Depth[l], lwd = 3, col = cols$cycles[i])
    points(st$Norg[l]/15*1e3, st$Depth[l], pch = 16)
    add.error.bars(x = st$Norg[l]/15*1e3, s.x = st$sNorg[l]/15*1e3, y = st$Depth[l], s.y = 0)
    text(x = st$Norg[min(l)]/15*1e3, st$Depth[min(l)], i, pos = 3, col = 'darkgrey', cex = 1.3)
  }
  
  par(new = T, plt = c(0.1,0.9,0.9,0.96))
  plot.new()
  legend('top', legend = c(1:5), col = cols$cycles, lwd = 3, cex = 1, horiz = T, bg = NA, box.col = 'white')
  #mtext('Cycle', adj = 0.15, line = -1.3)
  
  #dev.off()
}



pdf('_figures/Supplemental Figure 1 (part1).pdf')
{ 
  par(mfrow = c(1,2))
  ## Nitrate
  #par(plt = c(0.02, 0.9))
  plot(NULL, NULL, xlim = c(-2,1.8), ylim = c(200,0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)',
       xlab = 'Nitrate (uM N)', xaxt = 'n')
  add.log.axis(1, grid.major = T); box()
  mtext('a', side = 3, line = -1.6, adj = 0.05, cex = 1.3)
  
  for (i in 1:5) {
    l = which(nut$Cycle == i)
    ap = approx(nut$Depth[l], log10(nut$Nitrate[l]), xout = c(1:max(nut$Depth[l])), rule = 2)
    lines(ap$y, ap$x, col = cols$cycles[i], lwd = 2)
  }
  #add.error.bars.logx(data$Nitrate, data$sNitrate, data$Depth, 0, 10)
  
  
  ## Ammonium
  plot(NULL, NULL, xlim = c(-2,0), ylim = c(200,0), yaxs = 'i', xaxs = 'i', ylab = 'Depth (m)',
       xlab = 'Ammonium (uM N)', xaxt = 'n')
  grid(ny = 4, nx = NA); add.log.axis(1, grid.major = T, grid = T); box()
  mtext('b', side = 3, line = -1.6, adj = 0.05, cex = 1.3)
  
  for (i in 1:5) {
    l = which(data$Cycle == i)
    lines(log10(data$Ammonium[l]), data$Depth[l], col = cols$cycles[i], lwd = 2)
    text(log10(data$Ammonium[max(l)]), data$Depth[max(l)], i, pos = 1, col = 'darkgrey', cex = 1.3)
  }
  add.error.bars.logx(data$Ammonium, data$sAmmonium, data$Depth, 0, 10, col = '#cc5020')
  
}
dev.off()

#### Save Satellites
# Make satellite objects for use in the project.
# Include whether objects are inside or outside GoM.


## 2017
summary = rebuild.satellite.index('Data/Satellite/')
l = which(summary$Param == 'CHL' & summary$Year == 2017 & summary$Month == 5)[1]

box = data.frame(lat = c(30, 40, 30, 25, 21, 15, 15), lon = c(-100, -90, -82, -81, -87, -95, -100))
chl17 = read.satellite('Data/Satellite/', summary$File[l])
chl17$field[chl17$field == -9999] = NA
chl17 = trim.satellite(chl17, lon = c(-100, -70), lat = c(15,40))
chl17$grid = chl17$grid()
chl17$grid$inside = 0
for (i in 1:nrow(chl17$grid)) {
  chl17$grid$inside[i] = is.inside(c(chl17$grid$lat[i], chl17$grid$lon[i]), box)
}
chl17$grid$field = as.numeric(chl17$field)
chl17$grid$field[chl17$grid$inside == F] = NA

save(chl17, file = '_rdata/NF17_Satellite Chl.rdata')

## 2018
l = which(summary$Param == 'CHL' & summary$Year == 2018 & summary$Month == 5)[1]

chl18 = read.satellite('Data/Satellite/', summary$File[l])
chl18$field[chl18$field == -9999] = NA
chl18 = trim.satellite(chl18, lon = c(-100, -70), lat = c(15,40))
chl18$grid = chl18$grid()
chl18$grid$inside = 0
for (i in 1:nrow(chl18$grid)) {
  chl18$grid$inside[i] = is.inside(c(chl18$grid$lat[i], chl18$grid$lon[i]), box)
}
chl18$grid$field = as.numeric(chl18$field)
chl18$grid$field[chl18$grid$inside == F] = NA
save(chl18, file = '_rdata/NF17_Satellite Chl.rdata')


{
  par(mfrow = c(1,2))
  x = seq(0.01, 100, by = 0.01)
  y = ecdf(chl17$grid$field)
  plot(log10(x), y(x), type = 'l', xaxt = 'n', yaxs = 'i', ylim = c(0,1))
  add.log.axis()
  abline(h = y(0.1))
  
  y = ecdf(chl18$grid$field)
  plot(log10(x), y(x), type = 'l', xaxt = 'n', yaxs = 'i', ylim = c(0,1))
  add.log.axis()
  abline(h = y(0.1))
}









