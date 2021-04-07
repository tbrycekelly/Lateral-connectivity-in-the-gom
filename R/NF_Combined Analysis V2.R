library(TheSource)
library(ncdf4)
library(openxlsx)
#source('~/BTE/R/Offline Nemuro-GOM.R')
library(pals)
source('R/source.r')
#source('R/NF_Satellite Advection.source.r')
load('_rdata/NF_Drifter Positions.rdata')


cols = list(obs = '#4A9586',
            obs2 = '#A5D3CA',
            cycles = c('#4A958630', '#8DC7BB30', '#DCEDEA30', '#D54FD530', '#E697E630'),
            satellite = '#25A0C5',
            satellite2 = '#8ED6EA',
            nemuro = '#D73E68',
            nemuro2 = '#F0B9C8')


{## Plot map FIGURE 1
  par(plt = c(0.3,0.7,0.1,0.9), mfrow = c(1,1))
  map = make.map('coastlineWorldFine', lon.min = -91, lon.max = -85, lat.min = 22, lat.max = 30.5, land.col = '#414141',
                 p = '+proj=aea +lat_1=24 +lat_2=29 +lon_0=-90', dlon = 3, dlat = 3)
  add.map.bathy.shade(map, bathy.gom, zlim = c(-4e3, 200), refinement = 2)
  #add.map.bathy(map, bathy.gom, levels = c(-20, -55, -200, -1e3, -2e3), bathy.col = greyscale(7)[-1])
  add.map.bathy(map, bathy.gom, levels = c(-200, -2e3), bathy.col = greyscale(5)[-1], drawlabels = F, bathy.lty = c(3,1))
  
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
  add.colorbar(-200, 4e3, pal = ocean.ice, rev = T, x.pos = 0.75, width = 0.025, height = 0.9)
  #redraw.map(map)
}


#### Grid entire GOM
delta = 1
regions = expand.grid(lon = seq(-98, -80, by = delta), lat = seq(20, 30, by = delta))
regions$id = c(1:nrow(regions))

poc.files = list.files('Z:/Data/Taylor/_rdata', pattern = ' POC.8k', full.names = T)
doc.files = list.files('Z:/Data/Taylor/_rdata', pattern = 'CDOM.ADV', full.names = T)
nem.files = list.files('Z:/Data/Taylor/_rdata', pattern = 'NEMURO.', full.names = T)

## Initialize the main data object to hold all the flux data:
flux = data.frame(Time = make.time(), Zone = 0, Var = NA,
                  North = 0, South = 0, West = 0, East = 0,
                  North.cyclo = 0, South.cyclo = 0, West.cyclo = 0, East.cyclo = 0,
                  North.anti = 0, South.anti = 0, West.anti = 0, East.anti = 0,
                  Sink = 0, Vertical = 0, River = 0, Mix = 0, Total = 0, Inventory = 0)

{
  ## Load and calculate fluxes from NEMURO-GOM
  for (i in 1:length(nem.files)){
    message(Sys.time(), ': Loading NEMURO Advection file ', nem.files[i], ' (', i, '/', length(nem.files), ')')
    load(nem.files[i])
    d = c(0, 135)
    #d = c(55, 135)
    for (box in 1:nrow(regions)) {
      if (!is.null(nemuro$Kz)) {
        temp = calc.flux.nemuro(nemuro, lon = c(regions$lon[box], regions$lon[box]+delta), lat = c(regions$lat[box], regions$lat[box]+delta), depth = d, zone = box)
        flux = rbind(flux, temp) ## Concat
      } else {
        message(' NEMURO file w/o Kz: ', nem.files[i])
      }
    }
  }
  
  for (i in 1:length(poc.files)) {
      message(Sys.time(), ': Loading file ', poc.files[i], ' (', i, ' of ', length(poc.files), ')')
      load(poc.files[i])
      for (j in 1:nrow(regions)) {
        temp = calc.flux.poc(satellite, c(regions$lon[j], regions$lon[j] + delta), c(regions$lat[j], regions$lat[j] + delta), zone = j)
        flux = rbind(flux, temp) ## Concat
      }
  }
    
  flux = flux[-1,]
  flux = flux[flux$Inventory > 0,]
}

#save(flux, file = '_rdata/flux.nemuro.2020.09.15.rdata')


#load('_rdata/flux.2020.08.17.rdata')
#load('_rdata/flux.2020.08.17.rdata')


{ ## Perform integrations over the control region:
  l.north = c(142,143,163,164,146)
  l.south = c(123,124,125,126,146)
  l.west = c(163,142,123)
  l.east = c(164,146,127)
  l.vol = c(163,164,142,143,144,145,146,123,124,125,126)
  
  times = unique(flux$Time)
  
  SSPON = c()
  SSPONn = c()
  SSPONs = c()
  SSPONw = c()
  SSPONe = c()
  
  SSDON = c()
  NEMURO = c()
  NEMURO.up = c()
  NEMURO.pon = c()
  NEMURO.pon.cyclo = c()
  NEMURO.pon.anti = c()
  
  NEMURO.ponn = c()
  NEMURO.pons = c()
  NEMURO.pone = c()
  NEMURO.ponw = c()
  
  NEMURO.org = c()
  NEMURO.don = c()
  NEMURO.din = c()
  NEMURO.org.cyclo = c()
  NEMURO.don.cyclo = c()
  NEMURO.din.cyclo = c()
  NEMURO.org.anti = c()
  NEMURO.don.anti = c()
  NEMURO.din.anti = c()
  NEMURO.sink = c()
  
  for (i in 1:length(times)) {
    
    ## SSPON
    l = which(flux$Time == times[i] & flux$Var == 'SSPON')
    
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      temp = (s-n+w-e) / length(l.vol)
      SSPON = c(SSPON, temp)
    }
    
    
    ## NEMURO
    l = which(flux$Time == times[i] & flux$Var %in% c('Org', 'PON', 'DON', 'DIN'))
    
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      n.anti = sum(flux$North.anti[l[flux$Zone[l] %in% l.north]])
      s.anti = sum(flux$South.anti[l[flux$Zone[l] %in% l.south]])
      w.anti = sum(flux$West.anti[l[flux$Zone[l] %in% l.west]])
      e.anti = sum(flux$East.anti[l[flux$Zone[l] %in% l.east]])
      n.cyclo = sum(flux$North.cyclo[l[flux$Zone[l] %in% l.north]])
      s.cyclo = sum(flux$South.cyclo[l[flux$Zone[l] %in% l.south]])
      w.cyclo = sum(flux$West.cyclo[l[flux$Zone[l] %in% l.west]])
      e.cyclo = sum(flux$East.cyclo[l[flux$Zone[l] %in% l.east]])
      up = mean(flux$Vertical[l[flux$Zone[l] %in% l.vol]] + flux$Mix[l[flux$Zone[l] %in% l.vol]])
      sink = mean(flux$Sink[l[flux$Zone[l] %in% l.vol]])
      
      NEMURO = c(NEMURO, (s+n+w+e) / length(l.vol))
      NEMURO.up = c(NEMURO.up, up)
      NEMURO.sink = c(NEMURO.sink, sink)
    }
    
    ## DON
    l = which(flux$Time == times[i] & flux$Var == 'DON')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      n.anti = sum(flux$North.anti[l[flux$Zone[l] %in% l.north]])
      s.anti = sum(flux$South.anti[l[flux$Zone[l] %in% l.south]])
      w.anti = sum(flux$West.anti[l[flux$Zone[l] %in% l.west]])
      e.anti = sum(flux$East.anti[l[flux$Zone[l] %in% l.east]])
      n.cyclo = sum(flux$North.cyclo[l[flux$Zone[l] %in% l.north]])
      s.cyclo = sum(flux$South.cyclo[l[flux$Zone[l] %in% l.south]])
      w.cyclo = sum(flux$West.cyclo[l[flux$Zone[l] %in% l.west]])
      e.cyclo = sum(flux$East.cyclo[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.don = c(NEMURO.don, (s+n+w+e) / length(l.vol))
      NEMURO.don.cyclo = c(NEMURO.don.cyclo, (s.cyclo+n.cyclo+w.cyclo+e.cyclo) / length(l.vol))
      NEMURO.don.anti = c(NEMURO.don.anti, (s.anti+n.anti+w.anti+e.anti) / length(l.vol))
    }
    
    ## PON
    l = which(flux$Time == times[i] & flux$Var == 'PON')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      n.anti = sum(flux$North.anti[l[flux$Zone[l] %in% l.north]])
      s.anti = sum(flux$South.anti[l[flux$Zone[l] %in% l.south]])
      w.anti = sum(flux$West.anti[l[flux$Zone[l] %in% l.west]])
      e.anti = sum(flux$East.anti[l[flux$Zone[l] %in% l.east]])
      n.cyclo = sum(flux$North.cyclo[l[flux$Zone[l] %in% l.north]])
      s.cyclo = sum(flux$South.cyclo[l[flux$Zone[l] %in% l.south]])
      w.cyclo = sum(flux$West.cyclo[l[flux$Zone[l] %in% l.west]])
      e.cyclo = sum(flux$East.cyclo[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.pon = c(NEMURO.pon, (s+n+w+e) / length(l.vol))
      NEMURO.pon.anti = c(NEMURO.pon.anti, (s.anti+n.anti+w.anti+e.anti) / length(l.vol))
      NEMURO.pon.cyclo = c(NEMURO.pon.cyclo, (s.cyclo+n.cyclo+w.cyclo+e.cyclo) / length(l.vol))
      
      NEMURO.ponn = c(NEMURO.ponn, (s+n+w+e) * (s+n) / sqrt((s+n)^2+(w+e)^2))
      NEMURO.pone = c(NEMURO.pone, (s+n+w+e) * (w+e) / sqrt((s+n)^2+(w+e)^2))
    }
    
    ## Org
    l = which(flux$Time == times[i] & flux$Var == 'Org')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      n.anti = sum(flux$North.anti[l[flux$Zone[l] %in% l.north]])
      s.anti = sum(flux$South.anti[l[flux$Zone[l] %in% l.south]])
      w.anti = sum(flux$West.anti[l[flux$Zone[l] %in% l.west]])
      e.anti = sum(flux$East.anti[l[flux$Zone[l] %in% l.east]])
      n.cyclo = sum(flux$North.cyclo[l[flux$Zone[l] %in% l.north]])
      s.cyclo = sum(flux$South.cyclo[l[flux$Zone[l] %in% l.south]])
      w.cyclo = sum(flux$West.cyclo[l[flux$Zone[l] %in% l.west]])
      e.cyclo = sum(flux$East.cyclo[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.org = c(NEMURO.org, (s+n+w+e) / length(l.vol))
      NEMURO.org.cyclo = c(NEMURO.org.cyclo, (s.cyclo+n.cyclo+w.cyclo+e.cyclo) / length(l.vol))
      NEMURO.org.anti = c(NEMURO.org.anti, (s.anti+n.anti+w.anti+e.anti) / length(l.vol))
    }
    
    ## DIN
    l = which(flux$Time == times[i] & flux$Var == 'DIN')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      n.anti = sum(flux$North.anti[l[flux$Zone[l] %in% l.north]])
      s.anti = sum(flux$South.anti[l[flux$Zone[l] %in% l.south]])
      w.anti = sum(flux$West.anti[l[flux$Zone[l] %in% l.west]])
      e.anti = sum(flux$East.anti[l[flux$Zone[l] %in% l.east]])
      n.cyclo = sum(flux$North.cyclo[l[flux$Zone[l] %in% l.north]])
      s.cyclo = sum(flux$South.cyclo[l[flux$Zone[l] %in% l.south]])
      w.cyclo = sum(flux$West.cyclo[l[flux$Zone[l] %in% l.west]])
      e.cyclo = sum(flux$East.cyclo[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.din = c(NEMURO.din, (s+n+w+e) / length(l.vol))
      NEMURO.din.cyclo = c(NEMURO.din.cyclo, (s.cyclo+n.cyclo+w.cyclo+e.cyclo) / length(l.vol))
      NEMURO.din.anti = c(NEMURO.din.anti, (s.anti+n.anti+w.anti+e.anti) / length(l.vol))
    }
    
  }
}


f = function(x) {
  x = x[!is.na(x)]
  a = max(x[which(x < median(x) + 1.2 * IQR(x))])*55
  b = min(x[which(x > median(x) - 1 * IQR(x))])*55
  c = median(x)*55
  
  message(b, '\t', c, '\t', a)
}

f(NEMURO.org)
f(NEMURO.pon)
f(NEMURO.org + NEMURO.pon)
f(NEMURO.don)

f(NEMURO.org + NEMURO.pon + NEMURO.don + NEMURO.din)




{
  par(mfrow = c(2,2), plt = c(0.2,0.8, 0.25,0.8))
  
  ## DIN Flux
  plot(NEMURO.din*1e3*55, (NEMURO.din.cyclo + NEMURO.din.anti)*1e3*55, pch = 20, ylim = c(-100, 4e2), xlim = c(-100, 8e2),
       xlab = 'DIN Flux', ylab = 'Eddy-associated Flux'); grid(); box()
  
  abline(0,1,lty = 2); abline(0,0.1,lty = 2)
  
  ## DON Flux
  plot(NEMURO.don*1e3*55, (NEMURO.don.cyclo + NEMURO.don.anti) / NEMURO.don, pch = 20, ylim = c(-1, 1), xlim = c(-100, 8e3),
       xlab = 'DON Flux', ylab = 'Eddy-associated Flux'); grid(); box()
  abline(0,1,lty = 2); abline(0,0.1,lty = 2)
  
  ## PON Flux
  plot(NEMURO.pon*1e3*55 + NEMURO.org*1e3*55,
       (NEMURO.pon.cyclo + NEMURO.org.cyclo + NEMURO.pon.anti + NEMURO.org.anti)*1e3*55,
       pch = 20, ylim = c(-100, 4e3), xlim = c(-100, 8e3),
       xlab = 'PON Flux', ylab = 'Eddy-associated Flux'); grid(); box()
  abline(0,1,lty = 2); abline(0,0.1,lty = 2)
  
  ## Total Flux
  plot((NEMURO.pon + NEMURO.org + NEMURO.don + NEMURO.din)*1e3*55,
       (NEMURO.pon.cyclo + NEMURO.org.cyclo + NEMURO.pon.anti + NEMURO.org.anti + NEMURO.don.anti + NEMURO.don.cyclo + 
         NEMURO.din.anti + NEMURO.din.cyclo)*1e3*55, pch = 20, ylim = c(-100, 4e3), xlim = c(-100, 8e3),
       xlab = 'Total Flux', ylab = 'Eddy-associated Flux'); grid()
  abline(0,1,lty = 2); abline(0,0.1,lty = 2)
}


summary((NEMURO.org.cyclo + NEMURO.org.anti) / NEMURO.org)
summary((NEMURO.din.cyclo + NEMURO.din.anti) / NEMURO.din)
summary((NEMURO.don.cyclo + NEMURO.don.anti) / NEMURO.don)
summary((NEMURO.pon.cyclo + NEMURO.pon.anti) / NEMURO.pon)
summary((NEMURO.pon.cyclo + NEMURO.pon.anti + NEMURO.pon.cyclo + NEMURO.pon.anti + NEMURO.din.cyclo + NEMURO.din.anti + NEMURO.org.cyclo + NEMURO.org.anti)*1e3*55)
eddy = (NEMURO.pon.cyclo + NEMURO.pon.anti + NEMURO.pon.cyclo + NEMURO.pon.anti + NEMURO.din.cyclo + NEMURO.din.anti + NEMURO.org.cyclo + NEMURO.org.anti)

{
  par(mfrow = c(1,1), plt = c(0.3, 0.8, 0.2, 0.9))
  plot(NULL, NULL, xlim = c(-100,100), ylim = c(10, 0), xlab = 'Contribution to Flux (%)', ylab = '', yaxt = 'n')
  grid(); box()
  
  add.violin2(2, 100*(NEMURO.din.cyclo + NEMURO.din.anti) / NEMURO.din, '#30303020',
              make.pal(median(NEMURO.din.cyclo / (NEMURO.din.cyclo + NEMURO.din.anti)), min = 0, max = 1, pal = 'ocean.balance'))
  
  add.violin2(4, 100*(NEMURO.don.cyclo + NEMURO.don.anti) / NEMURO.don, '#30303020',
              make.pal(median(NEMURO.don.cyclo / (NEMURO.don.cyclo + NEMURO.don.anti)), min = 0, max = 1, pal = 'ocean.balance'))
  
  add.violin2(6, 100*(NEMURO.pon.cyclo + NEMURO.org.cyclo + NEMURO.pon.anti + NEMURO.org.anti) / (NEMURO.pon + NEMURO.org), '#30303020',
              make.pal(median((NEMURO.pon.cyclo + NEMURO.org.cyclo) / (NEMURO.pon.cyclo + NEMURO.org.cyclo + NEMURO.pon.anti + NEMURO.org.anti)), min = 0, max = 1, pal = 'ocean.balance'))
  
  add.violin2(8, 100*(NEMURO.pon.cyclo + NEMURO.org.cyclo + NEMURO.pon.anti + NEMURO.org.anti + NEMURO.don.anti + NEMURO.don.cyclo + 
                    NEMURO.din.anti + NEMURO.din.cyclo) / (NEMURO.pon + NEMURO.org + NEMURO.don + NEMURO.din), '#30303020',
              make.pal(median((NEMURO.pon.cyclo + NEMURO.org.cyclo + NEMURO.don.cyclo + NEMURO.din.cyclo) /
                              (NEMURO.pon.cyclo + NEMURO.org.cyclo + NEMURO.pon.anti + NEMURO.org.anti + NEMURO.don.anti +
                                 NEMURO.don.cyclo + NEMURO.din.anti + NEMURO.din.cyclo)), min = 0, max = 1, pal = 'ocean.balance'))
  
  add.colorbar(0,100,pal = 'ocean.balance', labels = c(0,50,100), width = 0.02, y.pos = 0.55)
}

obs = read.xlsx('Data/NF_SedTrap Data.xlsx', startRow = 2)
J = read.xlsx('Data/NF_Nitrogen Budget.xlsx', startRow = 2)

{
  #pdf('_figures/NF_Figure 3.pdf')
  par(plt = c(0.3,0.6, 0.2, 0.9), xpd = F)
  col = pals::brewer.paired(12)
  
  plot(NULL, NULL, xlim = c(-5e3, 10e3), xaxs = 'i', ylim = c(14,1), ylab = '', xlab = 'Flux (umol N m-2 d-1)', yaxt = 'n')
  axis(2, at = c(2,4,6,11, 13), labels = c('Obs Export', 'Satellite PON', 'NEMURO Lateral', 'Obs Vertical', 'NEMURO Vertical'), las = 1)
  abline(v = 0, lty = 1, col = 'darkblue')
  abline(v = c(-1,1,2,3,4)*1e3, lty = 3, col = 'grey')
  abline(h = c(2,4,6,11, 13), lty = 3, col = 'grey')
  
  density.SSPON = density(SSPON)
  density.NEMURO = density(NEMURO)
  density.NEMUROpon = density(NEMURO.pon + NEMURO.org)
  density.NEMUROdon = density(NEMURO.don)
  density.NEMUROdin = density(NEMURO.din)
  density.NEMUROup = density(NEMURO.up)
  density.NEMUROsink = density(NEMURO.sink)
  
  ## Obs
  l = which(obs$Depth < 70)
  y = 2 + obs$Cycle[l]/5-0.5
  points(obs$Norg[l]/14*1e3, y = y, pch = 15, cex = 1.2, col = cols$obs)
  points(obs$Norg[l]/14*1e3, y = y, pch = 22, cex = 1.3)
  add.error.bars(obs$Norg[l]/14*1e3, obs$σNorg[l]/14*1e3, y, 0)
  
  ## satellite
  add.violin(4, SSPON, density.SSPON, cols$satellite, cols$satellite2, 0.75)
  
  ## NEMURO
  add.violin(6, NEMURO, density.NEMURO, cols$nemuro, cols$nemuro2, 0.75)
  
  add.violin(7, NEMURO.pon+NEMURO.org, density.NEMUROpon, cols$nemuro, cols$nemuro2, 0.3)
  add.violin(8, NEMURO.don, density.NEMUROdon, cols$nemuro, cols$nemuro2, 0.3)
  add.violin(9, NEMURO.din, density.NEMUROdin, cols$nemuro, cols$nemuro2, 0.3)
  
  mtext(side = 2, text = 'Lateral PON', line = 1, at = 7, font = 3, col = 'darkgrey', las = 1)
  mtext(side = 2, text = 'Lateral DON', line = 1, at = 8, font = 3, col = 'darkgrey', las = 1)
  mtext(side = 2, text = 'Lateral DIN', line = 1, at = 9, font = 3, col = 'darkgrey', las = 1)
  
  ## Obs up
  points(J$J.Nitrate, y = y+9, pch = 15, cex = 1.2, col = cols$obs)
  points(J$J.Nitrate, y = y+9, pch = 22, cex = 1.3)
  add.error.bars(J$J.Nitrate, J$J.Nitrate, y+9, 0)
  
  ## NEMURO UP
  add.violin(13, NEMURO.up, density.NEMUROup,  cols$nemuro, cols$nemuro2, 0.75)
  par(xpd=TRUE)
  legend('top', legend = c('Field Data', '', 'Satellite-based', 'NEMURO-based'),
         col = c(cols$obs, cols$obs2, cols$satellite, cols$nemuro), pch = c(15,NA,-0x1f3bbL,-0x1f3bbL),
         box.col = NA, bg = NA, ncol = 2, inset = -0.13, cex = 0.9, pt.cex = 1.2)
  box()
  
  #dev.off()
}

add.violin = function(y, data, density, col1, col2, scale = 1, n = 1e4) {
  
  temp = approx(density$x, density$y, xout = seq(min(density$x), max(density$x), length.out = n))
  density$x = temp$x*1e3*55
  density$y = temp$y
  xx = c(density$x, rev(density$x))
  yy  = y + scale * c(density$y, -rev(density$y)) / mean(density$y) / 10
  
  polygon(x = xx, y = yy, col = col1, border = NA)
  add.boxplot.box2(y, data*1e3*55, col = col2, width = scale*0.8, outliers = F)
}

{
  pdf('_figures/NF_Supplemental Figure 5 (violin plot).pdf')
  par(plt = c(0.3,0.6, 0.2, 0.9), xpd = F)
  col = pals::brewer.paired(12)
  
  plot(NULL, NULL, xlim = c(-5e3, 10e3), xaxs = 'i', ylim = c(14,1), ylab = '', xlab = 'Flux (umol N m-2 d-1)', yaxt = 'n')
  axis(2, at = c(2,4,6,11, 13), labels = c('Obs Export', 'Satellite PON', 'NEMURO Lateral', 'Obs Vertical', 'NEMURO Vertical'), las = 1)
  abline(v = 0, lty = 1, col = 'darkblue')
  abline(v = c(-4,-2,2,4,6,8)*1e3, lty = 3, col = 'grey')
  abline(h = c(2,4,6,11, 13), lty = 3, col = 'grey')
  
  density.SSPON = density(SSPON)
  density.NEMURO = density(NEMURO)
  density.NEMUROpon = density(NEMURO.pon + NEMURO.org)
  density.NEMUROdon = density(NEMURO.don)
  density.NEMUROdin = density(NEMURO.din)
  density.NEMUROup = density(NEMURO.up)
  density.NEMUROsink = density(NEMURO.sink)
  
  ## Obs
  l = which(obs$Depth > 70 & obs$Depth < 150)
  y = 2 + obs$Cycle[l]/5-0.5
  points(obs$Norg[l]/14*1e3, y = y, pch = 15, cex = 1.2, col = cols$obs)
  points(obs$Norg[l]/14*1e3, y = y, pch = 22, cex = 1.3)
  add.error.bars(obs$Norg[l]/14*1e3, obs$σNorg[l]/14*1e3, y, 0)
  
  ## satellite
  add.violin(4, SSPON, density.SSPON, cols$satellite, cols$satellite2, 0.75)
  
  ## NEMURO
  add.violin(6, NEMURO, density.NEMURO, cols$nemuro, cols$nemuro2, 0.5)
  
  add.violin(7, NEMURO.pon+NEMURO.org, density.NEMUROpon, cols$nemuro, cols$nemuro2, 0.3)
  add.violin(8, NEMURO.don, density.NEMUROdon, cols$nemuro, cols$nemuro2, 0.3)
  add.violin(9, NEMURO.din, density.NEMUROdin, cols$nemuro, cols$nemuro2, 0.3)
  
  mtext(side = 2, text = 'Lateral PON', line = 1, at = 7, font = 3, col = 'darkgrey', las = 1)
  mtext(side = 2, text = 'Lateral DON', line = 1, at = 8, font = 3, col = 'darkgrey', las = 1)
  mtext(side = 2, text = 'Lateral DIN', line = 1, at = 9, font = 3, col = 'darkgrey', las = 1)
  
  ## Obs up
  points(J$J.Nitrate, y = y+9, pch = 15, cex = 1.2, col = cols$obs)
  points(J$J.Nitrate, y = y+9, pch = 22, cex = 1.3)
  add.error.bars(J$J.Nitrate, J$J.Nitrate, y+9, 0)
  
  ## NEMURO UP
  add.violin(13, NEMURO.up, density.NEMUROup,  cols$nemuro, cols$nemuro2, 0.75)
  par(xpd=TRUE)
  legend('top', legend = c('Field Data', '', 'Satellite-based', 'NEMURO-based'),
         col = c(cols$obs, cols$obs2, cols$satellite, cols$nemuro), pch = c(15,NA,-0x1f3bbL,-0x1f3bbL),
         box.col = NA, bg = NA, ncol = 2, inset = -0.13, cex = 0.9, pt.cex = 1.2)
  box()
  
  dev.off()
}


add.violin2 = function(y, data, col1, col2, scale = 1, n = 1e4) {
  density = density(data)
  temp = approx(density$x, density$y, xout = seq(min(density$x), max(density$x), length.out = n))
  density$x = temp$x
  density$y = temp$y
  xx = c(density$x, rev(density$x))
  yy  = y + scale * c(density$y, -rev(density$y)) / mean(density$y) / 100
  
  polygon(x = xx, y = yy, col = col1, border = NA)
  add.boxplot.box2(y, data, col = col2, width = scale*0.8, outliers = F)
}




#### 120 m depth horizon:

{
  l.north = c(142,143,163,164,146)
  l.south = c(123,124,125,126,146)
  l.west = c(163,142,123)
  l.east = c(164,146,127)
  l.vol = c(163,164,142,143,144,145,146,123,124,125,126)
  
  times = unique(flux$Time)
  
  SSPON = c()
  SSDON = c()
  NEMURO = c()
  NEMURO.up = c()
  NEMURO.pon = c()
  NEMURO.org = c()
  NEMURO.don = c()
  NEMURO.din = c()
  NEMURO.sink = c()
  
  for (i in 1:length(times)) {
    
    ## SSPON
    l = which(flux$Time == times[i] & flux$Var == 'SSPON')
    
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      
      SSPON = c(SSPON, (s-n+w-e) / length(l.vol))
    }
    
    ## SSDON
    l = which(flux$Time == times[i] & flux$Var == 'SSDON')
    
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      
      SSDON = c(SSDON, (s-n+w-e) / length(l.vol))
    }
    
    
    ## NEMURO
    l = which(flux$Time == times[i] & flux$Var %in% c('Org', 'PON', 'DON', 'DIN'))
    
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      up = mean(flux$Vertical[l[flux$Zone[l] %in% l.vol]] + flux$Mix[l[flux$Zone[l] %in% l.vol]])
      sink = mean(flux$Sink[l[flux$Zone[l] %in% l.vol]])
      
      NEMURO = c(NEMURO, (s-n+w-e) / length(l.vol))
      NEMURO.up = c(NEMURO.up, up)
      NEMURO.sink = c(NEMURO.sink, sink)
    }
    
    ## DON
    l = which(flux$Time == times[i] & flux$Var == 'DON')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.don = c(NEMURO.don, (s-n+w-e) / length(l.vol))
    }
    
    ## PON
    l = which(flux$Time == times[i] & flux$Var == 'PON')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.pon = c(NEMURO.pon, (s-n+w-e) / length(l.vol))
    }
    
    ## Org
    l = which(flux$Time == times[i] & flux$Var == 'Org')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.org = c(NEMURO.org, (s-n+w-e) / length(l.vol))
    }
    
    ## DIN
    l = which(flux$Time == times[i] & flux$Var == 'DIN')
    if (length(l) > 0) {
      n = sum(flux$North[l[flux$Zone[l] %in% l.north]])
      s = sum(flux$South[l[flux$Zone[l] %in% l.south]])
      w = sum(flux$West[l[flux$Zone[l] %in% l.west]])
      e = sum(flux$East[l[flux$Zone[l] %in% l.east]])
      
      NEMURO.din = c(NEMURO.din, (s-n+w-e) / length(l.vol))
    }
    
  }
}

{ ## Plot for 120m traps
  par(plt = c(0.3,0.8, 0.2, 0.9))
  col = pals::brewer.paired(10)
  
  plot(NULL, NULL, xlim = c(-0.3, 0.3), ylim = c(13,1), ylab = '', xlab = 'Flux (mmol N m-3 d-1)', yaxt = 'n')
  axis(2, at = c(2,4:8,10,12), labels = c('Obs Export', 'NEMURO Lateral', '"PON', '"Org',
                                          '"DON', '"DIN', 'NEMURO Vertical', 'NEMURO Sinking'), las = 1)
  abline(v = 0, lty = 1, col = 'darkblue')
  abline(v = c(-0.1,-0.05, 0.05, 0.1), lty = 3, col = 'grey')
  abline(h = c(2,4,6,12), lty = 3, col = 'grey')
  
  #density.SSPON = density(SSPON)
  #density.SSDON = density(SSDON)
  density.NEMURO = density(NEMURO)
  density.NEMUROpon = density(NEMURO.pon)
  density.NEMUROdon = density(NEMURO.don)
  density.NEMUROdin = density(NEMURO.din)
  density.NEMUROorg = density(NEMURO.org)
  density.NEMUROup = density(NEMURO.up)
  density.NEMUROsink = density(NEMURO.sink)
  
  ## Obs
  l = which(obs$Depth > 70 & obs$Depth < 180)
  y = runif(length(l), -0.25, 0.25) + rep(2, length(l))
  points(obs$Norg[l]/14/obs$Depth[l], y = y, pch = 1, cex = 1.2)
  add.error.bars(obs$Norg[l]/14/obs$Depth[l], obs$σNorg[l]/14/obs$Depth[l], y, 0)
  
  ## NEMURO
  add.violin(4, NEMURO, density.NEMURO, col[6], col[5], 0.75)
  
  add.violin(5, NEMURO.pon, density.NEMUROpon, 'grey', col[5], 0.5)
  add.violin(6, NEMURO.org, density.NEMUROorg, 'grey', col[5], 0.5)
  add.violin(7, NEMURO.don, density.NEMUROdon, 'grey', col[5], 0.5)
  add.violin(8, NEMURO.din, density.NEMUROdin, 'grey', col[5], 0.5)
  
  ## NEMURO UP
  add.violin(10, NEMURO.up, density.NEMUROup, col[6], col[5], 0.75)
  add.violin(12, NEMURO.sink, density.NEMUROsink, col[6], col[5], 0.75)
}



probs = 0.25
regional.flux = data.frame(Zone = regions$id, Lateral = NA, Vertical = NA, SSPON = NA, SSDON = NA,
                           North = NA, South = NA, East = NA, West = NA, SSNorth = NA, SSSouth = NA, SSEast = NA, SSWest = NA)

for (i in 1:nrow(regional.flux)) {
  ## NEMURO-GOM
  l = which(flux$Zone == regional.flux$Zone[i] & flux$Var %in% c('PON', 'DON', 'Org', 'DIN'))
  regional.flux$Lateral[i] = quantile(flux$North[l] + flux$South[l] + flux$West[l] + flux$East[l] + flux$River[l], probs = probs)
  regional.flux$Vertical[i] = quantile(flux$Vertical[l] + flux$Mix[l], probs = probs)
  
  regional.flux$North[i] = quantile(flux$North[l], probs = probs)
  regional.flux$South[i] = quantile(flux$South[l], probs = probs)
  regional.flux$East[i] = quantile(flux$East[l], probs = probs)
  regional.flux$West[i] = quantile(flux$West[l], probs = probs)
}


pdf('_figures/NF_NEMURO Meam State by Year (2x2).pdf', width = 12, height = 10)
par(mfrow = c(2,2), plt = c(0.1,0.9,0.1,0.8))
for (i in c(1993, 1998, 2005, 2010)) {
  load(paste0('Z:/Data/Taylor/_rdata/NEMUROAVG.', i, '.rdata')) ## Load nemuro.avg

## Supplemental Figure of N Inventories
  map = make.map('coastlineWorldFine', lon.min = -95, lon.max = -83, lat.min = 19, lat.max = 30.5,
                 land.col = '#414141', p = '+proj=aea +lat_1=24 +lat_2=29 +lon_0=-90',
                 dlon = 5, dlat = 5)
  #add.map.layer(regions$lon[regional.flux$Zone], regions$lat[regional.flux$Zone], regional.flux$Vertical, pal = 'ocean.balance', zlim = c(-1,1)/100)
  add.map.layer(nemuro.avg$XC, nemuro.avg$YC, log10(nemuro.avg$Org[,,1] + nemuro.avg$DIN[,,1] + nemuro.avg$DON[,,1] + nemuro.avg$PON[,,1]),
                pal = 'inferno', zlim = log10(c(0.1, 10)), indicate = F)
  
  redraw.map(map)
  
  ## Add River Points
  #l = which(nemuro.avg$river > 100000)
  #add.map.points(nemuro.avg$domain$grid$lon[l], nemuro.avg$domain$grid$lat[l], cex = make.cex(log10(as.numeric(nemuro.avg$river[l]))), col = 'black')
  
  add.map.line(lon = c(-90.5,-90.5,-86.5,-86.5,-85.5,-85.5,-86.5,-86.5,-88.5,-88.5,-90.5),
               lat = c(27.5,25.5,25.5,26.5,26.5,27.5,27.5,28.5,28.5,27.5,27.5), lwd = 2, col = '#FAFEFB')
  
  
  
  ## Labels and colorbar
  #mtext(side = 3, adj = 1, cex = 0.8, paste0('Surface N Concentration in May ', i, ' (mmol N m-3)'))
  add.colorbar(0.1, 10, pal = 'inferno', log = T, base = 10, horizontal = T, labels = c(1e-1, 1, 10),
               ticks= c(c(1:9)*1e-1, c(1:9)))
  
}
dev.off()


pdf('_figures/NF_NEMURO Meam Flux by Year (zint).pdf')
for (i in c(1993, 1994, 1996:2012)) {
  load(paste0('Z:/Data/Taylor/_rdata/NEMUROAVG.', i, '.rdata')) ## Load nemuro.avg
  
  ## Supplemental Figure of N Inventories
  par(mfrow = c(1,1), plt = c(0.1,0.9,0.1,0.8))
  map = make.map('coastlineWorldFine', lon.min = -95, lon.max = -83, lat.min = 19, lat.max = 30.5,
                 land.col = '#414141', p = '+proj=aea +lat_1=24 +lat_2=29 +lon_0=-90',
                 dlon = 5, dlat = 5)
  
  add.map.layer(nemuro.avg$XC, nemuro.avg$YC,
                sqrt(apply((nemuro.avg$PON.Flux.U + nemuro.avg$DON.Flux.U + nemuro.avg$DIN.Flux.U + nemuro.avg$Org.Flux.U)[,,1:6], c(1,2), sum)^2 +
                       apply((nemuro.avg$PON.Flux.V + nemuro.avg$DON.Flux.V + nemuro.avg$DIN.Flux.V + nemuro.avg$Org.Flux.V)[,,1:6], c(1,2), sum)^2) /
                  (nemuro.avg$domain$dx * nemuro.avg$domain$dy),
                pal = 'inferno', zlim = c(0, 1000), indicate = F)
  
  
  redraw.map(map)
  
  ## Add River Points
  #l = which(nemuro.avg$river > 100000)
  #add.map.points(nemuro.avg$domain$grid$lon[l], nemuro.avg$domain$grid$lat[l], cex = make.cex(log10(as.numeric(nemuro.avg$river[l]))), col = 'black')
  
  add.map.line(lon = c(-90.5,-90.5,-86.5,-86.5,-85.5,-85.5,-86.5,-86.5,-88.5,-88.5,-90.5),
               lat = c(27.5,25.5,25.5,26.5,26.5,27.5,27.5,28.5,28.5,27.5,27.5), lwd = 2, col = '#FAFEFB')
  
  
  
  ## Labels and colorbar
  mtext(side = 3, adj = 1, cex = 0.8, paste0('Surface N Flux in May ', i, ' (mmol N m-3)'))
  add.colorbar(0, 1000, pal = 'inferno', horizontal = T, labels = c(0, 500, 1000), ticks = c(1:20)*100)

}
dev.off()


#### Supplemental Figure 3 (new)
pdf('_figures/NF_NEMURO Meam Flux by Year (2x2).pdf', width = 12, height = 10)

par(mfrow = c(2,2), plt = c(0.1,0.9,0.1,0.8))
for (i in c(1993,1998,2005,2010)) {
  load(paste0('Z:/Data/Taylor/_rdata/NEMUROAVG.', i, '.rdata')) ## Load nemuro.avg
  
  ## Supplemental Figure of N Inventories
  map = make.map('coastlineWorldFine', lon.min = -95, lon.max = -83, lat.min = 19, lat.max = 30.5,
                 land.col = '#414141', p = '+proj=aea +lat_1=24 +lat_2=29 +lon_0=-90',
                 dlon = 5, dlat = 5)
  
  add.map.layer(nemuro.avg$XC, nemuro.avg$YC,
                sqrt(apply((nemuro.avg$PON.Flux.U + nemuro.avg$DON.Flux.U + nemuro.avg$DIN.Flux.U + nemuro.avg$Org.Flux.U)[,,1:6], c(1,2), sum)^2 +
                       apply((nemuro.avg$PON.Flux.V + nemuro.avg$DON.Flux.V + nemuro.avg$DIN.Flux.V + nemuro.avg$Org.Flux.V)[,,1:6], c(1,2), sum)^2) /
                  (nemuro.avg$domain$dx * nemuro.avg$domain$dy),
                pal = 'inferno', zlim = c(0, 1000), indicate = F, filled = T)
  
  
  redraw.map(map)
  
  ## Add River Points
  #l = which(nemuro.avg$river > 100000)
  #add.map.points(nemuro.avg$domain$grid$lon[l], nemuro.avg$domain$grid$lat[l], cex = make.cex(log10(as.numeric(nemuro.avg$river[l]))), col = 'black')
  
  add.map.line(lon = c(-90.5,-90.5,-86.5,-86.5,-85.5,-85.5,-86.5,-86.5,-88.5,-88.5,-90.5),
               lat = c(27.5,25.5,25.5,26.5,26.5,27.5,27.5,28.5,28.5,27.5,27.5), lwd = 2, col = '#FAFEFB')
  
  mtext(side = 3, adj = 1, cex = 0.8, paste0('Surface N Flux in May ', i, ' (mmol N m-3)'))
}
add.colorbar(0, 1000, pal = 'inferno', horizontal = T, labels = c(0, 500, 1000), ticks = c(1:20)*100)
dev.off()


#### Supplemental Figure 4 (new)
pdf('_figures/NF_Satellite Meam POC by Year (2x2).pdf', width = 12, height = 10)

par(mfrow = c(2,2), plt = c(0.1,0.9,0.1,0.8))
for (i in c(21,41,69,73)) {
  #load(paste0('Z:/Data/Taylor/_rdata/Satellite AVGPOC.8k.', i, '.rdata')) ## Load nemuro.avg
  load(poc.files[i])
  ## Supplemental Figure of N Inventories
  map = make.map('coastlineWorldFine', lon.min = -95, lon.max = -83, lat.min = 19, lat.max = 30.5,
                 land.col = '#414141', p = '+proj=aea +lat_1=24 +lat_2=29 +lon_0=-90',
                 dlon = 5, dlat = 5)
  
  add.map.layer(satellite$lon, satellite$lat, log10(satellite$field*6.625), pal = 'inferno', zlim = c(0, 1), indicate = F, filled = T)
  add.map.layer(satellite$lon, satellite$lat, log10(satellite$field*6.625), pal = 'inferno', zlim = c(0, 1), indicate = F, filled = T)
  add.map.layer(satellite$lon, satellite$lat, log10(satellite$field*6.625), pal = 'inferno', zlim = c(0, 1), indicate = F, filled = T)
  
  
  redraw.map(map)
  
  ## Add River Points
  #l = which(nemuro.avg$river > 100000)
  #add.map.points(nemuro.avg$domain$grid$lon[l], nemuro.avg$domain$grid$lat[l], cex = make.cex(log10(as.numeric(nemuro.avg$river[l]))), col = 'black')
  
  add.map.line(lon = c(-90.5,-90.5,-86.5,-86.5,-85.5,-85.5,-86.5,-86.5,-88.5,-88.5,-90.5),
               lat = c(27.5,25.5,25.5,26.5,26.5,27.5,27.5,28.5,28.5,27.5,27.5), lwd = 2, col = '#FAFEFB')
  
  #mtext(side = 3, adj = 1, cex = 0.8, paste0('Surface N Flux in May ', i, ' (mmol N m-3)'))
  
  add.colorbar(1, 10, log = T, base = 10, pal = 'inferno', horizontal = T, labels = c(1, 10), ticks = c(1:10))
}
dev.off()


#### Calculate annual mean state from NEMURO-GOM

nem.files = list.files('Z:/Data/Taylor/_rdata', pattern = 'NEMURO.', full.names = T)

for (year in 1994:2014) {
  l = which(grepl(year, nem.files))
  
  if (length(l) > 1) {
    message(Sys.time(), ': Starting average for ', year, ' (n = ', length(l),'). ', appendLF = F)
    
    load(nem.files[l[1]])
    nemuro.avg = nemuro
    message('1  ', appendLF = F)
    
    ## For each day
    for (i in 2:length(l)) {
      load(nem.files[l[i]])
      for (k in c(5:7,9:18,20:36)) {
        nemuro.avg[[k]] = nemuro.avg[[k]] + nemuro[[k]]
      }
      
      message(i, '  ', appendLF = F)
    }
    
    ## Divide by n
    for (k in c(5:7,9:18,20:36)) {
      nemuro.avg[[k]] = nemuro.avg[[k]] / length(l)
    }
    message()
    message(Sys.time(), ': Saving file.')
    save(nemuro.avg, file = paste0('Z:/Data/Taylor/_rdata/NEMUROAVG.', year, '.rdata'))
  }
}




for (year in 2000:2019) {
  l = which(grepl(year, poc.files))
  
  if (length(l) > 1) {
    message(Sys.time(), ': Starting average for ', year, ' (n = ', length(l),'). ', appendLF = F)
    
    load(poc.files[l[1]])
    satellite.avg = satellite
    message('1  ', appendLF = F)
    
    ## For each day
    for (i in 2:length(l)) {
      load(poc.files[l[i]])
      satellite.avg$field = satellite.avg$field + satellite$field
      
      message(i, '  ', appendLF = F)
    }
    
    satellite.avg$field = satellite.avg$field / length(l)
    message()
    message(Sys.time(), ': Saving file.')
    save(satellite.avg, file = paste0('Z:/Data/Taylor/_rdata/Satellite AVGPOC.8k.', year, '.rdata'))
  }
}

a = c()
for (y in unique(year)) {
  l = which(year == y)
  a = c(a, mean((NEMURO.don[l] / (NEMURO.don + NEMURO.pon + NEMURO.org + NEMURO.din)[l])))
}



box = data.frame(lon = c(-90.5,-90.5,-86.5,-86.5,-85.5,-85.5,-86.5,-86.5,-88.5,-88.5,-90.5),
                 lat = c(27.5,25.5,25.5,26.5,26.5,27.5,27.5,28.5,28.5,27.5,27.5))

l = which(is.inside(grid, box))
summary(as.numeric(sqrt(nemuro.avg$U[,,1]^2 + nemuro.avg$V[,,1]^2))[l])

box2 = box
box2$lon = box2$lon %% 360
l = which(is.inside(oscar$grid, box2))

summary(as.numeric(sqrt(oscar$u[,,33]^2 + oscar$v[,,33]^2))[l])



{
  plot(NULL, NULL, xlim = c(-0.1, 0.1), ylim = c(0,2000), yaxt = 'n', xaxt = 'n', ylab = '', xlab = '')
  n = density(SSPONn, n = 2^14)
  s = density(SSPONs, n = 2^14, bw = 0.00025)
  e = density(SSPONe, n = 2^14, bw = 0.00025)
  w = density(SSPONw, n = 2^14, bw = 0.00025)
  
  lines(n$x, n$y+1500, lwd = 2, col = 'darkred')
  lines(s$x, s$y, lwd = 2, col = 'darkgreen')
  lines(e$x, e$y+1000, lwd = 2, col = 'darkblue')
  lines(w$x, w$y+500, lwd = 2, col = 'black')
  abline(v = 0, lty = 3)
  abline(h = c(0, 500,1000,1500), lty = 3)
  
  text(-0.0045, 1800, 'North', col = 'darkred', cex = 2)
  text(-0.0045, 1300, 'East', col = 'darkgreen', cex = 2)
  text(-0.0045, 800, 'West', col = 'darkblue', cex = 2)
  text(-0.0045, 300, 'South', col = 'black', cex = 2)
}


{
  plot(SSPONe*55, SSPONn*55, col = cols$satellite,
       pch = 16, ylim = c(-10, 80), xlim = c(-10, 80), ylab = 'V Component', xlab= 'U Component')
  
  points(NEMURO.pone * 1000, NEMURO.ponn * 1000, col = paste0(cols$nemuro, '10'), pch = 16)
  grid()
  abline(v = 0, col = 'grey')
  abline(h = 0, col = 'grey')
  
}

