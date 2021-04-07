calc.flux.nemuro = function(nemuro, lon, lat, depth, zone = 0) {
  flux = data.frame(Time = conv.time.matlab(nemuro$TC[2]), Zone = zone, Var = c('PON', 'Org', 'DON', 'DIN', 'Water'),
                    North = 0, South = 0, West = 0, East = 0,
                    North.cyclo = 0, South.cyclo = 0, West.cyclo = 0, East.cyclo = 0,
                    North.anti = 0, South.anti = 0, West.anti = 0, East.anti = 0,
                    Sink = 0, Vertical = 0, River = 0, Mix = 0, Total = 0, Inventory = 0)
  
  lx = which(nemuro$XC >= lon[1] & nemuro$XC < lon[2])
  ly = which(nemuro$YC >= lat[1] & nemuro$YC < lat[2])
  ld = which(nemuro$ZC >= depth[1] & nemuro$ZC < depth[2])
  lx = lx[lx<length(nemuro$XC)]
  ly = ly[ly<length(nemuro$YC)]
  ld = ld[ld<length(nemuro$ZC)]
  dz = diff(nemuro$ZC)
  
  ## PON
  flux$North[1] = -sum(nemuro$PON.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South[1] = sum(nemuro$PON.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West[1] =  sum(nemuro$PON.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East[1] =  -sum(nemuro$PON.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.cyclo[1] = -sum((nemuro$Eddy[lx, max(ly)+1] == 2) * nemuro$PON.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.cyclo[1] = sum((nemuro$Eddy[lx, min(ly)] == 2) * nemuro$PON.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.cyclo[1] =  sum((nemuro$Eddy[min(lx), ly] == 2) * nemuro$PON.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.cyclo[1] =  -sum((nemuro$Eddy[max(lx)+1, ly] == 2) * nemuro$PON.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.anti[1] = -sum((nemuro$Eddy[lx, max(ly)+1] == -2) * nemuro$PON.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.anti[1] = sum((nemuro$Eddy[lx, min(ly)] == -2) * nemuro$PON.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.anti[1] =  sum((nemuro$Eddy[min(lx), ly] == -2) * nemuro$PON.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.anti[1] =  -sum((nemuro$Eddy[max(lx)+1, ly] == -2) * nemuro$PON.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  ## Org
  flux$North[2] = -sum(nemuro$Org.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South[2] = sum(nemuro$Org.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West[2] =  sum(nemuro$Org.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East[2] =  -sum(nemuro$Org.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.cyclo[2] = -sum((nemuro$Eddy[lx, max(ly)+1] == 2) * nemuro$Org.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.cyclo[2] = sum((nemuro$Eddy[lx, min(ly)] == 2) * nemuro$Org.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.cyclo[2] =  sum((nemuro$Eddy[min(lx), ly] == 2) * nemuro$Org.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.cyclo[2] =  -sum((nemuro$Eddy[max(lx)+1, ly] == 2) * nemuro$Org.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.anti[2] = -sum((nemuro$Eddy[lx, max(ly)+1] == -2) * nemuro$Org.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.anti[2] = sum((nemuro$Eddy[lx, min(ly)] == -2) * nemuro$Org.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.anti[2] =  sum((nemuro$Eddy[min(lx), ly] == -2) * nemuro$Org.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.anti[2] =  -sum((nemuro$Eddy[max(lx)+1, ly] == -2) * nemuro$Org.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  ## DON
  flux$North[3] = -sum(nemuro$DON.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South[3] = sum(nemuro$DON.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West[3] =  sum(nemuro$DON.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East[3] =  -sum(nemuro$DON.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.cyclo[3] = -sum((nemuro$Eddy[lx, max(ly)+1] == 2) * nemuro$DON.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.cyclo[3] = sum((nemuro$Eddy[lx, min(ly)] == 2) * nemuro$DON.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.cyclo[3] =  sum((nemuro$Eddy[min(lx), ly] == 2) * nemuro$DON.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.cyclo[3] =  -sum((nemuro$Eddy[max(lx)+1, ly] == 2) * nemuro$DON.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.anti[3] = -sum((nemuro$Eddy[lx, max(ly)+1] == -2) * nemuro$DON.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.anti[3] = sum((nemuro$Eddy[lx, min(ly)] == -2) * nemuro$DON.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.anti[3] =  sum((nemuro$Eddy[min(lx), ly] == -2) * nemuro$DON.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.anti[3] =  -sum((nemuro$Eddy[max(lx)+1, ly] == -2) * nemuro$DON.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  ## DIN
  flux$North[4] = -sum(nemuro$DIN.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South[4] = sum(nemuro$DIN.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West[4] =  sum(nemuro$DIN.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East[4] =  -sum(nemuro$DIN.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.cyclo[4] = -sum((nemuro$Eddy[lx, max(ly)+1] == 2) * nemuro$DIN.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.cyclo[4] = sum((nemuro$Eddy[lx, min(ly)] == 2) * nemuro$DIN.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.cyclo[4] =  sum((nemuro$Eddy[min(lx), ly] == 2) * nemuro$DIN.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.cyclo[4] =  -sum((nemuro$Eddy[max(lx)+1, ly] == 2) * nemuro$DIN.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  flux$North.anti[4] = -sum((nemuro$Eddy[lx, max(ly)+1] == -2) * nemuro$DIN.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South.anti[4] = sum((nemuro$Eddy[lx, min(ly)] == -2) * nemuro$DIN.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West.anti[4] =  sum((nemuro$Eddy[min(lx), ly] == -2) * nemuro$DIN.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East.anti[4] =  -sum((nemuro$Eddy[max(lx)+1, ly] == -2) * nemuro$DIN.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  vol = 0
  for (lld in ld) {
    flux$Inventory[1] = flux$Inventory[1] + sum(nemuro$PON[lx,ly,lld] * dz[lld] * nemuro$domain$dx[lx, ly] * nemuro$domain$dy[lx, ly], na.rm = T)
    flux$Inventory[2] = flux$Inventory[2] + sum(nemuro$Org[lx,ly,lld] * dz[lld] * nemuro$domain$dx[lx, ly] * nemuro$domain$dy[lx, ly], na.rm = T)
    flux$Inventory[3] = flux$Inventory[3] + sum(nemuro$DON[lx,ly,lld] * dz[lld] * nemuro$domain$dx[lx, ly] * nemuro$domain$dy[lx, ly], na.rm = T)
    flux$Inventory[4] = flux$Inventory[4] + sum(nemuro$DIN[lx,ly,lld] * dz[lld] * nemuro$domain$dx[lx, ly] * nemuro$domain$dy[lx, ly], na.rm = T)
    
    vol = vol + sum(nemuro$domain$dx[lx,ly] * nemuro$domain$dy[lx,ly] * dz[lld] * nemuro$Mask[lx, ly, lld])
  }
  
  ## Water
  flux$North[5] = -sum(nemuro$Water.Flux.V[lx, max(ly)+1, ld], na.rm = T)
  flux$South[5] = sum(nemuro$Water.Flux.V[lx, min(ly), ld], na.rm = T)
  flux$West[5] =  sum(nemuro$Water.Flux.U[min(lx), ly, ld], na.rm = T)
  flux$East[5] =  -sum(nemuro$Water.Flux.U[max(lx)+1, ly, ld], na.rm = T)
  
  ## Sinking
  flux$Sink[1] =  sum(nemuro$PON.Flux.S[lx, ly, max(ld)], na.rm = T)
  flux$Sink[4] =  sum(nemuro$DIN.Flux.S[lx, ly, max(ld)], na.rm = T)
  flux$River[4] = sum(nemuro$river[lx,ly])
  
  ## Vertical Advction
  flux$Vertical[1] = sum(nemuro$PON.Flux.W[lx, ly, max(ld)], na.rm = T)
  flux$Vertical[2] = sum(nemuro$Org.Flux.W[lx, ly, max(ld)], na.rm = T)
  flux$Vertical[3] = sum(nemuro$DON.Flux.W[lx, ly, max(ld)], na.rm = T)
  flux$Vertical[4] = sum(nemuro$DIN.Flux.W[lx, ly, max(ld)], na.rm = T)
  flux$Vertical[5] = sum(nemuro$Water.Flux.W[lx, ly, max(ld)], na.rm = T)
  
  ## Vertical Mixing
  flux$Mix[1] = sum(-nemuro$Kz[lx,ly,lld] * (nemuro$PON[lx,ly,lld] - nemuro$PON[lx,ly,lld+1]) / dz[lld] * (nemuro$domain$dx[lx,ly] * nemuro$domain$dy[lx,ly]), na.rm = T) * 86400
  flux$Mix[2] = sum(-nemuro$Kz[lx,ly,lld] * (nemuro$Org[lx,ly,lld] - nemuro$Org[lx,ly,lld+1]) / dz[lld] * (nemuro$domain$dx[lx,ly] * nemuro$domain$dy[lx,ly]), na.rm = T) * 86400
  flux$Mix[3] = sum(-nemuro$Kz[lx,ly,lld] * (nemuro$DON[lx,ly,lld] - nemuro$DON[lx,ly,lld+1]) / dz[lld] * (nemuro$domain$dx[lx,ly] * nemuro$domain$dy[lx,ly]), na.rm = T) * 86400
  flux$Mix[4] = sum(-nemuro$Kz[lx,ly,lld] * (nemuro$DIN[lx,ly,lld] - nemuro$DIN[lx,ly,lld+1]) / dz[lld] * (nemuro$domain$dx[lx,ly] * nemuro$domain$dy[lx,ly]), na.rm = T) * 86400
  
  ## Totals
  flux$Total = flux$North + flux$South + flux$West + flux$East + flux$Sink + flux$Vertical + flux$River + flux$Mix
  
  ## Normalize all values by vol
  flux[,c(4:ncol(flux))] = flux[,c(4:ncol(flux))] / (vol + 1e-6)
  
  flux
}




calc.flux.poc = function(satellite, lon, lat, zone = 0) {
  flux = data.frame(Time = satellite$times$mid, Zone = zone, Var = 'SSPON',
                    North = 0, South = 0, West = 0, East = 0,
                    North.cyclo = 0, South.cyclo = 0, West.cyclo = 0, East.cyclo = 0,
                    North.anti = 0, South.anti = 0, West.anti = 0, East.anti = 0,
                    Sink = 0, Vertical = 0, River = 0, Mix = 0, Total = 0, Inventory = 0)
  
  lx = which(satellite$lon >= lon[1] & satellite$lon < lon[2])
  ly = which(satellite$lat >= lat[1] & satellite$lat < lat[2])
  lx = lx[lx<length(satellite$lon)]
  ly = ly[ly<length(satellite$lat)]
  
  if (length(lx) == 0 | length(ly) == 0) {
    return(list(flux = flux, sum = sum(flux[1,2:9])))
  }
  
  ## PON
  flux$North[1] = -sum(satellite$adv$V[lx, min(ly)], na.rm = T)
  flux$South[1] = sum(satellite$adv$V[lx, max(ly)], na.rm = T)
  flux$West[1] =  sum(satellite$adv$U[min(lx), ly], na.rm = T)
  flux$East[1] =  -sum(satellite$adv$U[max(lx), ly], na.rm = T)
  flux$Inventory[1] = sum(satellite$field[lx,ly] * satellite$dx[lx, ly] * satellite$dy[lx, ly], na.rm = T)
  
  vol = sum(satellite$dx[lx,ly] * satellite$dy[lx,ly])
  flux[,c(4:ncol(flux))] = flux[,c(4:ncol(flux))] / vol
  flux$Total = sum(flux[4:10], na.rm = T)
  flux
}



calc.flux.doc = function(CDOM, lon, lat, zone = 0) {
  flux = data.frame(Time = CDOM$times$mid, Zone = zone, Var = 'SSDON', North = 0, South = 0, West = 0, East = 0,
                    Sink = 0, Vertical = 0, River = 0, Mix = 0, Total = 0, Inventory = 0)
  
  lx = which(CDOM$lon >= lon[1] & CDOM$lon < lon[2])
  ly = which(CDOM$lat >= lat[1] & CDOM$lat < lat[2])
  lx = lx[lx<length(CDOM$lon)]
  ly = ly[ly<length(CDOM$lat)]
  
  ## PON
  flux$North[1] = -sum(CDOM$adv$V[lx, min(ly)], na.rm = T) / 15
  flux$South[1] = sum(CDOM$adv$V[lx, max(ly)], na.rm = T) / 15
  flux$West[1] =  sum(CDOM$adv$U[min(lx), ly], na.rm = T) / 15
  flux$East[1] =  -sum(CDOM$adv$U[max(lx), ly], na.rm = T) / 15
  flux$Inventory[1] = sum(CDOM$MLR2[lx,ly] * CDOM$dx[lx, ly] * CDOM$dy[lx, ly], na.rm = T) / 15
  
  vol = sum(CDOM$dx[lx,ly] * CDOM$dy[lx,ly])
  flux[,c(4:ncol(flux))] = flux[,c(4:ncol(flux))] / vol
  flux$Total = sum(flux[4:10], na.rm = T)
  flux
}


add.map.quiver.sparse = function(x, y, u, v, N, zscale = 1, col = 'black') {
  x = as.numeric(x)
  y = as.numeric(y)
  u = as.numeric(u)
  v = as.numeric(v)
  if (length(x) < length(u)) {
    grid = expand.grid(x = x, y = y)
    if (nrow(grid) != length(u)) {
      stop('Incompatible lengths.')
    } 
    x = grid$x
    y = grid$y
  }
  
  l = sample(1:length(x), size = N)
  add.map.quiver(x[l], y[l], u[l], v[l], zscale = zscale, col = col)
}




add.boxplot.box2 = function(x, y, col = 'grey', border = 'black', width = 0.7, lty = 1, lcol = 'black',
                            lwd = 1, xlwd = NULL, outliers = T, pcol = 'black', pch = 1, cex = 1) {
  if (length(x) == 1) { x = rep(x, length(y))}
  for (xx in unique(x)) {
    l = which(x == xx)
    if(is.null(xlwd)) { xlwd = width / 3 }
    
    ## statistics
    q1 = quantile(y[l], probs = 0.25, na.rm = T)
    q3 = quantile(y[l], probs = 0.75, na.rm = T)
    iqr = IQR(y[l], na.rm = TRUE)
    m = median(y[l], na.rm = TRUE)
    
    ## Box
    rect(ybottom = xx - width/2, xleft = q1, ytop = xx + width/2, xright = q3, col = col, border = border)
    lines(y = c(xx - width/2, xx + width/2), x = rep(m,2)) # Horizontal
    
    ## Add outliers
    k = which(y[l] < q1 - 1.5 * iqr | y[l] > q3 + 1.5 * iqr)
    if (length(k) > 0 & outliers) { points(y = rep(xx, length(k)), x = y[l[k]], pch = pch, col = pcol, cex = cex) }
    
    ## Add whiskers
    if (length(k) > 0) {
      lines(y = rep(xx, 2), x = c(q1, min(y[l[-k]])), col = lcol, lwd = lwd)
      lines(y = rep(xx, 2), x = c(q3, max(y[l[-k]])), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(min(y[l[-k]]), 2), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(max(y[l[-k]]), 2), col = lcol, lwd = lwd)
    } else {
      lines(y = rep(xx, 2), x = c(q1, min(y[l])), col = lcol, lwd = lwd)
      lines(y = rep(xx, 2), x = c(q3, max(y[l])), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(min(y[l]), 2), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(max(y[l]), 2), col = lcol, lwd = lwd)
    }
  }
}
