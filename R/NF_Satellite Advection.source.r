
#### Oscar + Satellite Functions

load.oscar = function(file, lats = NULL, lons = NULL) {
  
  f = ncdf4::nc_open(file)
  oscar = f
  u = ncdf4::ncvar_get(f, 'u')
  v = ncdf4::ncvar_get(f, 'v')
  nc_close(f)
  
  ## Set names
  times = as.POSIXct(oscar$var$u$dim[[4]]$vals * 86400, origin = '1992-10-05 00:00:00', tz = 'GMT')
  dimnames(u) = list(longitude = oscar$var$u$dim[[1]]$vals, latitude = oscar$var$u$dim[[2]]$vals, time = times)
  dimnames(v) = list(longitude = oscar$var$u$dim[[1]]$vals, latitude = oscar$var$u$dim[[2]]$vals, time = times)
  
  oscar = list(u = u, v = v, lon = oscar$var$u$dim[[1]]$vals, lat = oscar$var$u$dim[[2]]$vals, time = times,
               grid = expand.grid(lon = oscar$var$u$dim[[1]]$vals, lat = oscar$var$u$dim[[2]]$vals))
  
  if (!is.null(lats)) {
    
    l = which(oscar$lat >= lats[1] & oscar$lat <= lats[2])
    oscar$u = oscar$u[,l,]
    oscar$v = oscar$v[,l,]
    oscar$lat = oscar$lat[l]
    oscar$grid = expand.grid(lon = oscar$lon, lat = oscar$lat)
  }
  if (!is.null(lons)) {
    lons[lons<0] = lons[lons<0] + 360
    ## Filter
    l = which(oscar$lon >= lons[1] & oscar$lon <= lons[2])
    oscar$u = oscar$u[l,,]
    oscar$v = oscar$v[l,,]
    oscar$lon = oscar$lon[l]
    oscar$grid = expand.grid(lon = oscar$lon, lat = oscar$lat)
  }
  oscar
}



get.vel = function(lon, lat, time, oscar) {
  
  grid = expand.grid(lon = lon, lat = lat)
  grid$lon[grid$lon < 0] = grid$lon[grid$lon < 0] + 360
  grid$u = NA
  grid$v = NA
  
  if (time < min(oscar$time) | time > max(oscar$time)) {
    stop(paste0('Improper OSCAR product loaded for time = ', time))
  }
  t2 = min(which(oscar$time > time))
  t1 = t2-1
  t1w = as.numeric(oscar$time[t2] - time) / as.numeric(oscar$time[t2] - oscar$time[t1])
  
  for (i in 1:nrow(grid)) {
    ## Find Indicies
    x2 = min(which(oscar$lon > grid$lon[i]))
    x1 = x2 - 1                         ## Makes assumptions about no boundaries
    y2 = max(which(oscar$lat > grid$lat[i]))
    y1 = y2 + 1
    
    
    if (x1 == 0 | y2 == 0 | t1 == 0 | x2 > length(oscar$lon) | y1 > length(oscar$lat) | t2 > length(oscar$time)){
      message('Point out of bounds! x1 = ', x1, ', y1 = ', y1, ', t1 = ', t1)
      grid$u[i] = 0
      grid$v[i] = 0
    } else {
    
      ## Calculate weights
      x1w = (oscar$lon[x2] - grid$lon[i]) / (oscar$lon[x2] - oscar$lon[x1])
      y1w = (oscar$lat[y2] - grid$lat[i]) / (oscar$lat[y2] - oscar$lat[y1])
      
      ## calculate u
      u1 = (oscar$u[x1,y1,t1] * x1w + oscar$u[x2,y1,t1] * (1 - x1w)) * y1w + (oscar$u[x1,y2,t1] * x1w + oscar$u[x2,y2,t1] * (1 - x1w)) * (1 - y1w)
      u2 = (oscar$u[x1,y1,t2] * x1w + oscar$u[x2,y1,t2] * (1 - x1w)) * y1w + (oscar$u[x1,y2,t2] * x1w + oscar$u[x2,y2,t2] * (1 - x1w)) * (1 - y1w)
      u = u1 * t1w + u2 * (1 - t1w)
      
      ## calculate v
      v1 = (oscar$v[x1,y1,t1] * x1w + oscar$v[x2,y1,t1] * (1 - x1w)) * y1w + (oscar$v[x1,y2,t1] * x1w + oscar$v[x2,y2,t1] * (1 - x1w)) * (1 - y1w)
      v2 = (oscar$v[x1,y1,t2] * x1w + oscar$v[x2,y1,t2] * (1 - x1w)) * y1w + (oscar$v[x1,y2,t2] * x1w + oscar$v[x2,y2,t2] * (1 - x1w)) * (1 - y1w)
      v = v1 * t1w + v2 * (1 - t1w)
      
      grid$u[i] = u
      grid$v[i] = v
      
      if (is.na(u)) { grid$u[i] = 0 }
      if (is.na(v)) { grid$v[i] = 0 }
    }
  }
  list(u = matrix(grid$u, nrow = length(lon), ncol = length(lat)),
       v = matrix(grid$v, nrow = length(lon), ncol = length(lat)))
}


advect = function(lat, lon, time, oscar, dt) {
  
  ## Initial Guess
  latt = lat[length(lat)]
  lont = lon[length(lon)]
  timet = time[length(time)]
  vel = get.vel(lont, latt, timet, oscar)
  
  ## p2 based on velocity at final pos at dt/2
  dlat = vel$v * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lat2 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lon2 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat2 / 180 * pi))
  time2 = timet + 86400 * dt/2
  vel2 = get.vel(lon2, lat2, time2, oscar)
  
  ## p3 based on velocity at p2 at dt/2 
  dlat = vel2$v * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lat3 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel2$u * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lon3 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat3 / 180 * pi))
  time3 = timet + 86400 * dt/2
  vel3 = get.vel(lon3, lat3, time3, oscar)
  
  ## p4 based on velocity at p2 at dt/2 
  dlat = vel$v * 86.4 * dt / 6378.14 # distance / earth radius
  lat4 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt / 6378.14 # distance / earth radius
  lon4 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat4 / 180 * pi))
  time4 = timet + 86400 * dt
  vel4 = get.vel(lon4, lat4, time4, oscar)
  
  ## RK4 velocity
  vel$u = (vel$u + 2 * vel2$u + 2 * vel3$u + vel4$u)/ 6
  vel$v = (vel$v + 2 * vel2$v + 2 * vel3$v + vel4$v)/ 6
  
  dlat = vel$v * 86.4 * dt / 6378.14 # distance / earth radius
  lat.final = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt / 6378.14 # distance / earth radius
  lon.final = lon[length(lon)] + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat.final / 180 * pi))
  
  list(lon = c(lon, lon.final), lat = c(lat, lat.final), time = c(time, time4))
}


## Particle advection routine
get.trajectory = function(cycle = 1, delta.t = 14, dt = -0.5, n = 100, oscar = oscar2017, with.poc = F) {
  out = list()
  
  ## Get relavent drifter entries
  l = which(drifter$Cycle == cycle & drifter$Device == 'OHM-MS-0001')
  
  for (k in l[seq(1, length(l), length.out = n)]) {
    message(k, ' ', appendLF = F)
    pos = list(lon = drifter$Lon[k], lat = drifter$Lat[k], time = drifter$Datetime[k])
    
    for (i in 1:(delta.t/abs(dt))) {
      pos = advect(lat = pos$lat, lon = pos$lon, time = pos$time, oscar = oscar, dt = dt)
    }
    out[[length(out) + 1]] = pos
  }
  message(' Done.')
  
  if (with.poc) {
    ## Add POC
    for (i in 1:length(out)) {
      out[[i]]$poc = rep(NA, length(out[[i]]$lon))
      for (j in 1:length(out[[i]]$poc)) {
        poc = get.field(out[[i]]$time[j], summary[summary$Param == 'POC.poc.4km',])
        k = which.min((out[[i]]$lon[j] - poc$grid()$lon)^2 + (out[[i]]$lat[j] - poc$grid()$lat)^2)
        out[[i]]$poc[j] = as.numeric(poc$field)[k]
      }
    }
  }
  
  ## Return
  out
}


add.trajectory = function(pos) {
  ## Assumes input object is output object of get.trajectory()
  for (i in 1:length(pos)) {
    
    ## mark starting point
    add.map.points(pos$lon, pos$lat, pch = 20, col = 'blue', cex = 0.6)
    
    ## Add path
    add.map.line(pos$lon, pos$lat, col = '#00008830', lwd = 0.5)
  }
}



get.field = function(time, summary, param = NULL, timeframe = NULL, conv = function(x){x/12/8.625}) {
  
  if (!is.null(param)) { summary = summary[summary$Param == param,] }
  if (!is.null(timeframe)) { summary = summary[summary$Timeframe == timeframe,] }
  
  if (nrow(summary) < 2) { stop('Too few satellite fields to use!') }
  
  if (!exists('index1')) { index1 <<- 0 }
  if (!exists('index2')) { index2 <<- 0 }
  
  summary = summary[order(summary$Datetime),] ## Sort by time
  l2 = which(summary$Datetime == summary$Datetime[min(which(summary$Datetime > time))])
  l1 = which(summary$Datetime == summary$Datetime[max(which(summary$Datetime <= time))])
  
  if (any(index1 != l1)) {
    index1 <<- l1
    message('Reading new satellite data, n = ', length(index1))
    satellite1 <<- read.satellite('Z:/Data/Satellite/Raw/', summary$File[l1], trim = T, lon = c(-99, -79), lat = c(18,31), conv = conv)
  }
  if (any(index2 != l2)) {
    index2 <<- l2
    message('Reading new satellite data, n = ', length(index2))
    satellite2 <<- read.satellite('Z:/Data/Satellite/Raw/', summary$File[l2], trim = T, lon = c(-99, -79), lat = c(18,31), conv = conv)
  }
  
  w = (as.numeric(time) - as.numeric(satellite1$times$mid)) / (as.numeric(satellite2$times$mid) - as.numeric(satellite1$times$mid))
  out = satellite1
  out$field = satellite1$field * w + satellite2$field * (1-w) # umol N m-2 (assuming redfield)
  
  l.lon = c(1:floor(length(out$lon)/2)) * 2
  l.lat = c(1:floor(length(out$lat)/2)) * 2
  field = matrix(0, nrow = length(l.lon), ncol = length(l.lat))
  
  for (i in 1:length(l.lon)) {
    for (j in 1:length(l.lat)) {
      field[i, j] = mean(out$field[l.lon[c(i,i,i-1,i-1)], l.lat[c(j,j-1,j,j-1)]], na.rm = T)
    }
  }
  out$field = field
  out$lon = out$lon[l.lon]
  out$lat = out$lat[l.lat]
  out$grid = function() { expand.grid(lon = out$lon, lat = out$lat) }
  
  out$cloud.free = length(which(is.na(out$field))) / length(out$field)
  message(Sys.time(), ': Percent cloud free: ', format(out$cloud.free * 100, digits = 1), '%')
  
  out$times$start = conv.time.unix(as.numeric(satellite2$times$start) * w + as.numeric(satellite2$times$start) * (1-w))
  out$times$mid = conv.time.unix(as.numeric(satellite2$times$mid) * w + as.numeric(satellite2$times$mid) * (1-w))
  out$times$end = conv.time.unix(as.numeric(satellite2$times$end) * w + as.numeric(satellite2$times$end) * (1-w))
  out
}





get.year = function (x) { as.POSIXlt(x)$year + 1900 }
get.month = function (x) { as.numeric(format(x, '%m')) }
get.day = function (x) { as.numeric(format(x, '%d')) }



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
    if (length(k) > 0) { points(y = rep(xx, length(k)), x = y[l[k]], pch = pch, col = pcol, cex = cex) }
    
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



