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



calc.oscar.advdiff = function(satellite) {
  
  year = get.year(satellite$times$start)
  
  message(Sys.time(), ': Loading Oscar for ', year)
  oscar = load.oscar(paste0('Z:/Data/OSCAR/oscar_vel', year, '.nc'), lons = c(-100, -80), lats = c(17, 32))
  oscar$Year = year
  
  message(Sys.time(), ': Calculating physical terms.')
  physics = get.vel(satellite$lon, satellite$lat, time = satellite$times$mid, oscar = oscar) # m s-1
  
  u = satellite$field * physics$u * 86400 * satellite$dy # umol N d-1 m-1
  v = satellite$field * physics$v * 86400 * satellite$dx #  N d-1 m-1

  message(Sys.time(), ': Done.')
  list(U = u, V = v)
}



