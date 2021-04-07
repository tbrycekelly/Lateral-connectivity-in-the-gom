load.river = function(dis.file, nut.file) {
  
  finfo = file.info(dis.file)
  dis.file = file(description = dis.file, 'rb')
  dis = readBin(dis.file, double(), size=4, n = finfo$size, endian="big")
  close(dis.file)
  
  finfo = file.info(nut.file)
  nut.file = file(description = nut.file, 'rb')
  nut = readBin(nut.file, double(), size=4, n = finfo$size, endian="big")
  close(nut.file)
  
  dis = array(dis, dim = c(536,380,31)) # m s-1?
  nut = array(nut, dim = c(536,380,31))
  dis = dis[,,] * 86400 # m3 d-1
  nut = nut[,,] # mmol m-3
  
  dis * nut
  
}


load.physics = function(model,
                        U = 'Z:/Data/Taylor/gom_data/input_off/uvel/U.0000000000.data',
                        V = 'Z:/Data/Taylor/gom_data/input_off/vvel/V.0000000000.data',
                        W = 'Z:/Data/Taylor/gom_data/input_off/wvel/W.0000000000.data',
                        S = 'Z:/Data/Taylor/gom_data/input_off/sal/S.0000000000.data',
                        T = 'Z:/Data/Taylor/gom_data/input_off/tmp/T.0000000000.data') {
  
  ## Get file info
  finfo = file.info(S)
  
  ## read file
  u.file = file(description = U, 'rb')
  U = readBin(u.file, double(), size=4, n = finfo$size, endian="big")
  close(u.file)
  
  v.file = file(description = V, 'rb')
  V = readBin(v.file, double(), size=4, n = finfo$size, endian="big")
  close(v.file)
  
  w.file = file(description = W, 'rb')
  W = readBin(w.file, double(), size=4, n = finfo$size, endian="big")
  close(w.file)
  
  s.file = file(description = S, 'rb')
  S = readBin(s.file, double(), size=4, n = finfo$size, endian="big")
  close(s.file)
  
  t.file = file(description = T, 'rb')
  T = readBin(t.file, double(), size=4, n = finfo$size, endian="big")
  close(t.file)
  
  ## format data
  list(U = array(U, dim = c(length(model$XC), length(model$YC), length(model$ZC))),
       V = array(V, dim = c(length(model$XC), length(model$YC), length(model$ZC))),
       W = array(W, dim = c(length(model$XC), length(model$YC), length(model$ZC))),
       S = array(S, dim = c(length(model$XC), length(model$YC), length(model$ZC))),
       T = array(T, dim = c(length(model$XC), length(model$YC), length(model$ZC)))
  )
}


## Light data
load.lightfield = function(model, file = 'Z:/Data/Taylor/gom_data/initial_condition_files/ini1993/light_data.bin') {
  ## Get file info
  finfo = file.info(file)
  
  ## read file
  par.file = file(description = file, 'rb')
  light = readBin(par.file, double(), size=4, n = finfo$size, endian="big")
  close(par.file)
  
  ## format data
  array(light, dim = c(length(model$XC), length(model$YC), length(light)/(length(model$XC) * length(model$YC))))
}


## Legacy code
load.nc = function(file) {
  file = ncdf4::nc_open(file)
  data = list()
  for (var in names(file$var)) {
    data[[var]] = ncdf4::ncvar_get(nc = file, varid = var)
  }
  ncdf4::nc_close(file)
  
  ##return
  data
}


calc.nemuro.advdiff = function(nemuro, d, mask) {
  
  dz = diff(nemuro$ZC)
  
  ## Set NA to zero s othey dont poison everything!
  nemuro$PON[is.na(nemuro$PON)] = 0
  nemuro$DON[is.na(nemuro$DON)] = 0
  nemuro$DIN[is.na(nemuro$DIN)] = 0
  nemuro$Org[is.na(nemuro$Org)] = 0
  
  # PON
  lx = c(1:length(nemuro$XC))
  llx = c(1, 1:length(nemuro$XC[-1]))
  ly = c(1:length(nemuro$YC))
  lly = c(1, 1:length(nemuro$YC[-1]))
  
  if (d == 1) {
    dd = d
  } else {
    dd = d-1
  }
  
  pon.u = 0.5 * (nemuro$PON[lx,ly,d] + nemuro$PON[lx,lly,d]) * 1e3 * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = umol N d-1
  pon.v = 0.5 * (nemuro$PON[lx,ly,d] + nemuro$PON[llx,ly,d]) * 1e3 * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = umol N d-1
  pon.w = 0.5 * (nemuro$PON[,,dd] + nemuro$PON[,,d]) * 1e3 * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = umol N d-1
  pon.s = nemuro$PON[,,d] * 1e3 * 15 * nemuro$domain$dy * nemuro$domain$dx * mask[,,d+1]
  # DON
  don.u = 0.5 * (nemuro$DON[lx,ly,d] + nemuro$DON[lx,lly,d]) * 1e3 * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = umol N d-1
  don.v = 0.5 * (nemuro$DON[lx,ly,d] + nemuro$DON[llx,ly,d]) * 1e3 * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = umol N d-1
  don.w = 0.5 * (nemuro$DON[,,dd] + nemuro$DON[,,d]) * 1e3 * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = umol N d-1
  # Org
  org.u = 0.5 * (nemuro$Org[lx,ly,d] + nemuro$Org[lx,lly,d]) * 1e3 * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = umol N d-1
  org.v = 0.5 * (nemuro$Org[lx,ly,d] + nemuro$Org[llx,ly,d]) * 1e3 * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = umol N d-1
  org.w = 0.5 * (nemuro$Org[,,dd] + nemuro$Org[,,d]) * 1e3 * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = umol N d-1
  # DIN
  din.u = 0.5 * (nemuro$DIN[lx,ly,d] + nemuro$DIN[lx,lly,d]) * 1e3 * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = umol N d-1
  din.v = 0.5 * (nemuro$DIN[lx,ly,d] + nemuro$DIN[llx,ly,d]) * 1e3 * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = umol N d-1
  din.w = 0.5 * (nemuro$DIN[,,dd] + nemuro$DIN[,,d]) * 1e3 * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = umol N d-1
  # water
  water.u = nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # m3 d-1
  water.v = nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # m3 d-1
  water.w = nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx
  
  list(PON = list(U = pon.u, V = pon.v, W = pon.w, S = pon.s),
       DON = list(U = don.u, V = don.v, W = don.w),
       Org = list(U = org.u, V = org.v, W = org.w),
       DIN = list(U = din.u, V = din.v, W = din.w),
       Water = list(U = water.u, V = water.v, W = water.w),
       Total = list(U = pon.u + don.u + din.u + org.u, V = don.v + pon.v + din.v + org.v, W = pon.w + don.w + din.w + org.w))
}