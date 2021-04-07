library(TheSource)
library(ncdf4)
source('RX/source.r')
source('RX/NEMURO.source.r')

## Load in a NEMURO File
raw.files = list.files('Z:/Data/Taylor/gom_data/may_runs_nemuro_gom', pattern = '.nc', full.names = T, recursive = T)
raw.files = raw.files[!grepl('raw', raw.files)]
raw.files = raw.files[grepl('2002', raw.files)]
  
for (i in 1:length(raw.files)) {
    message()
    
    ## Datetime Info
    year = strsplit(raw.files[i], '/')[[1]][6]
    month = substr(strsplit(raw.files[i], '/')[[1]][7], 8,9)
    day = substr(strsplit(raw.files[i], '/')[[1]][7], 10,11)
    
    out.file = paste0('Z:/Data/Taylor/_rdata/NEMURO.', year, '.', month, '.', day, '.rdata')
    
    if (file.exists(out.file)) {
      message(Sys.time(), ': File exists for ', raw.files[i], '. Skipping.')
    } else {
      
      message(Sys.time(), ': Loading Data for ', month, '/', day, '/', year, ' (', i, '/', length(raw.files), ')')
      
      ## Load and build domain
      nemuro = load.nc(raw.files[i])
      nemuro$domain = list(Year = as.numeric(year),
                           Month = as.numeric(month),
                           Day = as.numeric(day),
                           XC = nemuro$XC,
                           YC = nemuro$YC,
                           ZC = nemuro$ZC)
      delta = as.numeric(difftime(conv.time.matlab(nemuro$TC)[2], conv.time.matlab(nemuro$TC)[1], units = 'days'))
      
      ## Load external files
      nemuro.light = load.lightfield(nemuro, paste0('Z:/Data/Taylor/gom_data/initial_condition_files/ini', year, '/light_data.bin'))
      
      physics = load.physics(nemuro,
                             U = paste0('Z:/Data/Taylor/gom_data/may_runs_nemuro_gom/', year, '/forcing_', nemuro$domain$Year, '/uvel/u.00', delta, '.data'),
                             V = paste0('Z:/Data/Taylor/gom_data/may_runs_nemuro_gom/', year, '/forcing_', nemuro$domain$Year, '/vvel/v.00', delta, '.data'),
                             W = paste0('Z:/Data/Taylor/gom_data/may_runs_nemuro_gom/', year, '/forcing_', nemuro$domain$Year, '/wvel/w.00', delta, '.data'),
                             S = paste0('Z:/Data/Taylor/gom_data/may_runs_nemuro_gom/', year, '/forcing_', nemuro$domain$Year, '/sal/s.00', delta, '.data'),
                             T = paste0('Z:/Data/Taylor/gom_data/may_runs_nemuro_gom/', year, '/forcing_', nemuro$domain$Year, '/tmp/t.00', delta, '.data'))
      ## TODO Load river
      river = load.river(dis.file = paste0('Z:/Data/Taylor/gom_data/initial_condition_files/ini', year, '/riv_dis_data.bin'),
                         nut.file = paste0('Z:/Data/Taylor/gom_data/initial_condition_files/ini', year, '/riv_nutr_data.bin'))
      nemuro$river = river[,,delta] ## Already transposed!
      
      message(Sys.time(), ': Compiling dataset.  ', appendLF = F)
      k = c(1:20)
      #message('Unstaggering U, V, W... ', appendLF = F)
      nemuro$U = physics$U[,,k]
      nemuro$V = physics$V[,,k]
      nemuro$W = physics$W[,,k]
      nemuro$S = physics$S[,,k]
      nemuro$T = physics$T[,,k]
      rm(physics)
      
      ## DIN
      nemuro$DIN = nemuro$NO[,,k] + nemuro$NH[,,k]
      nemuro$NH = NULL
      nemuro$NH = NULL
      nemuro$SI = NULL
      ## DON
      nemuro$DON = nemuro$DON[,,k]
      ## PON
      nemuro$PON = nemuro$PON[,,k]
      nemuro$OP = NULL
      ##Org
      nemuro$Org = nemuro$SP[,,k] + nemuro$LP[,,k] + nemuro$SZ[,,k] + nemuro$LZ[,,k] + nemuro$PZ[,,k]
      nemuro$SP = NULL
      nemuro$LP = NULL
      nemuro$SZ = NULL
      nemuro$LZ = NULL
      nemuro$PZ = NULL
      nemuro$CHL = NULL
      ## Conditions
      nemuro$S = NULL
      
      ## Domain
      nemuro$domain$grid = expand.grid(lon = nemuro$XC, lat = nemuro$YC)
      nemuro$ZC = nemuro$ZC[k]
      nemuro$domain$ZC = nemuro$domain$ZC[k]
      
      dx = matrix(0, nrow = length(nemuro$XC), ncol = length(nemuro$YC))
      dy = matrix(0, nrow = length(nemuro$XC), ncol = length(nemuro$YC))
      
      for (i in 1:length(nemuro$XC)) {
        for (j in 1:length(nemuro$YC)) {
          if (i == 1) {
            dx[i,j] = geosphere::distCosine(p1 = c(nemuro$XC[2], nemuro$YC[j]), p2 = c(nemuro$XC[1], nemuro$YC[j]))
          } else {
            dx[i,j] = geosphere::distCosine(p1 = c(nemuro$XC[i], nemuro$YC[j]), p2 = c(nemuro$XC[i-1], nemuro$YC[j]))
          }
          if (j == 1) {
            dy[i,j] = geosphere::distCosine(p1 = c(nemuro$XC[i], nemuro$YC[2]), p2 = c(nemuro$XC[i], nemuro$YC[1]))
          } else {
            dy[i,j] = geosphere::distCosine(p1 = c(nemuro$XC[i], nemuro$YC[j]), p2 = c(nemuro$XC[i], nemuro$YC[j-1]))
          }
        }
      }
      
      nemuro$domain$dx = dx
      nemuro$domain$dy = dy
      nemuro$Mask = !is.na(nemuro$DIN)
      
      lx = c(1:length(nemuro$XC))
      llx = c(1, 1:length(nemuro$XC[-1]))
      ly = c(1:length(nemuro$YC))
      lly = c(1, 1:length(nemuro$YC[-1]))
      dz = diff(nemuro$ZC)
      
      # PON
      nemuro$PON.Flux.U = nemuro$PON
      nemuro$PON.Flux.V = nemuro$PON
      nemuro$PON.Flux.W = nemuro$PON
      nemuro$PON.Flux.S = nemuro$PON
      # DON
      nemuro$DON.Flux.U = nemuro$PON
      nemuro$DON.Flux.V = nemuro$PON
      nemuro$DON.Flux.W = nemuro$PON
      #DIN
      nemuro$DIN.Flux.U = nemuro$PON
      nemuro$DIN.Flux.V = nemuro$PON
      nemuro$DIN.Flux.W = nemuro$PON
      # Org
      nemuro$Org.Flux.U = nemuro$PON
      nemuro$Org.Flux.V = nemuro$PON
      nemuro$Org.Flux.W = nemuro$PON
      # Water
      nemuro$Water.Flux.U = nemuro$PON
      nemuro$Water.Flux.V = nemuro$PON
      nemuro$Water.Flux.W = nemuro$PON
      
      for (d in 1:length(nemuro$ZC)) {
        if (d == 1) { dd = d } else { dd = d-1 }
        
        # PON
        nemuro$PON.Flux.U[,,d] = 0.5 * (nemuro$PON[lx,ly,d] + nemuro$PON[lx,lly,d]) * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = mmol N d-1
        nemuro$PON.Flux.V[,,d] = 0.5 * (nemuro$PON[lx,ly,d] + nemuro$PON[llx,ly,d]) * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = mmol N d-1
        nemuro$PON.Flux.W[,,d] = 0.5 * (nemuro$PON[,,dd] + nemuro$PON[,,d]) * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = mmol N d-1
        if (d == length(nemuro$ZC)) {
          nemuro$PON.Flux.S[,,d] = 0 ## No sinking out of last cell.
        } else {
          nemuro$PON.Flux.S[,,d] = nemuro$PON[,,d] * 15 * nemuro$domain$dy * nemuro$domain$dx * nemuro$Mask[,,d+1]
        }
        # DON
        nemuro$DON.Flux.U[,,d] = 0.5 * (nemuro$DON[lx,ly,d] + nemuro$DON[lx,lly,d])  * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = umol N d-1
        nemuro$DON.Flux.V[,,d] = 0.5 * (nemuro$DON[lx,ly,d] + nemuro$DON[llx,ly,d])  * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = umol N d-1
        nemuro$DON.Flux.W[,,d] = 0.5 * (nemuro$DON[,,dd] + nemuro$DON[,,d]) * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = umol N d-1
        # Org
        nemuro$Org.Flux.U[,,d] = 0.5 * (nemuro$Org[lx,ly,d] + nemuro$Org[lx,lly,d])  * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = umol N d-1
        nemuro$Org.Flux.V[,,d] = 0.5 * (nemuro$Org[lx,ly,d] + nemuro$Org[llx,ly,d])  * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = umol N d-1
        nemuro$Org.Flux.W[,,d] = 0.5 * (nemuro$Org[,,dd] + nemuro$Org[,,d]) * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = umol N d-1
        # DIN
        nemuro$DIN.Flux.U[,,d] = 0.5 * (nemuro$DIN[lx,ly,d] + nemuro$DIN[lx,lly,d])  * nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # C * vel * area = umol N d-1
        nemuro$DIN.Flux.V[,,d] = 0.5 * (nemuro$DIN[lx,ly,d] + nemuro$DIN[llx,ly,d])  * nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # C * vel * area = umol N d-1
        nemuro$DIN.Flux.W[,,d] = 0.5 * (nemuro$DIN[,,dd] + nemuro$DIN[,,d]) * nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx # C * vel * area = umol N d-1
        # water
        nemuro$Water.Flux.U[,,d] = nemuro$U[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dy[lx,ly] + nemuro$domain$dy[lx,lly]) # m3 d-1
        nemuro$Water.Flux.V[,,d] = nemuro$V[,,d] * 86400 * dz[d] * 0.5 * (nemuro$domain$dx[lx,ly] + nemuro$domain$dx[llx,ly]) # m3 d-1
        nemuro$Water.Flux.W[,,d] = nemuro$W[,,d] * 86400 * nemuro$domain$dy * nemuro$domain$dx
      }
      
      
      message(Sys.time(), ': Saving to file ', out.file)
      save(nemuro, file = out.file, compression_level = 1)
    }
}

