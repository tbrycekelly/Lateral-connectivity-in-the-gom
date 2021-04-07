
param = list(
  ## Small phytoplankton
  V.sp = 0.4, # d-1
  KNO2.sp = 0.5, # mmol N m-3
  KNH2.sp = 0.1, # mmol N m-3
  alpha.sp = 0.1, # m2 W-1 d-1
  beta.sp = 4.5e-4, # m2 W-1 d-1
  R.sp = 0.03, # d-1
  M.sp = 0.002, # d-1
  E.sp = 0.135, # ND
  QV.sp = 0.0693, # C-1
  QR.sp = 0.0693,
  QM.sp = 0.0693,
  ## Large Phytoplankton
  V.lp = 0.8, # d-1
  KNO2.lp = 3.0,
  KNH2.lp = 0.3,
  KSI.lp = 6.0,
  alpha.lp = 0.1,
  beta.lp = 4.5e-4,
  R.lp = 0.03,
  M.lp = 0.001,
  E.lp = 0.135,
  QV.lp = 0.0693,
  QR.lp = 0.0693,
  QM.lp = 0.0693,
  ## Small Zooplankton
  G.spsz = 0.6,
  phi.spsz = 1.4,
  tau.spsz = 0.04,
  AE.sz = 0.7,
  GGE.sz = 0.3,
  M.sz = 0.022,
  QG.sz = 0.0693,
  QM.sz = 0.0693,
  ## Large zooplankton
  G.splz = 0,
  G.lplz = 0.3,
  G.szlz = 0.3,
  phi.lplz = 1.4,
  phi.szlz = 1.4,
  tau.lplz = 0.04,
  tau.szlz = 0.04,
  AE.lz = 0.7,
  GGE.lz = 0.3,
  M.lz = 0.022,
  QG.lz = 0.0693,
  QM.lz = 0.0693,
  ## Predatory Zooplankton
  G.lppz = 0.1,
  G.szpz = 0.1,
  G.lzpz = 0.3,
  phi.lppz = 1.4,
  phi.szpz = 1.4,
  phi.lzpz = 1.4,
  tau.lppz = 0.04,
  tau.szpz = 0.04,
  tau.lzpz = 0.04,
  AE.pz = 0.7,
  GGE.pz = 0.3,
  M.pz = 0.12,
  QG.pz = 0.0693,
  QM.pz = 0.0693,
  lambda.pz = 4.605,
  lambda.lzpz = 3.01,
  ## Implicit Bacterial Rates
  B.nit = 0.003,
  B.ammon = 0.01,
  B.don = 0.05,
  B.dammon = 0.02,
  B.si = 0.01,
  Q = 0.0693,
  ## other Parameters
  e.w = 0.03,
  e.sp = 0.03,
  e.lp = 0.03,
  par = 0.43,
  Si.N = 1.0,
  Omega = 15,
  Si.Nriver = 2.0,
  ## Chl:C Sub Model
  ChlC.spmin = 0,
  ChlC.spmax = 0.015,
  ChlC.lpmin = 0.005,
  ChlC.lpmax = 0.03,
  alpha.chl = 0.28
)


state = list(
  SP = NULL,
  LP = NULL, 
  SZ = NULL,
  LZ = NULL,
  PZ = NULL,
  NO = NULL,
  NH = NULL,
  SI = NULL,
  DON = NULL,
  PON = NULL,
  OP = NULL,
  PAR = NULL,
  T = NULL
)


make.state = function(model, x, y, z, t) {
  
  list(
    SP = model$SP[x,y,z],
    LP = model$LP[x,y,z], 
    SZ = model$SZ[x,y,z],
    LZ = model$LZ[x,y,z],
    PZ = model$PZ[x,y,z],
    NO = model$NO[x,y,z],
    NH = model$NH[x,y,z],
    SI = model$SI[x,y,z],
    DON = model$DON[x,y,z],
    PON = model$PON[x,y,z],
    OP = model$OP[x,y,z],
    PAR = model$PAR[x,y,z],
    T = model$T[x,y,z]
  )
}


make.state.1d = function(model, x, y) {
  state = list()
  
  for (i in 1:20) {
    state[[i]] = make.state(model, x, y, i)
  }
  
  state
}


nemuro.1d = function(param, state, meta) {
  dz = c(diff(meta$ZC), 1)
  diagnostics = list(sinking = rep(0, length(state)))
  
  ## Loop for each time step
  for (i in 1:(meta$delta.t / meta$dt)) {
    meta2 = meta
    meta2$delta.t = meta$dt
    message('   ', i, ' time')
    ## Calculate light field
    #state[[1]]$PAR = state[[1]]$Surface.PAR * param$par
    #for (d in 2:length(state)) {
    #  state[[d]]$PAR = state[[d-1]]$PAR * exp(-(param$e.w + param$e.sp * state[[d-1]]$SP + param$e.lp * state[[d-1]]$LP) * dz[d-1]) 
    #}
    
    ## Run for each depth
    for (d in 1:length(state)) {
      ## Run biology time step
      state[[d]] = nemuro.0d(param, state[[d]], meta2)$state
      
      ## sinking
      if (d < length(state)) {
        d.pon = state[[d]]$PON * param$Omega * meta$dt / dz[d]
        state[[d]]$PON = max(0, state[[d]]$PON - d.pon)
        state[[d+1]]$PON = max(0, state[[d+1]]$PON + d.pon)
      }
      if (d == length(state)) {
        d.pon = 0 ## No sedimentation
      }
      diagnostics$sinking[d] = diagnostics$sinking[d] + d.pon
      
    }
  }
  
  state$diagnostics = diagnostics
  state
}


nemuro.0d = function(param, state, meta) {
  
  dt = meta$dt
  light = state$PAR 
  
  state2 = state
  diagnostic = list()
  
  ## Run model forward for delta.t
  for (i in 1:(meta$delta.t / meta$dt)) {
    
    ## Starting state for backwards euler stepping:
    state1 = state2
    
    
    ## Start backward euler time stepping:
    for (j in c(1:5)) {
      
      ## Small Phytoplankton uptake terms
      l.lim.sp = (1 - exp(-1 * param$alpha.sp * light / param$V.sp)) * exp(-1 * param$beta.sp * light / param$V.sp)
      t.lim.sp = exp(param$QV.sp * state$T) ## TODO * scale_tmp_phyto
      NO.lim.sp = state2$NO / (state2$NO + param$KNO2.sp) * (1 / (1 + param$KNH2.sp/param$KNO2.sp))
      NH.lim.sp = state2$NH / (state2$NH + param$KNH2.sp)
      # Deltas
      delta.no = state1$NO - (state1$NO / (1 + dt / state2$NO) * (param$V.sp * l.lim.sp * t.lim.sp * NO.lim.sp * state2$SP))
      delta.nh = state1$NH - (state1$NH / (1 + dt / state2$NH) * (param$V.sp * l.lim.sp * t.lim.sp * NH.lim.sp * state2$SP))
      delta.sp = param$E.sp * (delta.no + delta.nh)
      # Update
      state2$NO = state1$NO - delta.no
      state2$NH = state1$NH - delta.nh
      state2$DON = state1$DON + delta.sp
      state2$SP = state1$SP + delta.no + delta.nh - delta.sp
      # Diagnostics
      diagnostic$np.frac.sp = delta.no / (delta.no + delta.nh + 1e-12)
      diagnostic$npp.sp = delta.no + delta.nh
      
      ## Large Phytoplankton
      l.lim.lp = (1 - exp(-1 * param$alpha.lp * light / param$V.lp)) * exp(-1 * param$beta.lp * light / param$V.lp)
      t.lim.lp = exp(param$QV.lp * state$T) ## TODO * scale_tmp_phyto
      NO.lim.lp = state2$NO / (state2$NO + param$KNO2.lp) * (1 / (1 + param$KNH2.lp/param$KNO2.lp))
      NH.lim.lp = state2$NH / (state2$NH + param$KNH2.lp)
      Si.lim.lp = state2$SI / (state2$SI + param$KSI.lp)
      # Deltas
      delta.no = state2$NO - (state2$NO / (1 + dt / state2$NO) * (param$V.lp * l.lim.lp * t.lim.lp * NO.lim.lp * Si.lim.lp * state2$LP))
      delta.nh = state2$NH - (state2$NH / (1 + dt / state2$NH) * (param$V.lp * l.lim.lp * t.lim.lp * NH.lim.lp * Si.lim.lp * state2$LP))
      delta.sp = param$E.sp * (delta.no + delta.nh)
      # Update
      state2$NO = state2$NO - delta.no
      state2$NH = state2$NH - delta.nh
      state2$DON = state2$DON + delta.sp
      state2$LP = state1$LP + delta.no + delta.nh - delta.sp
      # Diagnostics
      diagnostic$np.frac.lp = delta.no / (delta.no + delta.nh + 1e-12)
      diagnostic$npp.lp = delta.no + delta.nh
      
      ## SP Respiration
      delta.sp = state2$SP - (state2$SP / (1 + dt / state2$SP * (param$R.sp * t.lim.sp * state2$SP)))
      state2$SP = state2$SP - delta.sp
      if (diagnostic$npp.sp < delta.sp) { ## If NPP < R then don't just make nitrate
        state2$NH = state2$NH + delta.sp
      } else { ## If NPP > R then respire nitrate and ammonium in proportion to uptake
        state2$NO = state2$NO + delta.sp * diagnostic$np.frac.sp
        state2$NH = state2$NH + delta.sp * (1 - diagnostic$np.frac.sp)
      }
      # Diagnostics
      diagnostic$npp.sp = diagnostic$npp.sp - delta.sp
      
      ## LP Respiration
      delta.lp = state2$LP - (state2$LP / (1 + dt / state2$LP * (param$R.lp * t.lim.lp * state2$LP)))
      state2$LP = state2$LP - delta.lp
      if (diagnostic$npp.lp < delta.lp) { ## Same for large phytoplankton
        state2$NH = state2$NH + delta.lp
      } else {
        state2$NO = state2$NO + delta.lp * diagnostic$np.frac.lp
        state2$NH = state2$NH + delta.lp * (1 - diagnostic$np.frac.lp)
      }
      # Diagnostics
      diagnostic$npp.lp = diagnostic$npp.lp - delta.lp
      
      ## Mortality
      delta.sp = state2$SP - (state2$SP / (1 + dt / state2$SP * (param$M.sp * t.lim.sp * state2$SP)))
      delta.lp = state2$LP - (state2$LP / (1 + dt / state2$LP * (param$M.lp * t.lim.lp * state2$LP)))
      # Update
      state2$SP = state2$SP - delta.sp
      state2$LP = state2$LP - delta.lp
      state2$PON = state2$PON + delta.sp + delta.lp
      state2$OP = state2$OP + delta.lp * param$Si.N
      
      ## Small Zooplankton
      t.lim.sz = exp(param$QG.sz * state$T)
      g.lim.sz.sp = max(c(0, 1 - exp(param$phi.spsz * (param$tau.spsz - state2$SP))))
      delta.sp = state2$SP - (state2$SP / (1 + dt / state2$SP * (param$G.spsz * t.lim.sz * g.lim.sz.sp * state2$SZ)))
      # Update
      state2$SP = state2$SP - delta.sp
      state2$SZ = state2$SZ + delta.sp * param$GGE.sz
      state2$NH = state2$NH + (param$AE.sz - param$GGE.sz) * delta.sp
      state2$PON = state2$PON + (1 - param$AE.sz) * delta.sp
      # Diagnostics
      diagnostic$sz.graz.sp = delta.sp * param$GGE.sz
      
      ## Large Zooplankton
      t.lim.lz = exp(param$QG.lz * state$T)
      g.lim.lz.lp = max(c(0, 1 - exp(param$phi.lplz * (param$tau.lplz - state2$LP))))
      g.lim.lz.sz = max(c(0, 1 - exp(param$phi.szlz * (param$tau.szlz - state2$SZ))))
      delta.lp = state2$SP - (state2$SP / (1 + dt / state2$SP * (param$G.lplz * t.lim.lp * g.lim.lz.lp * state2$LZ)))
      delta.sz = state2$SZ - (state2$SZ / (1 + dt / state2$SZ * (param$G.szlz * t.lim.sz * g.lim.lz.sz * state2$LZ)))
      # Update
      state2$LP = state2$LP - delta.lp
      state2$SZ = state2$SZ - delta.sz
      state2$LZ = state1$LZ + param$GGE.lz * (delta.lp + delta.sz)
      state2$NH = state2$NH + (param$AE.lz - param$GGE.lz) * (delta.lp + delta.sz)
      state2$PON = state2$PON + (1 - param$AE.lz) * (delta.lp + delta.sz)
      state2$OP = state2$OP + delta.lp * param$Si.N
      # Diagnostics
      diagnostic$lz.graz.lp = delta.lp * param$GGE.lz
      diagnostic$lz.graz.sz = delta.sz * param$GGE.sz ## Wrong in taylor's code
      
      ## Predatory Zooplankton
      t.lim.pz = exp(param$QG.pz * state$T)
      g.lim.pz.lp = max(c(0, (1 - exp(param$phi.lppz * (param$tau.lppz - state2$LP))) * exp(-param$lambda.pz * (state2$SZ + state2$LZ))))
      g.lim.pz.sz = max(c(0, (1 - exp(param$phi.szpz * (param$tau.szpz - state2$SZ))) * exp(-param$lambda.lzpz * state2$LZ)))
      g.lim.pz.lz = max(c(0, 1 - exp(param$phi.lzpz * (param$tau.lzpz - state2$LZ))))
      delta.lp = state2$LP - (state2$LP / (1 + dt / state2$LP * (param$G.lppz * t.lim.pz * g.lim.pz.lp * state2$PZ)))
      delta.sz = state2$SZ - (state2$SZ / (1 + dt / state2$SZ * (param$G.szpz * t.lim.pz * g.lim.pz.sz * state2$PZ)))
      delta.lz = state2$LZ - (state2$LZ / (1 + dt / state2$LZ * (param$G.lzpz * t.lim.pz * g.lim.pz.lz * state2$PZ)))
      # Update
      state2$LP = state2$LP - delta.lp
      state2$SZ = state2$SZ - delta.sz
      state2$LZ = state2$LZ - delta.lz
      state2$PZ = state2$PZ + param$GGE.pz * (delta.lp + delta.sz + delta.lz)
      state2$NH = state2$NH + (param$AE.pz - param$GGE.pz) * (delta.lp + delta.lz + delta.sz)
      state2$PON = state2$PON + (1 - param$AE.pz) * (delta.lp + delta.lz + delta.sz)
      state2$OP = state2$OP + delta.lp * param$Si.N
      # Diagnostics
      diagnostic$pz.graz.lp = delta.lp * param$GGE.pz
      diagnostic$pz.graz.sz = delta.sz * param$GGE.pz
      diagnostic$pz.graz.lz = delta.lz * param$GGE.pz
      
      ## Zooplankton Mortality
      delta.sz = state2$SZ - (state2$SZ / (1 + dt / state2$SZ * (param$M.sz * t.lim.sz * state2$SZ)))
      delta.lz = state2$LZ - (state2$LZ / (1 + dt / state2$LZ * (param$M.lz * t.lim.lz * state2$LZ)))
      delta.pz = state2$PZ - (state2$PZ / (1 + dt / state2$PZ * (param$M.p * t.lim.pz * state2$PZ^2)))
      # Update
      state2$SZ = state2$SZ - delta.sz
      state2$LZ = state2$LZ - delta.lz
      state2$PZ = state2$PZ - delta.pz
      state2$PON = state2$PON + delta.sz + delta.lz + delta.pz
      
      ## Chemical transforms
      # Nitrification
      t.lim.nit = exp(param$Q * state$T)
      delta.nh = state2$NH - (state2$NH / (1 + dt / state2$NH * (param$B.nit * t.lim.nit * state2$NH)))
      state2$NH = state2$NH - delta.nh
      state2$NO = state2$NO + delta.nh
      diagnostic$nitrification = delta.nh
      # Decomp of PON to NH
      delta.pon = state2$PON - (state2$PON / (1 + dt / state$PON * (t.lim.nit * param$B.ammon * state2$PON)))
      state2$PON = state2$PON - delta.pon
      state2$NH = state2$NH + delta.pon
      diagnostic$pon.2.nh = delta.pon
      # Decomp of DOn to NH
      delta.don = state2$DON - (state2$DON / (1 + dt / state$DON * (t.lim.nit * param$B.dammon * state2$DON)))
      state2$DON = state2$DON - delta.don
      state2$NH = state2$NH + delta.don
      diagnostic$don.2.nh = delta.don
      # Decomp of PON to DON
      delta.pon = state2$PON - (state2$PON / (1 + dt / state$PON * (t.lim.nit * param$B.don * state2$PON)))
      state2$PON = state2$PON - delta.pon
      state2$DON = state2$DON + delta.pon
      diagnostic$pon.to.don = delta.pon
      # Decomp of OP to SI
      delta.op = state2$OP - (state2$OP / (1 + dt / state$OP * (t.lim.nit * param$B.si * state2$OP)))
      state2$OP = state2$OP - delta.op
      state2$SI = state2$SI + delta.op
      diagnostic$op.2.si = delta.op
    }
  }
  list(state = state2, diagnostic = diagnostic)
}

meta = list(ZC = nemuro$ZC, dt = 1/1000, delta.t = 1)
profile = make.state.1d(nemuro, x = 100, y = 100)
out = nemuro.1d(param, profile, meta)
state = profile[[1]]
