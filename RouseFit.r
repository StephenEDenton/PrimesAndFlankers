# Rouse Model to fit all conditions
# Stephen Denton, Nov. 2010
rm(list = ls())

# Set toggles
saving = F
fitting = T
restricted = F

# Set initial values for parameters
# Parameter matrices are set up as [ Short P-Short F, Short P-Long F ; 
#                                    Long P-Short F , Long P-Long F  ]
N = 20
a.alpha = 0.147
e.alpha = matrix(c(0.0791, 0.213), 2, 1)
a.phi = 0.125
e.phi = matrix(c( 0.0783, 0.106), 1, 2)  
if (restricted) { a.beta = 0.0652 } else {
a.beta = matrix(c( 0.0621, 0.0980, 0.0351, 0.0647 ), 2, 2) }
e.beta = a.beta
e.beta_is_a.beta = all(e.beta == a.beta)
a.gamma = 0.03
e.gamma = a.gamma
e.gamma_is_a.gamma = all(e.gamma == a.gamma )
rho.prime = 1
tau.prime = 1

# Construct a single vector of parameters that can be handed to the optimizer
params2fit = c("a.alpha","e.alpha","a.phi","e.phi","a.beta","a.gamma")
init.params = NULL
param.vec.len = NULL
params.fit.str = NULL
for (idx in 1:length(params2fit) ) {
    param.append = eval(parse(text=params2fit[idx]))
    init.params = c(init.params, param.append)
    param.vec.len = c(param.vec.len, length(param.append))
    if (length(param.append) == 1 ) {
      params.fit.str = c(params.fit.str, params2fit[idx])  
    } else {
      params.fit.str = c(params.fit.str, paste(params2fit[idx],1:length(param.append),sep=""))
      # paste(params2fit[idx],which(param.append>0, arr.ind=TRUE),sep="")
    } 
} 

# Manuscript paramters 
if (restricted) {
  # init.params = c(0.149, 0.0789, 0.213, 0.125, 0.0789, 0.108, 0.0652, 0.0322)
  init.params = c(0.1494925, 0.07890128, 0.2129554, 0.1247343, 0.07888727, 0.1076305, 0.0652048, 0.0321632) # NLL = 28699.34
} else {
  # init.params = c(0.145, 0.0793, 0.214, 0.124, 0.0778, 0.104, 0.0621, 0.098, 0.0351, 0.0647, 0.0275)
  init.params = c(0.1447892, 0.07928391, 0.2137624, 0.1238942, 0.07778932, 0.1044196, 0.06208081, 0.09797519, 0.03507875, 0.06468665, 0.02750193) # NLL = 28171.35 
}

# Load the data (contains meanAcc, sumCorr, and totTrls)
data.filename = "PrimingObsData.Rdata"
load(data.filename)
obs.acc = mean.acc
obs.n.corr = sum.corr
obs.N = tot.trls 
obs.n.inc = obs.N-obs.n.corr
pred = pred.acc
n.prDur = dim(obs.acc)[3]
n.flDur = dim(obs.acc)[4]

# Get dimension of all parameter matrices
a.alpha.dim = dim(as.matrix(a.alpha))
e.alpha.dim = dim(as.matrix(e.alpha))
a.phi.dim = dim(as.matrix(a.phi))
e.phi.dim = dim(as.matrix(e.phi))
a.beta.dim = dim(as.matrix(a.beta))
e.beta.dim = dim(as.matrix(e.beta))
a.gamma.dim = dim(as.matrix(a.gamma))
e.gamma.dim = dim(as.matrix(e.gamma))
rho.prime.dim = dim(as.matrix(rho.prime))
tau.prime.dim = dim(as.matrix(tau.prime))

# Define trinomial function
dtrinom = function(n, N, prb) {
    probVec = choose(N,n[,1])*choose(N-n[,1],n[,2]) * 
                  prb[1]^n[,1] * prb[2]^n[,2] * prb[3]^n[,3]
    return(probVec)
}

# Feature state counts n[,1]=n_off, n[,2]=n_on, n[,3]=n_dis 
### Can be put inside the below function to make it self contained (faster this way)
n.bi = cbind(0:N, N:0)
n.tri = as.matrix(expand.grid(0:N, 0:N, 0:N))
n.tri = n.tri[rowSums(n.tri) == N,]

# Analytic ROUSE model function
rouse.model = function( a.alpha, e.alpha, a.phi, e.phi, a.beta, e.beta, a.gamma, e.gamma, rho.prime, tau.prime, N ) { 
   
  # Calculate the probability of an off feature in an UNPRIMED, UNFLANKED FOIL 
  p_off.uuf = 1 - a.gamma #; p_on.uf = a.gamma
  # Calculate the prob of all possible on and off feature combinations for an UNPRIMED FOIL
  p_uuf = dbinom(n.bi[,1], N, p_off.uuf)
  
  # Do the same for an UNPRIMED, UNFLANKED TARGET
  p_off.uut = (1 - a.gamma) * (1 - a.beta) #; p_on.ut = 1 - p_off.ut
  p_uut = dbinom(n.bi[,1], N, p_off.uut)
  
  # PRIMED, UNFLANKED FOIL 
  p_on.puf = (1 - rho.prime) * a.gamma
  p_dis.puf = rho.prime * (1 - (1 - a.alpha)*(1 - a.gamma))
  p_off.puf = 1 - (p_on.puf + p_dis.puf)
  p_puf = dtrinom(n.tri, N, c(p_off.puf,p_on.puf,p_dis.puf) )
  
  # PRIMED, UNFLANKED TARGET
  p_on.put = (1 - rho.prime) * (1 - (1 - a.gamma)*(1 - a.beta))
  p_dis.put = rho.prime * (1 - (1 - a.alpha)*(1 - a.beta)*(1 - a.gamma))
  p_off.put = 1 - (p_on.put + p_dis.put)
  p_put = dtrinom(n.tri, N, c(p_off.put,p_on.put,p_dis.put) )
   
  # UNPRIMED, FLANKED FOIL 
  p_on.upf = (1 - tau.prime) * a.gamma
  p_dis.upf = tau.prime * (1 - (1 - a.phi)*(1 - a.gamma))
  p_off.upf = 1 - (p_on.upf + p_dis.upf)
  p_upf = dtrinom(n.tri, N, c(p_off.upf,p_on.upf,p_dis.upf) )
  
  # UNPRIMED, FLANKED TARGET
  p_on.upt = (1 - tau.prime) * (1 - (1 - a.gamma)*(1 - a.beta))
  p_dis.upt = tau.prime * (1 - (1 - a.phi)*(1 - a.beta)*(1 - a.gamma))
  p_off.upt = 1 - (p_on.upt + p_dis.upt)
  p_upt = dtrinom(n.tri, N, c(p_off.upt,p_on.upt,p_dis.upt) )
  
  # PRIMED, FLANKED FOIL 
  p_on.ppf = (1 - rho.prime) * (1 - tau.prime) * a.gamma
  p_dis.ppf = rho.prime * tau.prime * (1 - (1 - a.alpha)*(1 - a.phi)*(1 - a.gamma))
  p_off.ppf = 1 - (p_on.ppf + p_dis.ppf)
  p_ppf = dtrinom(n.tri, N, c(p_off.ppf,p_on.ppf,p_dis.ppf) )
  
  # PRIMED, FLANKED TARGET
  p_on.ppt = (1 - tau.prime) * (1 - (1 - a.gamma)*(1 - a.beta))
  p_dis.ppt = rho.prime * tau.prime * (1 - (1 - a.alpha)*(1 - a.phi)*(1 - a.beta)*(1 - a.gamma))
  p_off.ppt = 1 - (p_on.ppt + p_dis.ppt)
  p_ppt = dtrinom(n.tri, N, c(p_off.ppt,p_on.ppt,p_dis.ppt) )
  
  # Calculate the five unique feature likelihood ratios
  F_off = (1 - e.beta)
  F_on = (1 - (1 - e.gamma) * (1 - e.beta) ) / ( 1 - ( 1 - e.gamma) )
  F_d.pr = (1 - (1 - e.gamma) * (1 - e.beta) * (1 - e.alpha)) /
            (1 - (1 - e.gamma) * (1 - e.alpha))
  F_d.fl = (1 - (1 - e.gamma) * (1 - e.beta) * (1 - e.phi)) /
            (1 - (1 - e.gamma) * (1 - e.phi))
  F_d.b = (1 - (1 - e.gamma) * (1 - e.beta) * (1 - e.alpha) * (1 - e.phi)) /
            (1 - (1 - e.gamma) * (1 - e.alpha) * (1 - e.phi))
  
  # Calculate log word likelihood ratios for all combinations of feature types (on, off, dis)
  log.F_off = log(F_off)
  log.F_on  = log(F_on)
  log_W.bi = n.bi[,1] * log.F_off + n.bi[,2] * log.F_on
  log_W.tri.pr = n.tri[,1] * log.F_off + n.tri[,2] * log.F_on + n.tri[,3] * log(F_d.pr)
  log_W.tri.fl = n.tri[,1] * log.F_off + n.tri[,2] * log.F_on + n.tri[,3] * log(F_d.fl)
  log_W.tri.b = n.tri[,1] * log.F_off + n.tri[,2] * log.F_on + n.tri[,3] * log(F_d.b)
  
  # Combine possible log word likelihood ratios for two uprimed/unflanked words
  log_W.UU.UU = as.matrix(expand.grid(log_W.bi, log_W.bi))
  # Combine possible log word likelihood ratios for primed and unprimed words
  log_W.PU.UU = as.matrix(expand.grid(log_W.tri.pr, log_W.bi))
  log_W.UP.UU = as.matrix(expand.grid(log_W.tri.fl, log_W.bi))
  log_W.PP.UU = as.matrix(expand.grid(log_W.tri.b , log_W.bi))
  log_W.PU.UP = as.matrix(expand.grid(log_W.tri.pr, log_W.tri.fl))
  
  # Compute the joint probablity of all possible unprimed targets and foils 
  # Across all possible combinations of flanked and unflanked targets and foils 
  p_uut.X.p_uuf = c(outer(p_uut, p_uuf))
  p_corr.NPNF = sum( p_uut.X.p_uuf[ log_W.UU.UU[,1] > log_W.UU.UU[,2] ] ) + 
             0.5 * sum( p_uut.X.p_uuf[ log_W.UU.UU[,1] == log_W.UU.UU[,2] ] )
             
  p_upt.X.p_uuf = c(outer(p_upt, p_uuf))
  p_corr.NPTF = sum( p_upt.X.p_uuf[ log_W.UP.UU[,1] > log_W.UP.UU[,2] ] ) + 
             0.5 * sum( p_upt.X.p_uuf[ log_W.UP.UU[,1] == log_W.UP.UU[,2] ] )
             
  p_uut.X.p_upf = c(outer(p_upf, p_uut))
  p_corr.NPFF = sum( p_uut.X.p_upf[ log_W.UP.UU[,1] < log_W.UP.UU[,2] ] ) + 
             0.5 * sum( p_uut.X.p_upf[ log_W.UP.UU[,1] == log_W.UP.UU[,2] ] )           
             
  # Compute the joint probablity of all possible primed targets and unprimed foils
  # Across all possible combinations of flanked and unflanked targets and foils     
  p_put.X.p_uuf = c(outer(p_put, p_uuf))
  p_corr.TPNF = sum( p_put.X.p_uuf[ log_W.PU.UU[,1] > log_W.PU.UU[,2] ] ) + 
             0.5 * sum( p_put.X.p_uuf[ log_W.PU.UU[,1] == log_W.PU.UU[,2] ] )
             
  p_ppt.X.p_uuf = c(outer(p_ppt, p_uuf))
  p_corr.TPTF = sum( p_ppt.X.p_uuf[ log_W.PP.UU[,1] > log_W.PP.UU[,2] ] ) + 
             0.5 * sum( p_ppt.X.p_uuf[ log_W.PP.UU[,1] == log_W.PP.UU[,2] ] )
  
  p_put.X.p_upf = c(outer(p_put, p_upf))
  p_corr.TPFF = sum( p_put.X.p_upf[ log_W.PU.UP[,1] > log_W.PU.UP[,2] ] ) + 
             0.5 * sum( p_put.X.p_upf[ log_W.PU.UP[,1] == log_W.PU.UP[,2] ] )
             
  # Compute the joint probablity of all possible unprimed targets and primed foils
  # Across all possible combinations of flanked and unflanked targets and foils             
  p_uut.X.p_puf = c(outer(p_puf, p_uut))
  p_corr.FPNF = sum( p_uut.X.p_puf[ log_W.PU.UU[,1] < log_W.PU.UU[,2] ] ) + 
             0.5 * sum( p_uut.X.p_puf[ log_W.PU.UU[,1] == log_W.PU.UU[,2] ] )
  
  p_upt.X.p_puf = c(outer(p_puf, p_upt))
  p_corr.FPTF = sum( p_upt.X.p_puf[ log_W.PU.UP[,1] < log_W.PU.UP[,2] ] ) + 
             0.5 * sum( p_upt.X.p_puf[ log_W.PU.UP[,1] == log_W.PU.UP[,2] ] )
  
  p_uut.X.p_ppf = c(outer(p_ppf, p_uut))
  p_corr.FPFF = sum( p_uut.X.p_ppf[ log_W.PP.UU[,1] < log_W.PP.UU[,2] ] ) + 
             0.5 * sum( p_uut.X.p_ppf[ log_W.PP.UU[,1] == log_W.PP.UU[,2] ] )
  
  cond.pred = matrix( c(p_corr.NPNF, p_corr.NPTF, p_corr.NPFF,
                  p_corr.TPNF, p_corr.TPTF, p_corr.TPFF, 
                  p_corr.FPNF, p_corr.FPTF, p_corr.FPFF), 3, byrow=T)
  
  return(cond.pred)
}


# Function to run a fit 
rouse.fit = function( params ) {
  if ( any(params<0) | any(params>1) ) {
      return( NaN )
  }
  # Assign parameter values
  cat("\n***********************************************************************\n") 
  cat("<", paste(params.fit.str, collapse = ", "), ">\n< " )
  cat(params, sep=", ") ; cat(" >\n")
  param.idx = 0
  for ( i in 1:length(params2fit) ) {
    param.vec = NULL
    for ( idx in 1:param.vec.len[i] ) {
      param.idx = param.idx + 1
      param.vec = c( param.vec, params[param.idx] ) 
    }
    assign(params2fit[i], param.vec, envir = .GlobalEnv  )
  }
  if ( e.beta_is_a.beta ) { e.beta = a.beta }
  if ( e.gamma_is_a.gamma ) { e.gamma = a.gamma }  
  a.alpha.mat = matrix( a.alpha, n.prDur, n.flDur, byrow=(a.alpha.dim[1]==1) )
  e.alpha.mat = matrix( e.alpha, n.prDur, n.flDur, byrow=(e.alpha.dim[1]==1) )
  a.phi.mat   = matrix( a.phi  , n.prDur, n.flDur, byrow=(a.phi.dim[1]==1) )
  e.phi.mat   = matrix( e.phi  , n.prDur, n.flDur, byrow=(e.phi.dim[1]==1) )
  a.beta.mat  = matrix( a.beta , n.prDur, n.flDur, byrow=(a.beta.dim[1]==1) )
  e.beta.mat  = matrix( e.beta , n.prDur, n.flDur, byrow=(e.beta.dim[1]==1) )
  a.gamma.mat = matrix( a.gamma, n.prDur, n.flDur, byrow=(a.gamma.dim[1]==1) )
  e.gamma.mat = matrix( e.gamma, n.prDur, n.flDur, byrow=(e.gamma.dim[1]==1) )
  rho.mat     = matrix( rho.prime, n.prDur, n.flDur, byrow=(rho.prime.dim[1]==1) )
  tau.mat     = matrix( tau.prime, n.prDur, n.flDur, byrow=(tau.prime.dim[1]==1) )
  
  # Get predictions from ROUSE model
  for ( prDurIdx in 1:n.prDur ) {
    for ( flDurIdx in 1:n.flDur ) {
        pred[,,prDurIdx,flDurIdx]  <<- rouse.model( a.alpha.mat[prDurIdx,flDurIdx], 
                                                e.alpha.mat[prDurIdx,flDurIdx],
                                                a.phi.mat[prDurIdx,flDurIdx], 
                                                e.phi.mat[prDurIdx,flDurIdx], 
                                                a.beta.mat[prDurIdx,flDurIdx], 
                                                e.beta.mat[prDurIdx,flDurIdx], 
                                                a.gamma.mat[prDurIdx,flDurIdx], 
                                                e.gamma.mat[prDurIdx,flDurIdx], 
                                                rho.mat[prDurIdx,flDurIdx],
                                                tau.mat[prDurIdx,flDurIdx], N)
    }
  }
  
  # Print some output
  cat("Pred. Probs:\n")
  names(pred) = names(obs.acc)
  print(pred)
  
  # Calculate the negative log likelihood 
  NLL <<- -sum( obs.n.corr*log(pred) + (obs.n.inc)*log(1-pred) ) 
  
  pred.corr = obs.N * pred
  pred.inc  = obs.N * (1-pred)

  # Chi.Sqr <<- sum((obs.n.corr - pred.corr)^2/pred.corr) + sum((obs.n.inc - pred.inc)^2/pred.inc)
  G.Sqr <<- 2 * ( sum(obs.n.corr*log(obs.n.corr/pred.corr)) + sum(obs.n.inc*log(obs.n.inc/pred.inc)) )
  # return(G.Sqr)

  return(NLL)
}

if (fitting) {
  bestFit = optim(init.params, rouse.fit, control=list(trace=1))
  bestFitParams = bestFit$par
  
  cat("\nOptim output\n")
  print(bestFit)
} else {
  bestFitParams = init.params
}

# Use best fit parameters to run the model to get best fit predictions
cat("ROUSE Predictions")      
model.pred = rouse.fit(bestFitParams)
cat("Obs. Probs:\n")
print(obs.acc)
cat("NLL =", NLL, "\n\n")
cat("G-Sqr =", G.Sqr, "\n\n")
#cat("Chi-Sqr =", Chi.Sqr, "\n\n")

# Save predictions for plotting
pred.acc = pred
if (saving) { save( mean.acc, sum.corr, tot.trls, pred.acc, sem, subjNums, 
                    mean.acc.pr, sem.pr, mean.acc.fl, sem.fl, NLL,
                    file=data.filename) }
