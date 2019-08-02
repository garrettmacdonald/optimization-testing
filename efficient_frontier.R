# libraries ---------------------------------------------------------------
req_packages <- 
  c("quantmod", # to get return data from ETFs 
    "lpSolve",  # to get max return portfolio
    "quadprog") # to get min vol portfolio and frontier points

new_packages <- 
  req_packages[!(req_packages %in% installed.packages()[, "Package"])]

if(length(new_packages)) install.packages(new_packages)

sapply(req_packages, library, c = T)

# get etf returns ---------------------------------------------------------
tickers <- c("IVV", # US stocks
             "EFA", # DM stocks
             "LQD", # investment grade corporate bonds
             "SHY", # 1-3y UST
             "IEF", # 7-10y UST
             "TLT") # 20+y UST

daily_px <- NULL
for(i in tickers) {
  df <- getSymbols(i, auto.assign = FALSE)
  close <- df[as.character(2003:2018), 4]
  daily_px <- cbind(daily_px, close)
}

daily_ret <- log(as.matrix(daily_px[2:nrow(daily_px), ])) - 
  log(as.matrix(daily_px[1:(nrow(daily_px) - 1), ]))

# optimization parameters -------------------------------------------------
num_assets <- ncol(daily_ret)
num_group_constraints <- 1 # must be between 30-70% equities
num_periods <- nrow(daily_ret)
period_type <- 250 # for annualization

# for bvec
asset_class_min <- rep(0, num_assets) # 0-100% in each portfolio
asset_class_max <- rep(1, num_assets)
group_constraint_min <- rep(.3, num_group_constraints) # 30%-70% in equities
group_constraint_max <- rep(.7, num_group_constraints)

# for Amat
group_constraint_amat <- matrix(data = c(1, 1, 0, 0, 0, 0)) # equity assets

# optimization set-up -----------------------------------------------------
# bvec
bvec <- c(1, # must add up to 100%
          -asset_class_max,
          asset_class_min)

if (num_group_constraints > 0) {
  bvec <- c(bvec,
            -group_constraint_max,
            group_constraint_min)
}

# Amat
Amat <- cbind(rep(1, num_assets), diag(-1, num_assets), diag(1, num_assets))

if (num_group_constraints > 0) {
  Amat <- cbind(Amat,
                -group_constraint_amat,
                group_constraint_amat)
}

# dvec
dvec <- apply(daily_ret, 2, mean) * period_type

# Dmat
Dmat <- cov(daily_ret)

# optimization ------------------------------------------------------------
# solve for minimum volatility portfolio
mv_dvec <- dvec * 0

solution <- solve.QP(Dmat, mv_dvec, Amat, bvec, meq = 1)$solution
names(solution) <- colnames(daily_ret)

port_vol <- sqrt((t(solution) %*% Dmat %*% solution)) * sqrt(period_type)
port_ret <- t(solution) %*% dvec

sol_summary <- c(solution, 0, port_vol, port_ret)

# solve for maximum return portfolio
max_sol <- 
  lp("max", dvec, t(Amat), c("=", rep(">=", ncol(Amat) - 1)), bvec)$solution

max_vol <- sqrt((t(max_sol) %*% Dmat %*% max_sol)) * sqrt(period_type)
max_ret <- dvec %*% max_sol

# solve for incremental portfolios
num_sol <- 1000
tgt_ret <- port_ret
Amat_incr <- cbind(dvec, Amat)
range <- max_ret - port_ret
increment <- range / num_sol

opt_port_table <- matrix(NA, nrow = num_sol + 1, ncol = length(sol_summary))
opt_port_table[1, ] <- sol_summary

for (i in 1:num_sol) {
  tgt_ret <- tgt_ret + increment
  bvec_incr <- c(tgt_ret, bvec)
  solution <- 
    try(solve.QP(Dmat, dvec, Amat_incr, bvec_incr, meq = 2)$solution, TRUE)
  
  if (length(attr(solution, "class")) != 1) {
    solution <- as.array(solution)
    port_vol <- sqrt((t(solution) %*% Dmat %*% solution)) * sqrt(period_type)
    port_ret <- t(solution) %*% dvec
    opt_port_table[i + 1, ] <- c(solution, i, port_vol, port_ret)
  }
}

opt_port_table[num_sol + 1, ] <- c(max_sol, num_sol, max_vol, max_ret)

# clean up and plot
colnames(opt_port_table) <- c(colnames(daily_ret), "obs_num", "vol", "ret")
opt_port_table <- round(opt_port_table, 4)
opt_port_table <- as.matrix(opt_port_table, 4)
plot(opt_port_table[, "vol"], opt_port_table[, "ret"])