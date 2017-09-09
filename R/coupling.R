
near <- function(x, y, tol=1e-5) {
  abs(x - y) <= tol
}

clear_objective <- function(model) {
  model$setattr("Obj", setNames(numeric(model$get_sizes()$NumVars), model$get_names()$VarName))
}

flux_variability <- function(model, obj_frac=0.999, return_fluxes=FALSE) {
  prev_obj <- model$getattr("Obj")
  prev_sense <- model$getattr("ModelSense")
  if (!is.na(obj_frac) && !near(obj_frac, 0)) {
    model$optimize()
    obj_val <- model$getattr('ObjVal')
    # add the objective
    obj_con_name <- "OPT_CON"
    model$addconstr(model$getattr('Obj'), '>', obj_frac*obj_val, obj_con_name)
  } else {
    obj_con_name <- NA
  }
  n <- model$get_sizes()$NumVars
  maxflux <- numeric(n)
  minflux <- numeric(n)
  if (return_fluxes) {
    fluxes <- matrix(0.0, nrow=2*n, ncol=n)
  }
  vars <- model$get_names()$VarName
  for (i in 1:n) {
    print(i)
    clear_objective(model)
    model$setattr("Obj", setNames(1, vars[[i]]))
    model$set_model_sense(maximize=TRUE)
    model$optimize()
    sol <- model$get_solution()
    maxflux[i] <- sol$ObjVal
    if (return_fluxes) {
      fluxes[i, ] <- sol$X
    }
  }
  for (i in 1:n) {
    clear_objective(model)
    model$setattr("Obj", setNames(1, vars[[i]]))
    model$set_model_sense(minimize=TRUE)
    model$optimize()
    sol <- model$get_solution()
    minflux[i] <- sol$ObjVal
    if (return_fluxes) {
      fluxes[i+n, ] <- sol$X
    }
  }

  if (!is.na(obj_con_name)) {
    .Call("GRB_delconstrs", model$exptr, 1L, model$get_sizes()$NumConstrs - 1L)
  }
  model$setattr("Obj", prev_obj)
  model$setattr("ModelSense", prev_sense)

  retval <- list(
    minflux = minflux,
    maxflux = maxflux
  )
  if (return_fluxes) {
    retval$fluxes <- fluxes
  }

  return(retval)
}

flux_coupling <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01) {
  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName
  prev_obj <- model$getattr("Obj")
  model$setattr("Obj", setNames(numeric(n), vars)) # clear the objective
  prev_sense <- model$getattr("ModelSense")

  fva <- flux_variability(model, obj_frac=NA, return_fluxes=TRUE)
  use_min_fva_cor <- min_fva_cor != 0.0
  if (use_min_fva_cor) {
    fva_cor <- cor(fva$fluxes)
  }

  coupled <- matrix(FALSE, nrow=n, ncol=n, dimnames=list(vars, vars))
  fixed <- near(fva$minflux, fva$maxflux)
  blocked <- near(fva$minflux, 0) & near(fva$maxflux, 0)
  active <- !(fixed | blocked)

  not_fixed <- function(x,y) {
    !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
  }

  lp_calls <- 0
  for (i in 1:(n-1)) {
    if (!active[i]) next
    sub_max <- rep(-Inf, n)
    sub_min <- rep(Inf, n)
    prev_ub <- model$getattr("UB")[vars[i]]
    prev_lb <- model$getattr("LB")[vars[i]]
    fixed_val <- fva$minflux[i] + fix_frac*(fva$maxflux[i] - fva$minflux[i])
    model$setattr("UB", setNames(fixed_val + 0.5*fix_tol_frac*abs(fixed_val), vars[i]))
    model$setattr("LB", setNames(fixed_val - 0.5*fix_tol_frac*abs(fixed_val), vars[i]))
    for (j in (i+1):n) {
      if (!active[j]) next
      if (use_min_fva_cor && !is.na(fva_cor[i,j]) && abs(fva_cor[i,j]) < min_fva_cor) next
      if (not_fixed(sub_max[j], sub_min[j])) next

      skip <- FALSE
      model$setattr("Obj", setNames(1.0, vars[j]))
      if (!skip) {
        model$set_model_sense(maximize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)
        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (!skip) {
        model$set_model_sense(minimize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)
        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (!skip) {
        coupled[i,j] <- TRUE
        active[j] <- FALSE
      }

      model$setattr("Obj", setNames(0.0, vars[j]))
    }
    # unfix i
    model$setattr("UB", prev_ub)
    model$setattr("LB", prev_lb)
  }

  model$setattr("Obj", prev_obj)
  model$setattr("ModelSense", prev_sense)

  list(
    coupled = coupled,
    lp_calls = lp_calls
  )
}
