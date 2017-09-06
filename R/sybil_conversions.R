
S <- Matrix::sparseMatrix(i=c(1, 1, 2, 2, 3), j=c(1, 2, 3, 1, 2), x=c(1, 1, 2, 1, 3))
S <- Matrix::sparseMatrix(i=c(1, 1, 1, 2, 2), j=c(1, 2, 3, 1, 2), x=c(1, 2, 3, 1, 1))

S_to_GRB_CSC <- function(S) {
  list(
    beg = S@p[1:(length(S@p)-1)],
    len = S@p[2:length(S@p)] - S@p[1:(length(S@p)-1)],
    ind = S@i,
    val = S@x
  )
}

S_to_GRB_CSR <- function(S) {
  S_to_GRB_CSC(Matrix::t(S))
}

as_GRBmodel <- function(x) UseMethod("as_GRBmodel")
as_GRBmodel.modelorg <- function(sybil_model, add_S=TRUE, add_obj=TRUE,
                                 constrnames=sybil::met_id(sybil_model),
                                 varnames=sybil::react_id(sybil_model)) {
  model <- grb::GRBmodel$new(
    env=grb::GRBenv$new(),
    name=sybil::mod_name(sybil_model),
    numvars=sybil::react_num(sybil_model),
    obj=if (add_obj) sybil::obj_coef(sybil_model) else NULL,
    lb=sybil::lowbnd(sybil_model),
    ub=sybil::uppbnd(sybil_model),
    vnames=sybil::react_id(sybil_model)
  )

  model$set_model_sense(maximize=TRUE)

  if (add_S) {
    csr <- S_to_GRB_CSR(sybil::S(sybil_model))
    .Call(
      grb:::GRB_addconstrs,
      model$exptr,
      sybil::met_num(sybil_model),
      length(csr$val),
      csr$beg,
      csr$ind,
      csr$val,
      paste(rep('=', sybil::met_num(sybil_model)), collapse = ""),
      NULL,
      constrnames
    )
    model$update()
  }

  return(model)
}
