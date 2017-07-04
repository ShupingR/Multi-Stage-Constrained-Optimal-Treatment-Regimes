#
eval_f <- function(tau) {
  return( list( "objective" = values(tau)[[1]],
                "gradient" = Null ))
}
# constraint functions
# inequalities
eval_g_ineq <- function(tau) {
  constr <- values(tau)[[2]] - 10
  return( list( "constraints"=constr, "jacobian"=Null ) )
}

# equalities
eval_g_eq <- function( x ) {
  constr <- Null
  grad   <- Null
  return( list( "constraints"=constr, "jacobian"=grad ) )
}

# initial values
x0 <- c( 1, 5, 5, 1, 2, 4 )
# lower and upper bounds of control
lb <- c( 1, 1, 1, 1, 1, 1 )
ub <- c( 5, 5, 5, 5, 5, 5 )
local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel"  = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel"  = 1.0e-7,
              "maxeval"   = 1000,
              "local_opts" = local_opts)

res <- nloptr( x0=x0,
               eval_f=eval_f,
               lb=lb,
               ub=ub,
               eval_g_ineq=eval_g_ineq,
               eval_g_eq=eval_g_eq,
               opts=opts )
