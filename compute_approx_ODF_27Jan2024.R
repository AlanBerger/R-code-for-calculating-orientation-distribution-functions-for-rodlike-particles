

compute_approx_ODF <- function(N.theta.intervals, N.phi.intervals, W, 
                                                       V = function(theta) {0}) { 
# December 8, 2023  Alan E. Berger    edits January 27, 2024

print("run of   compute_approx_ODF   January 27, 2024 version")

# Jan 27, 2024 allow an applied field V(theta) 
#                    V is to be a scalar function (only takes a single value as its argument)
#                    the default for V is  V(theta) = 0 
#                    if you don't want to include an applied field
#                    just leave V out of the argument list when call compute_approx_ODF 
# Jan 27, 2024 protect against 0 divide from complete convergence
#              before the minimum number of iterations

# Use the iterative method indicated by the calculus of variations, discretized 
# as in J. Herzfeld, A. E. Berger, J. W. Wingate, 
# A highly convergent algorithm for computing the orientation distribution functions 
# of rodlike particles, Macromolecules 1984, v17, 1718-1723
# [HBW 1984], 
# to calculate approximations to Orientation Distribution Functions 

# abbreviation used below: ODF is orientation distribution function

# The code below is generally using the notation and setup in the 
# Judith Herzfeld, Alan E. Berger and John W. Wingate paper, except, 
# doing calculations on [0, pi] 
# Also, rather than use Richardson extrapolation as done in
# [HBW 1984] to obtain increased accuracy for quantities of interest, here the 
# intent is to simply use larger numbers of grid points to obtain desired accuracy. 

# Use the  lexical scoping  and  function closure  properties of the 
# R programming language to do, once, various setup calculations, and then have the 
# results available. And define and return the function   calculate_ODF   that for 
# a given value of the constant B in the free energy function F, and initial 
# vector of values initial.f at the  grid points {Pj = j * pi / N.theta.intervals}, 
# j = 1, 2, ..., N.theta.intervals - 1), and function V(theta),
# calculates an approximate ODF f (a vector of values {fj} approximating the values
# of an ODF at the grid points {Pj}). 
# f generally should approximate the values of a local minimum of F at 
# the grid points {Pj}. Convergence of the iteration is tested as in [HBW 1984], 
# using as the convergence bound   tau.conv   (an argument of the 
# function  calculate_ODF  defined below).

# Note also: R. F. Kayser Jr. and H. J. Raveche', 
# Bifurcation in Onsager's model of the isotropic-nematic transition,
# Physical Review A 1978, v17, 2067-2072 
# [KR 1978]

# and 

# A. E. Berger, Analysis of a constrained minimization problem modeling the 
# orientation distribution of rod-like particles, Nonlinear Analysis, Theory, 
# Methods & Applications 1987, v11, 719-731 
# [B 1987].

# For the hard core reference system, the function W in the free energy function
# (W is a function whose argument is denoted by gamma) equals sin(gamma)
# gamma will be equal
#     acos(sin(theta1)*sin(theta2)*cos(phi) + cos(theta1)*cos(theta2))
# (see for example [KR 1978], [B 1987])

# Trapezoidal numerical integration is used to approximate the value of 
# integrals, and the presence of sin(theta) factors in the free energy function 
# results in the values of f at the points corresponding to theta = 0 and 
# theta = pi having no effect on the value of the discretized free energy. 
# Note f (the vector of values defined at the grid points Pj) has a natural
# extension to a function defined for all theta in [0, pi]: see equation (2.8) 
# in [B 1987]; and for the hard core reference system, W(gamma) equal sin(gamma)
# and V(theta) = 0, by Theorem 1.2 or Lemma 3.1 in [B 1987] any ODF that is a  
# local minimum of the actual (not discretized) free energy function F will 
# have (at least) 2 continuous derivatives 
# on [0, pi] and its first derivative will be 0 at theta = 0, pi/2 and pi

# Note for the free energy functions being considered here, if W(pi - gamma) = W(gamma)
# and V(pi - theta) = V(theta), local minima will be 
# symmetric about pi/2 

# In a number of places, we choose simpler code over faster running code.

# An approximation for the value f(0) is obtained as the value at theta = 0 of   
# the quadratic function that has the values f1, f2, f3 at the grid points P1,  
# P2, P3, respectively (as done in [HBW 1984]), and similarly for obtaining f(pi) 

# While values of an approximate (discretized) ODF at 
# theta = 0 and theta = pi do not matter for calculating integrals involving
# a sin(theta) factor via the trapezoidal formula, these values are used in 
# doing plots.


#23456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890
# try to keep to max of 80 characters per line
                                                                                

# Calculate various vectors and the matrix for the discretized version of the  
# integral operator in the free energy. These only need to be calculated once, 
# and the setup here ensures that the values are available and preserved between 
# calls to   calculate_ODF 

delta.theta <- pi / N.theta.intervals  # length of theta subintervals

# The Pj points (located in the vector theta.vec) are interior to [0, pi])
# These are the theta locations for the values in discrete ODF vectors.
# The sin(theta) factor in theta integrals results in values at theta
# equal 0 and theta equal pi having no effect on trapezoidal approximation
# values

num.theta.points <- N.theta.intervals - 1
theta.indices <- 1:num.theta.points
theta.vec <- theta.indices * delta.theta  # theta values for numerical 
#                                           integration

# Calculating values of the kernel of the integral operator in F 
# requires integrals over [0, 2*pi] which are calculated using 
# trapezoidal numerical integration. These integrals do involve values
# at phi = 0 and phi = 2*pi since there is no sin(phi) factor in the 
# integrals over phi

delta.phi <- 2*pi / N.phi.intervals  # length of phi subintervals

# phi example   4 phi intervals   N.phi.intervals equals 4

# phi points     c(0,  1,    2,              3,               4) * delta.phi  
#              2*pi equals N.phi.intervals * delta.phi 
# phi indices      1   2     3               4                5 
# 
# point location   0   1*delta.phi ... (N.phi.intervals - 1)*delta.phi  
#                                                   N.phi.intervals*delta.phi   

num.phi.points <- N.phi.intervals + 1
phi.indices <- 1:num.phi.points  # phi.indices[2:(num.phi.points - 1)] are
#           the phi indices for phi values in the interior of [0, 2*pi]
#           the interior phi points have trapezoidal weight 1
#           the phi points at 0 and 2*pi have trapezoidal weight 1/2
phi.vec <- (phi.indices - 1) * delta.phi  # locations of the phi points

sin.theta.vec <- sin(theta.vec)
cos.theta.vec <- cos(theta.vec)

# get the vector Vvec of the values from V evaluated at the theta.vec points
# V is expected to be a scalar function (argument is to be a single value)
Vvec <- unname(sapply(X = theta.vec, FUN = V))
max_mag_Vvec <- max(abs(Vvec))


# sin.phi.vec <- sin(phi.vec)  # not needed
cos.phi.vec <- cos(phi.vec)


# for theta1 defined as theta.vec[i] and 
# theta2 defined as theta.vec[j], i and j running
# from 1 through num.theta.points, we want to have the matrix of values K[i,j]
# defined by
# (1/(2*pi) integral from phi = 0 to 2*pi of W(gamma) 
# where gamma = acos(sin(theta1)*sin(theta2)*cos(phi) + cos(theta1)*cos(theta2))
# This integral will be calculated using trapezoidal numerical integration
# over phi

K <- matrix(0, nrow = num.theta.points, ncol = num.theta.points) # initialize K

# K will be symmetric, but it will be conceptually simpler just to calculate all 
# of the entries (this is only done once for a given run of B values).
# Having the full matrix K, and K being symmetric, means we can use column j
# of K for row j of K (and column entries have consecutive storage locations).

for (j in theta.indices) {
   for (i in theta.indices) {

      acos.arg.vec <- sin.theta.vec[i] * sin.theta.vec[j] * cos.phi.vec +
                           cos.theta.vec[i] * cos.theta.vec[j]

#          need to check that each entry of acos.arg.vec is in [0,1]
#          this could fail to be the case due to finite precision
#          for example, in one case got a couple values v equal to 
#          1.0000000000000002
#          which would result in acos(v) giving not a number (NaN)

      wm1 <- which(acos.arg.vec < -1)
      if(length(wm1) > 0) acos.arg.vec[wm1] <-  -1  

      w1 <- which(acos.arg.vec > 1)
      if(length(w1) > 0)  acos.arg.vec[w1] <-   1

      trapezoid.values <- W(acos(acos.arg.vec))  # W needs to handle vectors
#            (return a vector of values when called with a vector of values)

      trap.sum <- (trapezoid.values[1] + trapezoid.values[num.phi.points]) / 2 +
                        sum(trapezoid.values[2:(num.phi.points - 1)])

      K[i, j] <- trap.sum * delta.phi / (2 * pi)  
#             the denominator two * pi is in the definition of K
   }
}

# Now have finished the preliminary calculations, 
# next define the function  calculate_ODF  which calculates
# an ODF, given various required input values.
# The quantities that were calculated above are available
# to calculate_ODF since compute_approx_ODF returns 
# calculate_ODF in the form of a function closure
# (this is a property of the R programming language)


#23456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890
######################################################################


calculate_ODF <- function(initial.f, B, tau.conv, max.num.iter, min.iterations) {
# December 8, 2023  Alan E. Berger    edits January 27, 2024
# Jan 27, 2024 allow an applied field V(theta) 
#                    V is to be a scalar function (only takes a single value as its argument)
#                    the default for V is  V(theta) = 0 
#                    if you don't want to include an applied field
#                    just leave V out of the argument list when call compute_approx_ODF 
# Jan 27, 2024 protect against 0 divide from complete convergence
#              before the minimum number of iterations

print("run of    calculate_ODF    January 27, 2023  version")


# calculate what should be an approximation f to a local minimum of the
# free energy F with parameter value B, given the vector of initial
# values initial.f (corresponding to values of an ODF at the Pj points).
# The vector f gives the (approximate) values at the Pj points.
# 
# tau.conv is the convergence bound stopping condition as in [HBW 1984]
# max.num.iter is the maximum number of iterations allowed - a warning message
# will be printed if the iteration didn't converge (didn't satisfy
# the stopping condition within max.num.iter iterations)
# min.iterations    require having done  min.iterations  iterations before
#                   start testing for convergence           

# use the iterations suggested by the calculus of variations - see [HBW 1984]

# W(gamma) is the function used in calculating the matrix K which is the kernel
# of the discrete integral operator in the discrete form of F (K is obtained in
# compute_approx_ODF  above).   
# For the hard core reference system, W(gamma) = sin(gamma)

# gamma will be equal
#     acos(sin(theta1)*sin(theta2)*cos(phi) + cos(theta1)*cos(theta2))

# use the various calculated values preserved in the 
# function closure created by  compute_approx_ODF 

# will check on whether each iteration decreases 
# the (discretized / numerically calculated) free energy.
 
# When W satisfies appropriate conditions, it has been proven (AEB) that 
# each iteration, in the case of continuous ODFs,
# will decrease the free energy
# unless one was at a fixed point of the iteration 
# (manuscript giving the proof is in preparation).


# Return a list of desired quantities including the discrete converged 
# ODF f (so can plot it if desired) and various values associated with it,
# including f(0) (calculated using the quadratic polynomial 
# passing through f[1], f[2], f[3]), and f(pi) similarly calculated
# 
# an estimate of the local Lipschitz constant L,
# and calculated values for  pi/4 - <<W(gamma)>>  (see [HBW 1984] and the
# calculations below for how this and similar quantities are defined);
# <ln(f(theta))>, <P2(cos(theta))>, <P4(cos(theta))> where here P2 and P4 here
# are the second and fourth Legendre polynomial.
# Put the scalar quantities into a vector that will be returned in the list.

# and also return the calculated values of the free energy for the initial value f0,
# and for all the ODF except for the ODF from last iteration 
# (the way the calulation is set up, to simplify the code, at each iteration after 
#  the first one,one calculates the free energy for the previous ODF) 
# If the iteration converges, there shouldn't be much difference between the final
# ODF and the next to last one.
# 

# Also return the number of iterations it took to converge, and the number
# of iterative steps for which the calculated free energy failed to decrease

# first convert initial.f to be a (discrete) ODF that is positive
# (set any values <= 0 to be slightly above 0, and then scale (normalize) so 
# the integral over theta in [0, pi] of  (1/2) * f(theta) * sin(theta) = 1
# as calculated using trapezoidal numerical integration)

#23456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890

fn <- pmax(initial.f, 0.0001) # fn is now > 0
# having the initial fn > 0 avoids any issue with ln(f) since by the nature of
# the iteration, all successive iterates will also be > 0

theta.integral.of.fn <- sum(fn * sin.theta.vec) * delta.theta / 2  
# this is trapezoidal integration for 
# (1/2) integral over theta in [0, pi] of fn * sin(theta)

fn <- fn / theta.integral.of.fn  #  fn is now normalized
# so (1/2) integral over theta in [0, pi] of fn(theta) sin(theta) equals 1
# (as calculated by the trapezoidal formula)


######### do the first iteration "by hand" since need three 
######### successive approximate f values
######### to estimate the approximation error - see [HBW 1984]

# use the iteration suggested by the calculus of variations - see [HBW 1984]

# Also, calculate the value of the free energy F for fn

free.energy.for.each.iteration <- numeric(0)  
# successively append each calculated F, use for 
# plotting F vs. interation number

fnp1 <- numeric(num.theta.points)  
# initialize "container" vector for the next iterate f^{n+1}                                   

vector.for.free.energy <- numeric(num.theta.points)  # container for vector
#                           used for calculating F(fn) and for getting fnp1

 
for (i in theta.indices) {  # calculate each entry of fnp1

#  use trapezoidal formula to approximate 
#  integral over theta2 in [0, pi] of 
#                                      K(theta1, theta2) fn(theta2) sin(theta2)
   Ki <- K[, i]  # since K[i,j] is symmetric, this equals the vector 
#                  equal to the ith row of K 
#                  (K[i, ] whose entries corr. to values of theta2)
   trapez.sum <- sum(Ki * fn * sin.theta.vec) 
   value.of.integral.over.theta2 <- trapez.sum * delta.theta
   value.for.free.energy <-  B * value.of.integral.over.theta2 / 2
   vector.for.free.energy[i] <- value.for.free.energy  # for getting F for fn
   fnp1[i] <- exp(-Vvec[i] - value.for.free.energy)  # fnp1 is not yet normalized
}

# Deal with the limits of finite precision in computations:
# if the argument v of the exponential function is too large
# (this is above 709 and below 710 on my computer), exp(v) will, in R, 
# return Inf (a reserved word in R).
# This would lead to error messages from downstream calculations.
# fnp1 is equal exp(-Vvec - vector.for.free.energy); note the minus signs
   maxv <- max(abs(vector.for.free.energy)) + max_mag_Vvec
   if(maxv > 700) cat("an argument of exp will be ", maxv, "  at iteration ", 1)
   if(maxv > 700) {
      stop("have situation of abs of argument of exp > 700; ending execution")
   }


#  now normalize fnp1
   theta.integral.of.fnp1 <- sum(fnp1 * sin.theta.vec) * delta.theta / 2  
   fnp1 <- fnp1 / theta.integral.of.fnp1  #  fnp1 is now normalized

#  finish calculation of integral operator K contribution to the discrete
#  free energy (for fn) 
   integral.operator.K.part.of.free.energy <- 
        sum(fn * vector.for.free.energy * sin.theta.vec) * delta.theta / 4  

#  calculate the contribution to the discrete free energy from 
#  the natural log term (for fn)

   log.part.for.F <- sum(fn * log(fn) * sin.theta.vec) * delta.theta / 2  

  V.part.for.F <- sum(fn * Vvec * sin.theta.vec) * delta.theta / 2  

   free.energy.for.fn <- 
               log.part.for.F + V.part.for.F + integral.operator.K.part.of.free.energy 
   free.energy.for.each.iteration  <-
               c(free.energy.for.each.iteration, free.energy.for.fn)

# end of first iteration (done individually)

    
# Now can do successive iterations, estimating the Lipschitz constant L for the
# iteration and the approximation error (relative to max(fn)) as in [HBW 1984].
# Use tau.conv as the convergence criterion on the approximation error.

num_iter_F_did.NOT.decrease <- 0

for (iter in 2:max.num.iter) {
   iteration.num <- iter

   fnm1 <- fn  # fnm1 is the variable name for f^{n-1} 
#                when iter is 2, fn is from the calculations above
   free.energy.for.fnm1 <- free.energy.for.fn  
#       when iter is 2, free.energy.for.fn   is from the calculations above
   fn   <- fnp1  # when iter is 2, fnp1 is from the calculations above

#  calculate next iteration value fnp1 using fn, and then use fnm1, fn, fnp1 
#  to estimate the error between fnp1 and what the iteration 
#  supposedly converges to

# calculate the free energy for fn (while calculating fnp1) and test if
# it is < the free energy for fnm1


#  from previous calculations, fnp1 is an existing vector of the correct  
#  length so can use it as a container vector for the next iteration values

#23456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890

   for (i in theta.indices) {  # calculate each entry of fnp1

#  use trapezoidal formula to approximate 
#  integral over theta2 in [0, pi] of 
#                                      K(theta1, theta2) fn(theta2) sin(theta2)
   Ki <- K[, i]  # since K[i,j] is symmetric, this equals the vector 
#                  equal to the ith row of K 
#                  (K[i, ] whose entries corr. to values of theta2)
   trapez.sum <- sum(Ki * fn * sin.theta.vec) 
   value.of.integral.over.theta2 <- trapez.sum * delta.theta
   value.for.free.energy <-  B * value.of.integral.over.theta2 / 2
   vector.for.free.energy[i] <- value.for.free.energy  # for fn
   fnp1[i] <- exp(-Vvec[i] - value.for.free.energy)  # not yet normalized

   }   # end of loop to calculate each entry of fnp1

# Deal with the limits of finite precision in computations:
# if the argument v of the exponential function is too large
# (this is above 709 and below 710 on my computer), exp(v) will, in R, 
# return Inf (a reserved word in R).
# This would lead to error messages from downstream calculations.
# fnp1 is equal exp(-Vvec - vector.for.free.energy); note the minus signs
   maxv <- max(abs(vector.for.free.energy)) + max_mag_Vvec 
   if(maxv > 700) cat("an argument of exp will be ", maxv, "  at iteration ", iteration.num)
   if(maxv > 700) {
      stop("have situation of abs of argument of exp > 700; ending execution")
   }


#  now normalize fnp1
   theta.integral.of.fnp1 <- sum(fnp1 * sin.theta.vec) * delta.theta / 2  
   fnp1 <- fnp1 / theta.integral.of.fnp1  #  fnp1 is now normalized

#  finish calculation of integral operator K contribution to the discrete
#  free energy (for fn) 
   integral.operator.K.part.of.free.energy <- 
        sum(fn * vector.for.free.energy * sin.theta.vec) * delta.theta / 4  

#  calculate the contribution to the discrete free energy from 
#  the natural log term (for fn)

   log.part.for.F <- sum(fn * log(fn) * sin.theta.vec) * delta.theta / 2  

 V.part.for.F <- sum(fn * Vvec * sin.theta.vec) * delta.theta / 2  

   free.energy.for.fn <- 
               log.part.for.F + V.part.for.F + integral.operator.K.part.of.free.energy 
   free.energy.for.each.iteration  <-
               c(free.energy.for.each.iteration, free.energy.for.fn)

# don't flag if test for increase in free energy is close to 
# rounding error  
# flag only the first 12 non-decreases in F
   if (free.energy.for.fn >= free.energy.for.fnm1 + 1.0e-12) { 
      if(num_iter_F_did.NOT.decrease <= 11) {
        cat("iter = ", iter, 
            "  F did not decrease   only printing for first 12 instances", "\n")
      }
            num_iter_F_did.NOT.decrease <- num_iter_F_did.NOT.decrease + 1
   }

############################  test for convergence; if converged then 
#                             break   out of the iteration loop

############  estimate the approximation error as done in [HBW 1984] 
   
# do min.iterations iterations before start testing for convergence 

# protect if have complete convergence

maxdiff.for.iter <- max(abs(fn - fnm1))
L <- 0  # for exit at next line
if (maxdiff.for.iter < 1.0e-10) break

maxdiff.for.iter <- max(abs(fnp1 - fn)) 
L <- 0  # for exit at next line
if (maxdiff.for.iter < 1.0e-10) break

   L <- max(abs(fnp1 - fn)) / max(abs(fn - fnm1))
   err.est.denom <- abs(1 - L) + 1.0e-4  # protect against L = 1

   estimated.error <- (L / err.est.denom) * max(abs(fnp1 - fn)) / max(abs(fnp1))

   if ((estimated.error <= tau.conv) && (iter > min.iterations)) break

}  # end of iteration loop

#23456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890
###########  have finished iterations 

if (iteration.num == max.num.iter) print("WARNING iteration did not converge")

############################### calculate desired quantities 
############################### and return them in a list

f <- fnp1
f0 <- 3 * f[1] - 3 * f[2] + f[3]  # quadratic extrapolation

# do fpi by extrapolation on right end
len.fnp1 <- length(fnp1)

fpi <- fnp1[len.fnp1 - 2] - 3*fnp1[len.fnp1 - 1]  + 3*fnp1[len.fnp1]

fmid <- fnp1[N.theta.intervals / 2]

f.for.plotting <- c(f0, fnp1, fpi)  
# length is:   num.theta.points + 2 = (N.theta.intervals - 1) + 2

thetas.for.plotting <- c(0., theta.vec, pi)

len.F.m1 <- length(free.energy.for.each.iteration) - 1L
iterations.for.plotting <- 0L:len.F.m1

# P2(x)  the 2nd Legendre polynomial is  (3*x^2 - 1) / 2
# P4(x)  the 4th Legendre polynomial is  (35*x^4 - 30*x^2 + 3) / 8

pi.d4.minus.ave.W.gamma <- pi/4 - 
                          integral.operator.K.part.of.free.energy * 2 / B 
# note this is for the previous iterate, not the final iterate fnp1, 
# but should be quite close enough if the iteration converged

ave.ln.f <- log.part.for.F   # also for previous iterate

ave.P2.cos.theta <- 0.5 * sum((3*cos.theta.vec^2 - 1) * fnp1 * sin.theta.vec) *
                         delta.theta / 2

ave.P4.cos.theta <- 0.125 * sum((35*cos.theta.vec^4 - 30*cos.theta.vec^2  + 3) * 
                         fnp1 * sin.theta.vec) * delta.theta / 2

vector.of.values <- c(B, N.theta.intervals, N.phi.intervals, f0, fmid, fpi,  
                      L, pi.d4.minus.ave.W.gamma,
                      ave.ln.f, ave.P2.cos.theta, ave.P4.cos.theta, 
                      free.energy.for.fn,   
                      iteration.num, num_iter_F_did.NOT.decrease )

# use notation in [HBW 1984]
names.for.vector.of.values <- c("B", "num_theta_intervals", "num.phi.intervals",
         "f(0)", "fmid", "fpi", "L", "pi/4 - <<W(gamma)>>", "<ln(f(theta))>", "<P2(cos(theta))>", 
         "<P4(cos(theta))>", "free.energy",
         "num_iterations", "num_iter_F_did_NOT_decrease")

list.to.return <- list(fnp1 = fnp1, 
          f.for.plotting = f.for.plotting, thetas.for.plotting = thetas.for.plotting, theta.vec = theta.vec,
          free.energy.for.each.iteration = free.energy.for.each.iteration,
          iterations.for.plotting = iterations.for.plotting,
          vector.of.values = vector.of.values, 
          names.for.vector.of.values = names.for.vector.of.values)

return(list.to.return) 
      
}  # end of  calculate_ODF   function

return(calculate_ODF)  # calculate_ODF is defined when call compute_approx_ODF  

}  # end of  compute_approx_ODF
