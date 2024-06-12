# R-code-for-calculating-orientation-distribution-functions-for-rodlike-particles

README for the GitHub repository
R-code-for-calculating-orientation-distribution-functions-for-rodlike-particles
in the GitHub site:  https://github.com/AlanBerger/   

Last revision of this README: June 11, 2024

This repository contains R code for calculating orientation distribution functions
for rodlike particles using successive substitution and trapezoidal numerical integration as 
in J. Herzfeld, A. E. Berger, J. W. Wingate, A Highly Convergent Algorithm
for Computing the Orientation Distribution Functions of Rodlike Particles,
Macromolecules 1984, v17, 1718-1723,  
https://pubs.acs.org/doi/pdf/10.1021/ma00139a014
[HBW 1984]

The R code provided here can handle more general particle interaction terms W(gamma), 
and an optional applied field V(theta).

Note also

R. F. Kayser Jr. and H. J. Raveche', Bifurcation in Onsager's model of the 
isotropic-nematic transition, Physical Review A 1978, v17, 2067-2072 
[KR 1978]

## abbreviation used below: ODF is orientation distribution function. 


The R code provided in this repository: 

                  compute_approx_ODF_27Jan2024.R

is generally using the notation and 
setup in the Judith Herzfeld, Alan E. Berger and John W. Wingate paper,
except this code allows for an applied field V(theta) 
and for more general W(gamma), including W not satisfying W(pi - gamma) = W(gamma),
so calculations are done on theta in [0, pi] rather than assuming symmetry about pi/2 

Examples of using the functions in: compute_approx_ODF_27Jan2024.R  
are also contained in this GitHub repository, showing how figures in the manuscript:

Alan E. Berger, Conditions under which a natural iterative method for calculating the 
orientation distribution of rodlike particles decreases the free energy at each step, May 29, 2024
[Berger, 2024, submitted]

were generated.

Current computing power relative
to that in 1984 means that for these calculations one can favor simpler code 
over running speed. Rather than use Richardson extrapolation as done in
[HBW 1984] to obtain increased accuracy for quantities of interest, here the 
intent is to simply use larger numbers of grid points to obtain desired accuracy,
and can test for accuracy by successive doubling of the number of (equally spaced)
grid points. 

The R function:  compute_approx_ODF(N.theta.intervals, N.phi.intervals, W, 
                                                        V = function(theta) {0})
has the default for the applied field V being   V(theta) = 0 
and 
uses the function closure properties of the R programming 
language to do once, various setup calculations, and then have those results  
available to the function 
calculate_ODF(initial.f, B, tau.conv, max.num.iter, min.iterations)  
that compute_approx_ODF defines and returns. 

Here one is using the iterative method to find approximations to ODFs f(theta) 
that are local minima of the free energy function, denoted by F, and one is seaching 
and minimizing over non-negative functions f (defined at a set of discrete theta values) 
that are axisymmetric (so only depend on the spherical coordinate theta (not on phi)), 
and that satisfy the normalization: 

the integral over theta in [0, pi] of  (1/2) * f(theta) * sin(theta) = 1 
(as calculated using trapezoidal numerical integration).

#     Explanation of the arguments of these two functions.

N.theta.intervals is the number of equally spaced subintervals that [0, pi] is
divided into for doing trapezoidal numerical integration over theta in [0, pi].

Approximate values for ODFs are determined at the theta grid points in the interior 
of [0, pi].  sin(theta) factors in the integrals in the free energy F mean that values of an
ODF at theta = 0 and theta = pi have no effect on integrals over theta calculated using
trapezoidal numerical integration.

N.phi.intervals is the number of equally spaced subintervals that [0, 2pi] is
divided into for doing trapezoidal numerical integration over phi in [0, 2pi].

For a given value of the constant B in the free energy function F, and initial 
vector of values initial.f at the  grid points {Pj = j * pi / N.theta.intervals}, 
j = 1, 2, ..., N.theta.intervals - 1, the function calculate_ODF 
calculates an approximate ODF f  (a vector of values {fj} approximating the values
of an ODF at the grid points {Pj}). Under appropriate conditions on $W$,
f generally should approximate the values, at the grid points {Pj}, of an ODF that
is a local minimum of F. Convergence of the iteration is tested as in [HBW 1984], 
using as the convergence bound   tau.conv  

initial.f (a vector of values at the grid points {Pj}) is converted to be a (discrete) ODF 
that is positive (any initial.f values <= 0 are redefined to be slightly above 0, 
to avoid issues with the log function), and then initial.f is scaled (normalized)  
so that the integral over theta in [0, pi] of  (1/2) * initial.f(theta) * sin(theta) = 1
as calculated using trapezoidal numerical integration on the grid points {Pj}
(details are in the file compute_approx_ODF_27Jan2024.R in this repository).

max.num.iter is the maximum number of iterations allowed, and min.iterations
is the number of iterations to be done before testing for convergence.

W(gamma) is the function in the particle interaction term integral in the free energy 
function F (see [HBW 1984] for definitions and notation). The code here allows for
more general functions W and an applied field V(theta).
Under appropriate conditions* on W, the iteration
as applied to continuous ODFs strictly decreases the free energy at each step 
(unless one was at a fixed point of the iteration) [Berger, 2024, submitted].
This result has strong implications for 
convergence of the sequence or a sub-sequence of the iterations 
to a solution of the equation, given by the calculus of variations, 
that a local minimum of the free energy will satisfy. 

The condition * is:  If one expands W(gamma) as a sum of Legendre polynomials, i.e.,  
sum from m=0 to infinity of w_m P_m(cos gamma); the conditions on W are that: 
w_m <= 0 for m > 0 and sum |w_m| is finite.

The Rmd (R markdown) file  R-code-for-calculating-ODFs-and-example-runs.Rmd 
in this repository uses the R function
compute_approx_ODF (which returns the function calculate_ODF), 
and contains R code that will generate output and plots for example runs.

Its output (the corresponding PDF file) is also in this repository.

## Note ([Berger, 2024, submitted]):

If W(pi - gamma) = - W(gamma), and V = 0, and
the initial value f0 for the ODF to start the iteration
is symmetric about theta = pi / 2 (i.e., f0(pi - theta) = f0(theta)),
then all the iterates f1, f2, ... will equal the isotropic solution since the 
integral from the particle interaction term will be 0. Note, adding a constant to 
the free energy $F$  will not change its local minima, or change the ODFs from
the iteration, so if V = 0, and  Y = (W - a fixed constant) satisfies,
Y(pi - gamma) = - Y(gamma), then the note in 
this paragraph applies.

Note also that if V(pi - theta) = V(theta) and f(theta) is a local minimum of F, 
then g(theta) defined equal to f(pi - theta) will also be a local minimum. In the same vein, 
if V(pi - theta) = V(theta) and the initial (starting) ODF f0 satisfies 
f0(pi - theta) = f0(theta), then all the iterates will satisfy this symmetry.

In general, one should use several 
different initial ODFs to start the iteration, (including initial ODFs that are 
not symmetric about theta = pi / 2) to have a better chance of finding
all relevant local minima of the free energy. 

Rmd files that generate tables of values of quantities of interest for 
many values of B (1 row for each value of B) for the hard core
reference system (W(gamma) = sin(gamma) and V(theta) = 0) 
starting with an axial shaped initial ODF (maxima at theta = 0 and pi) 
and (2nd Rmd file) starting with a planar 
shaped initial ODF (maximum at theta = pi/2) are also in this repository, 
along with the results (that are Excel files derived from the tab delimited text 
files output by these Rmd files). These files are dated February 8, 2024

To download an R file or an Rmd file (as a text file) to get the R code,
click on the file name in this repository, then in the display one sees, 
(with a Windows computer) right click on "raw" and select "save link as", 
or in general, click on the downarrow icon (to the right of "raw") to download
the file.

If one wants to copy any of the R code, do that 
from the R or Rmd file (which are plain text files).

When any R markdown file in this repository is run through R's knitr, 
it produces the corresponding
PDF file that contains output from the R code including plots if present 
(and tab delimited output text files if coded for)
(the pdf file will be included in this repository for the Rmd files that generated
plots for [Berger, 2024, submitted]). 


For the hard core reference system, the function W in the free energy function
equals sin(gamma) (and V = 0). This $W$ satisfies the conditions 
w_m <= 0 for m > 0 and sum |w_m| is finite ([KR 1978]).

gamma is the angle between two spherical coordinates (theta1, phi1) and (theta2, phi2),
which by the law of cosines equals
    acos(sin(theta1) * sin(theta2) * cos(phi1 - phi2) + cos(theta1) * cos(theta2))
(see for example [KR 1978], [Berger, 2024, submitted]).
(Rather than using Greek symbols, I am using variable names that are in the R code.)

Trapezoidal numerical integration is used to approximate the value of 
integrals, and the presence of sin(theta) factors in the free energy function 
results in the values of an ODF f at the points corresponding to theta = 0 and 
theta = pi having no effect on the value of the discretized free energy. 

An approximation for the value at theta = 0 of a discrete ODF f  
(defined at the {Pj} points) is obtained as the value at theta = 0 of   
the quadratic function that has the values f1, f2, f3 at the 
first three theta grid points on the left inside (0, pi), (as done in [HBW 1984]),
and similarly for the value of f at theta = pi.
While values of an approximate (discretized) ODF at 
theta = 0 and theta = pi do not matter for calculating integrals involving
a sin(theta) factor via the trapezoidal formula, these values are used in 
doing plots.
