# R-code-for-calculating-orientation-distribution-functions-for-rodlike-particles

R code for calculating orientation distribution functions for rodlike particles 
using successive substitution and trapezoidal numerical integration as 
in Herzfeld, Berger, Wingate, Macromolecules 1984, v17, 1718-1723, 
https://pubs.acs.org/doi/pdf/10.1021/ma00139a014

README for the GitHub repository
R-code-for-calculating-orientation-distribution-functions-for-rodlike-particles

Alan E. Berger December 30, 2023    revised Jan 28  2024

The R code provided here uses the iterative method indicated by 
the calculus of variations, discretized 
as in J. Herzfeld, A. E. Berger, J. W. Wingate, A Highly Convergent Algorithm
for Computing the Orientation Distribution Functions of Rodlike Particles,
Macromolecules 1984, v17, 1718-1723, [HBW 1984], 
to calculate approximations to Orientation Distribution Functions 

Note also

R. F. Kayser Jr. and H. J. Raveche', Bifurcation in Onsager's model of the 
isotropic-nematic transition, Physical Review A 1978, v17, 2067-2072 [KR 1978]


abbreviation used below: ODF is orientation distribution function


The R code provided in this repository 

                  compute_approx_ODF_27Jan2024.R

is generally using the notation and 
setup in the Judith Herzfeld, Alan E. Berger and John W. Wingate paper,
except this code allows for an applied field V(theta) 
and for more general W(gamma), including W NOT satisfying W(pi - gamma) = W(gamma),
so calculations are done on [0, pi] rather than assuming symmetry about pi/2 

Current computing power relative
to that in 1984 means that for these calculations one can favor simpler code 
over running speed. Rather than use Richardson extrapolation as done in
[HBW 1984] to obtain increased accuracy for quantities of interest, here the 
intent is to simply use larger numbers of grid points to obtain desired accuracy,
and can test for accuracy by successive doubling of the number of (equally spaced)
grid points. 

The R function  compute_approx_ODF(N.theta.intervals, N.phi.intervals, W, 
                                                        V = function(theta) {0})
has the default for V being   V(theta) = 0 
and 
uses the function closure properties of the R programming 
language to do once, various setup calculations, and then have those results  
available to the function 
calculate_ODF(initial.f, B, tau.conv, max.num.iter, min.iterations)  
that compute_approx_ODF defines and returns. 
For a given value of the constant B in the free energy function F, and initial 
vector of values initial.f at the  grid points {Pj = j * pi / N.theta.intervals}, 
j = 1, 2, ..., N.theta.intervals - 1, the function calculate_ODF 
calculates an approximate ODF f (a vector of values {fj} approximating the values
of an ODF at the grid points {Pj}). 
f generally should approximate the values, at the grid points {Pj}, of an ODF that
is a local minimum of F. Convergence of the iteration is tested as in [HBW 1984], 
using as the convergence bound   tau.conv  

max.num.iter is the maximum number of iterations allowed, and min.iterations
is the number of iterations to be done before testing for convergence.

W is the function in the particle interaction term integral in the free energy 
function F (see [HBW 1984] for definitions and notation). The code here allows for
more general functions W, and also can show examples where the iteration does not
always decrease the free energy. Under appropriate conditions on W, the iteration
as applied to continuous ODFs strictly decreases the free energy at each step 
(unless one was at a fixed point of the iteration) Alan E. Berger, 
manuscript in preparation. This result has strong implications for 
convergence of the sequence or a sub-sequence of the iterations 
to a solution of the equation, given by the calculus of variations, 
that a local minimum of the free energy will satisfy.

The Rmd (R markdown) file  R-code-for-calculating-ODFs-and-example-runs.Rmd 
in this repository uses the R code for
the function  compute_approx_ODF (which returns the function calculate_ODF) 
and contains R code that will generate output and plots for example runs.

At present (Jan 28, 2024)
R-code-for-calculating-ODFs-and-example-runs.Rmd 
is self contained, but uses an older version of  
of  compute_approx_ODF  and   calculate_ODF

To download the Rmd file (as a text file) to get the R code,
click on the Rmd file name in this repository, then in the display one sees, 
(with a Windows computer) right click on "raw" and select "save link as", 
or in general, click on the downarrow icon (to the right of "raw") to download
the Rmd file.

If one wants to copy any of the R code, do that 
from the R or Rmd file (which are plain text files).

When this R markdown file is run through R's knitr, it produces the corresponding
pdf file that contains output from the R code including plots
(the pdf file will be included in this repository).

Note also: 

A. E. Berger, Analysis of a constrained minimization problem modeling the orientation 
distribution of rod-like particles, Nonlinear Analysis, Theory, Methods & Applications 1987, 
v11, 719-731 [B 1987].

For the hard core reference system, the function W in the free energy function
(W is a function whose argument is denoted by gamma) equals sin(gamma)
gamma will be equal
    acos(sin(theta1)*sin(theta2)*cos(phi) + cos(theta1)*cos(theta2))
(see for example [KR 1978], [B 1987])

Trapezoidal numerical integration is used to approximate the value of 
integrals, and the presence of sin(theta) factors in the free energy function 
results in the values of an ODF f at the points corresponding to theta = 0 and 
theta = pi having no effect on the value of the discretized free energy. 
When W(gamma) equals sin(gamma), by Theorem 1.2 or Lemma 3.1 
in [B 1987] any ODF that is a  
local minimum of the actual (not discretized) free energy function F will 
have (at least) 2 continuous derivatives 
on [0, pi] and its first derivative will be 0 at theta = 0, pi/2 and pi

An approximation for the value of the ODF f at theta = 0 
is obtained as the value at theta = 0 of   
the quadratic function that has the values f1, f2, f3 at the 
first three grid points on the left inside (0, pi), (as done in [HBW 1984]),
and similarly for the value of f at theta = pi.
While values of an approximate (discretized) ODF at 
theta = 0 and theta = pi do not matter for calculating integrals involving
a sin(theta) factor via the trapezoidal formula, these values are used in 
doing plots.
