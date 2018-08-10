## Matrix related computational methods implemented functionally
This repository contains source files for a program designed to perform the following computations:
* Solving SLE using Gaussian elemination
* Solving SLE using Gaussian elemination with leading element
* Solving SLE using successive over-relaxation
* Computing matrix determinant using Gaussian elemination
* Computing inverse matrix using Gauss-Jordan elimination
* Computing the condition number of a linear map
<br>
  The project was initially an academic task, but other than that I consider it a decent example of my aproach to dealing with Haskell's clumsiness when it comes to interacting with system calls. I solve this by interfacing haskell with a small snipet of C code. In fact any imperative language will do, it's just that C seems a good choice, since the interface module I've developed to implement such an interface is universal, simple in exploitation and can be reused in future projects or even made into a built-in library. That makes it an easy to set up and an inspiringly beneficial fusion.
