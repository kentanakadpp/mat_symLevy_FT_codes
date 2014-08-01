Matlab codes for the symmetric Levy processes
====================

These are Matlab R2013a (8.1.0.604) codes for numerically solving the Kolmogorov forward equations for the symmetric L´evy processes. These codes were used to yield the results of the numerical experiments in my article

"A fast and accurate numerical method for the symmetric L´evy processes 
 based on the Fourier transform and sinc-Gauss sampling formula"
(http://arxiv.org/abs/ ).

See this article for the technical details of the numerical method. In particular, Steps 1, 2, and 3 below are explained in the article.  

The roles of the source files and their relationships are explained below. 

[Functions]
- LevyPDE_SolFunc_IIone.m --- Function to compute the numerical solutions using the following functions as subroutines. In this function, SG_IndefInt_sym.m is executed once. 
  - DE_NFFT.m --- Function to compute the Fourier transform in Step 1 through the DE formula for the Fourier transforms and the nonuniform FFT. The DE transform is implemented in the following file.
    - phi.m   
  - SG_IndefInt_sym.m --- Function to compute the indefinite integral in Step 2 through the sinc-Gauss indefinite integration formula. This function uses the fractional FFT, which is implemented in the following file.
    - FFFT.m 
  - ContEulerFFFT.m --- Function to compute the inverse Fourier transform in Step 3 through the formula for the Fourier transform with continuous Euler transform and the fractional FFT. Thus, this function also uses the following file.
    - FFFT.m 
- LevyPDE_SolFunc_IItwo.m --- Function to compute the numerical solutions using the same functions as LevyPDE_SolFunc_IIone.m. In this function, SG_IndefInt_sym.m is executed twice. 



