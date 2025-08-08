This project uses code developed by Martin Cox and David Borchers for density estimation in the presence of a density gradient.

You can install ```nupoint``` using:

```{r installnupoint,eval=FALSE}
devtools::install_github('martinjamescox/nupoint')
```

There is a vignette available describing analysis of a non-uniform density gradient and there are also two pdfs with worked examples of the uniform and non-uniform density gradients in nupoint/inst/extdata/

Publications related to this package are:

M.J. Cox, D.L. Borchers, D.A. Demer, G.R. Cutter, and A.S. Brierley. 2011. Estimating the density of Antarctic krill (Euphausia superba) from multi-beam echo-sounder observations using distance sampling methods. Journal of the Royal Statistical Society: Series C (Applied Statistics), 60(2):301-316.

M.J. Cox, D.L. Borchers and N. Kelly. 2013. nupoint: An R package for density estimation from point transects in the presence of non-uniform animal density Methods in Ecology and Evolution 4(6):589-594. DOI: 10.1111/2041-210X.12058

The project also uses a subset of data from:

Arranz, P., Borchers, D.L., de Soto, N.A., Johnson, M.P. and Cox, M.J., 2014. A new method to study inshore whale cue distribution from land-based observations. Marine Mammal Science, 30(2), pp.810-818.