# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.

# How to use:
    P0<double> p;
    SimpleVector<double> b(range);
    ...
    pred = p.next(b);
    // Or we can use:
    P0B<double> p(range);
    ...
    xnext = p.next(x);

# How to use (commandline):
    ./p0 <range> < data.txt
    # N.B. range > 0 for random walk (C0), range < 0 for random value itself.

# pf
If original function is in C1, there exists f(z+\bar{z}) in C1 on z in C.
So rotate &pi;/4, f is holomorphic because the function is described as f(z), so d/d\bar{z}f = 0.
So with gutzmer's inequation, f have laurent series and as the upper bound of coefficients,
we can cut them as error in numerical calculation on some cut off.  
And in sampled meaning, we can't suppose real cut off on coefficients.
So the cut off error can be configured with P0 class variable.

And, we can weaken this with cauchy's integrate theorem on ja.wikipedia.org link (doi:10.1090/S0002-9947-1900-1500519-7), C1 condition to C0 condition (epsilon delta).

And, if there's laurent series isn't essential singularity point on the domain, with f(x_0+\[0,1\[) case, the series can be described as taylor series around there (1+x+...+x^n+..., x in \[0,1\[). So we can suppose if, the case below and, around x=x_0, D:={|x-x_0|<1}, there exists taylor series descriptions near them if there's a real function.
