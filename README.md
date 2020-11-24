# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.  
This suppose original data stream is in C^&infin;(R) in before sampling meaning,
and in sampled meaning, suppose original data stream is in L2(R).

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

# pf
If original function is in C^&infin;, there f(z+\bar{z}) in C^&infin; on z in C.
So rotate &pi;/4, f is holomorphic because the function is described as f(z), so d/d\bar{z}f = 0.
So with gutzmer's inequation, f have laurent series and as the upper bound of coefficients,
we can cut them as error in numerical calculation on some cut off.  
And in sampled meaning, we can't suppose real cut off on coefficients.  
And, the cut off error can be configured with P0 class variable.
