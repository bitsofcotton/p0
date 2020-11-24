# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.  
This suppose original data stream is in C1(R) in before sampling meaning,
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
If original function is in C1, there exists f(z+\bar{z}) in C1 on z in C.
So rotate &pi;/4, f is holomorphic because the function is described as f(z), so d/d\bar{z}f = 0.
So with gutzmer's inequation, f have laurent series and as the upper bound of coefficients,
we can cut them as error in numerical calculation on some cut off.  
And in sampled meaning, we can't suppose real cut off on coefficients.
So the cut off error can be configured with P0 class variable.

# Tips
f in C1 supposes : exists! f': continuous.  
This is constructive with using discrete {a\_k}\_k, f(x)=Sum\_{n=0}^{n=x}Sum\_{k=0}^{k=n}a\_k,
so this is the predictor on Sum\_{k=0}^{k=x}a\_k sequence of points.
