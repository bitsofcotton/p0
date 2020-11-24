# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.
This suppose original data stream is in L2(C^&infty;(R)).

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

# Configurable variables:
P0's integrate number is configurable. It is ideally around ceil(lg(dim pole)).

# pf
If original function is in C^&infty;, there exists f(z+\bar{z}) in C^&infty; on z in C.
So rotate &pi;/4, f is holomorphic because the function is described as f(z), so d/d\bar{z}f = 0.
So with gutzmer's inequation, f have laurent series and as the upper bound of coefficients,
we can cut them as error in numerical calculation.
