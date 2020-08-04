# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.
This suppose original data stream is in L2(R).
And there exists the sampling theorem.

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

