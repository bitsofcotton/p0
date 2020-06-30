# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.
This suppose original data stream is in L2(R).
And there exists the sampling theorem.

# How to use:
    ...
    P0<double> p(range);
    SimpleVector<double> b(range);
    // ... b operation ...
    pred = p.next(b);
    ...

# How to use (commandline):
    ./p0 <range> <loop> < data.txt

