# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.
This suppose original data stream is in L2(R).
And there exists sampling theorem.

# How to use:
    ...
    p0_t p(range);
    pred = p.next(d);
    ...

# How to use (commandline):
    ./p0 <range> < data.txt