# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.
This suppose original data stream is in L2(R).
And there exists the sampling theorem.

# How to use:
    ...
    P0B<double> p(range);
    xnext = p.next(x);
    ...

# How to use (commandline):
    ./p0 <range> < data.txt

# Tips
This better predicts with 2 step pairs on the complement meaning. This stands:
xnext = (t ++) & 1 ? p.next(x) : q.next(x)
case, delta from the next step's p's or q's complement on x is better predicted.
So the sampling theorem is the matter.
