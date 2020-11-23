# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.
This suppose original data stream is in L2(R),
and original data stream isn't a essential singularity function in sampling meaning.
(There's a typical case that sampling frequency isn't enough, nor,
 original function is increasing or decreasing speed like sin(1/x) case. )

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
