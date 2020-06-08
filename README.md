# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.
This suppose original data stream is in L2(R).
And there exists the sampling theorem.

# How to use:
    ...
    P0<double> p(range);
    SimpleVector<double> b(range);
    // ... b operation ...
    pred = p.nextP(b);
    ...

# How to use (commandline):
    ./p0 <range> < data.txt

# Tips
P0C implements circular error correction utility, this stands:
for any k, x_k in \[-1/2,1/2\], we set y_k from x_k into circular geometry, then, we predict with cos(y_k), sin(y_k).
If each of the error of them is ec_k, es_k, set the geometric mean of ec_k and es_k to e_k, then, thanks to maxima,
the maximum error correction ratio is : x^2+y^2=1^2 case, x's error ratio into r : 1 - r, then, atan(r/sqrt(1-r^2))/atan(sqrt(1-r^2)/r) ratio this converges around &pi;/2r-1 if |r| << 1.
