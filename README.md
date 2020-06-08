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
If each of the error of them is ec_k, es_k, set the geometric mean of ec_k and es_k to r/(1\pm r), then,
the maximum error correction ratio is : x^2+y^2=1^2 case, x's error ratio into r : 1 \pm r, then, atan(r/sqrt(1\pm r^2))/atan(sqrt(1\pm r^2)/r) ratio, then, thanks to maxima, this converges around &pi;/2\*r+4/&pi;^2\*r^2+... if |r| < 1.  
So this concludes if p0 is stable prediction for the series in this meaning, P0C corrects some of the prediction.
