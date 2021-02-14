# bitsofcotton/p0
Generic predictor that not depends on data itself but doesn't win good randoms.  
This suppose original function in the domain needs in L2(R) to guarantee no essential singularity points on there.

# Contexts
There exists discrete fourier transform on given (same interval) series (This is well described on everywhere.).
And, if we make DFT and IDFT on them, there exists differential on them in DFT meaning.
(In continuous meaning, this is described a little on the books I refered,
 but, in discrete meaning, I can't find preceding results but might be exists.)

# How to use:
    P0<double> p;
    SimpleVector<double> b(range);
    ...
    pred = p.next(b.size());
    // Or we can use:
    P0B<double> p(range);
    ...
    xnext = p.next(x);

# How to use (commandline):
    ./p0 <range> < data.txt
    # range < 0 for arctan(input) prediction,
    # 0 < range for plain prediction.
    # the range actually used is range! .

# Proof
If original function is in C1, there exists F(z,&theta;) := f(z+\bar{z})+i\*f(z-\bar{z})\*tan(&theta;) in C1 on z in C.
So each &theta;, exists F: holomorphic function at some axis that real axis is same as f.
So &theta; = &pi;/4 \pm &epsilon;, if we have a area, with gutzmer's inequation and edit integrate path,
f have laurent series and as the upper bound of coefficients,
we can cut them as error in numerical calculation on some cut off.
If f has a laurent series around {x| |x-a|<1} and no essential singular point on there,
taylor series description exists in weak differential meaning because 1+x+x^2+...
converges and x^k -&gt; 0 and k -&gt; \infty on there.  
There's possibility to have area with lagrange multipliers with limit operations, but if then, it has to use the condition of f
other than the condition below C0.

And, we can weaken this condition with cauchy's integrate theorem on ja.wikipedia.org link (doi:10.1090/S0002-9947-1900-1500519-7), C1 condition to C0 condition (uses epsilon delta on multiple of different differential value with limit needs smaller than continuous condition).  

And, if there's no C0 condition, with the condition x_next:=integrate^x_now(f(x) - some &alpha;), the prediction is valid because each of them are structure of the sum between first point and each of them. So if we can define integrate(original f) (and if it's continuous enough), the prediction is valid a.e..  

N.B. DFT differential itself is ideal for trigometric function multiply sums. So applying tilt itself returns curved result, but this is reasonable one on the range with IDFT * DFT condition.

# Tips
If we don't predict well, please refer p2.

# Tips
This p0 uses weak differential, so integrate(diff) != id because of complex part is omitted (but Re diff(diff) is no matter if we calculate in R or C).

# Another Download Sites
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/

# Archive
actually freeze with this.

