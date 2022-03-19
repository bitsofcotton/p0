# bitsofcotton/p0
Generic 1-variable predictor when input values are ordered same as the order in R.

There's a plenty of the room to make this into n-variable predictor, but this repository won't implement such of them.

# Contexts
There exists discrete fourier transform on given (same interval) series (This is well described on everywhere.).
And, if we make DFT and IDFT on them, there exists differential on them in DFT meaning.
(In continuous meaning, this is described a little on the books I refered,
 but, in discrete meaning, I can't find preceding results but might be exists.)

# How to use:
    northPole<double, P0<double, /* feeder */> > p(P0<double, /* feeder */>(range, step));
    // pnext<T>(range, step) for prediction vector.
    ...
      xnext = p.next(x);

# How to use (commandline):
    ./p0 <avg> < stream-next-sign-unknown.txt
    # avg < 0 for average origin on pred range, otherwise, 0 origin.
    # abs(avg) specifies middle frequency band size.

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
other than the condition below C0. And if then, tan(&theta;) should not have pole or infinite value.
exp(Sum log(z-|z|cis(&pi;/4+t_k)))\*(Sum((f(z+\bar{z})+i\*f(z-\bar{z})\*tan(&pi;/4+t\_k))/(z-|z|cis(&pi;/4+t\_k)) =: f(z+\bar{z})\*g(z)+f(z-\bar{z})\*h(z) = f(z+\bar{z})\*(g\_real(z+\bar{z})+i\*g\_imag(z-\bar{z}))+f(z-\bar{z})\*(h\_real(z+\bar{z})+h\_imag(z+\bar{z})) == (f(z+\bar{z})\*g\_real(z+\bar{z}) - f(z-\bar{z})\*h\_imag) + imaginary part, so replace on (1+i)t, F(x):=(f(x)\*g\_real(x)-f(x)\*h\_imag(x))==f(x)\*(g\_real(x)-h\_imag(x)), F is holomorphic on some small range (f in C0) and Re F = f(x)\*some G(x), G(x) not depends f, so f(x) / G(x) =: f2(x), same method, f(x) will be holomorphic on some small range.

And, we can weaken this condition with cauchy's integrate theorem on ja.wikipedia.org link (doi:10.1090/S0002-9947-1900-1500519-7), C1 condition to C0 condition (uses epsilon delta on multiple of different differential value with limit needs smaller than continuous condition).  

N.B. if there's no C0 condition, with the condition x_next:=integrate^x_now(f(x) - some &alpha;), the prediction is valid because each of them are structure of the sum between first point and each of them. So if we can define integrate(original f) (and if it's continuous enough), the prediction is valid if original stream probability can be riemann integrable (and it's continuous enough).

N.B. DFT differential itself is ideal for trigometric function multiply sums. So applying tilt itself returns curved result, but this is reasonable one on the range with IDFT * DFT condition.

N.B. p0 uses northPole class to hack essential point whether or not if near the defined region.

# Tips
This p0 uses weak differential, so integrate(diff) != id because of complex part is omitted (but Re diff(diff) is no matter if we calculate in R or C).

# Tips
Frequency space prediction is equivalent to differential/integratial prediction in this. But it is equivalent to plain prediction in this.

# Tips
P0 calculates left differential by periodical one. To vanish right hand differential, we take sum of each range on next step vector. Whether this works well or not depends on left hand side distribution. But with some of differential sequence meaning, if we have some differential sequence average as 0 value, we can use them bared prediction.

# Tips
northPole class is formal laurent series essential point hack, if not both side is essential, multiplication inverse series works well, otherwise, the essential point hack enforces them into north pole near the defined region.
We need to northPole twice because only once sometimes doesn't converges bothside coefficients to 0 in general. But with the range they converges (non infinite values on any of samplepoints we have on the range), we only need to northPole once.
Because of this, we can handle any function which we can riemann integrable. Since we take probability, so if lebesgue integrate on original function is in C0, we can handle any function.

# Tips
There's a plenty of a space to extend this algorithm with &uparrow;, &downarrow; definition (they concludes exp(exp(...)) in C, inverse function is defined in randtools meaning inverse limit, but it's also in H if second operand is in C\\R.), they causes 3-variable function causes all of (x, f(x), status) handled on the series. With them, we can handle n-variable on reasonable function (because (n-1)-status block causes only the (x, f(x), status) triplet with extreme accuracy.).
But with this 1-variable predictor : (x, f(x)) any correspondence, if the parameter is large enough to input stream (as all of the status is on the taylor size), it's equivalent to the one.

-&gt; Knuth tower symbol seems to make some symmetric result on some of imaginary axises. This causes only the symmetric series on some of taylor series on H (with original series as j axis plotted). So this causes as the almost same prediction result via quarternion discrete fourier transform on symtax meaning, and they causes only the frequency-phase phase shifted results on the transformation. So we can define the differential on them on 2 of axises, but they seems to non information-rich one. And on information-rich meaning, we can apply this p0 to frequency series on complex space, they shall be equivalent to the one.

# Tips
With referring wikipedia article: https://ja.wikipedia.org/wiki/%E5%AE%9F%E9%96%89%E4%BD%93 (2022/03/19), any weakly o-minimal structure is real closed field causes dim A = 1, 2, 4, 8 (https://ja.wikipedia.org/wiki/%E5%A4%9A%E5%85%83%E4%BD%93 (2022/03/19)), so if we can define {f(x) \| x in R} as semi-ordered one, octanion is enough dimension to calculate.

# Another Download Sites (Closed)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/

