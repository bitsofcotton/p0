# bitsofcotton/p0 (closed)
Generic predictor for Riemann measurable input streams and their variants (on aleph_1 cf. C).

There's a plenty of the room to make this into n-variable predictor, but this repository won't implement such of them but pseudo discrete one is.

Also there's another plenty of the room to make th. Peter-Weyl based transformations as some of the normal orthogonal base but we cannot define differential of them in generic ones.

Also, there's some of the L-function based transformation s.t. \[\[Sum chi(k)(0/n)^k\*1, ...\],\[..., Sum chi(k)(1/n)^k\*2, ...\],...\] based matrix but this is obscure. With lacking foundation but they could be transformation on Archimedean helicoid also a single loop on P0DFT about prediction.

# Contexts
There exists discrete fourier transform on given (same interval) series (This is well described on everywhere.).
And, if we make DFT and IDFT on them, there exists differential on them in DFT meaning.
(In continuous meaning, this is described a little on the books I refered,
 but, in discrete meaning, I can't find preceding results but might be exists.)

# How to use:
    SimpleVector<double> buf;
    ...
    xnext = p0maxNext<double>(buf);

# How to use (commandline):
    ./p0(-(32|64)) <length>? <step>? < stream.txt
    # 0 == length for using whole length shallow.
    # 0 <  length for using recent length input numbers.
    # step < 0    for chain condition.

# Tips
Some of the important tips also implanted into lieon.hh as a comment.

# Another Download Sites (Closed)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/

# Leave
2023/02/28
2023/03/09 bug fix after close #1.
2023/03/13 integrate all files into lieonn.hh after close #2.
2023/03/18 offset pnext, eliminate recursive pnext. after close #3.
2023/03/24 code clean, after close #4.
2023/03/29 pnext needs 3 &lt;= r, after close #5.
2023/03/31 merge prand.
2023/04/02 merge catg fix.
2023/04/03 merge.
2023/04/04 update readme.
2023/04/21 only single side northPole method, we don't need them on stability in such case.
2023/06/24 fix to avoid observation matters.
2023/07/01 argv range minimum value change.
2023/07/07 update .cc comment.
2023/10/30 add p0l to add some continuity. update readme.
2024/05/06 integrate p0l into p0, p0's argv\[1\]\<0 is integrated into bitsofcotton/p1.
2024/05/31 compile jammer.
2024/06/01 fix JAM.
2024/06/02 fix _JAM replace with _JAM_, this caused pseudo correct results.
2024/06/02 JAM into p2/cr.py.
2024/06/13 add p0i.cc.
2024/06/15 punch p0i.cc, it's no use.
2024/06/15 add progression.
2024/06/16 add progression \<0 arg.
2024/06/17 fix progression.
2024/06/19 merge latest lieonn.
2024/06/21 merge latest lieonn. INCLUDES command line argument change.
2024/06/22 update readme, merge latest lieonn.
2024/06/23 large change around class instance initializer, also have progression short range fix.
2024/06/23 fatal fix around last update on lieonn.hh, readme.md.
2024/06/24 fix addp == true progression case.
2024/06/26 fix Ppersistent.
2024/07/07 Pprogression uses shorter range but enough internal states.
2024/08/27 add plenty of the room another directions.
2024/08/29 add low quality notice around plenty of the room.
2024/09/21 add p04.cc . elim p04.cc, it's from mistake.
2024/09/24 revert p0.cc as to use only raw P0maxRank. Pprogression is duplicate and exhaust of the resources.
2024/09/28 add step option.
2024/11/16 add p0 0 ... command.
2024/11/19 speed/accuracy fix on taylor on lieonn.hh. this change leads us doubles pnext r variable.
2024/12/02 taylor function improvement, reclose with this.
2024/12/05 exchanged argv[1] and argv[2] meanings.
2024/12/08 new p0 0 0 command line option is now default.
2024/12/09 fix readme calculation order as to match implementation.
2024/12/13 should really leave, close readme.md.
2024/12/14 shoud last update readme.md.
2025/02/07 add _CHAIN_ compile option.
2025/02/07 absent and integrate lastup into p2/persistent.cc.
2025/02/13 compat with latest lieonn.hh.
2025/02/27 elim step parameter.
2025/03/08 improves argv[1] == 0 memory usage instead of caching pnext results, this slow downs another argv[1] > 0 cases.
2025/03/09 merge latest lieonn.
2025/03/17 revert step param.
2025/04/17 merge latest lieonn, brush up length == 0 case move into lieonn.hh.
2025/04/18 eliminate step param, they doesn't improves result.
2025/04/19 merge latest lieonn.
2025/04/20 compile out exp-log scale to -D_NONLIN_ option. argc change.
2025/05/17 add p0p.cc for persistent.
2025/05/19 fix p0p function to have correct meaning.
2025/05/19 fix and now works p0p function with {-1, 1} input stream.
2025/05/19 p0p.cc with -D_CHAIN_ compiler option to have p0, p1, p2 series compat output.
2025/05/23 elim p0p.cc, it's integrated to 'p0 d' argv1, either refactoring.
2025/05/25 merge latest lieonn, de-select pbond class, argv meaning change.
2025/06/06 catch excluded p0 \[12\] arg into p0.cc.
2025/06/07 merge latest lieonn.
2025/06/08 argv meaning change because we don't need pnext, diff real vector output. freeze p0 command.
2025/06/11 compat compile option to gcc4.2.1.
2025/06/12 compat compile option with one variant of gcc2.95.3.
2025/06/17 merge latest ddpmopt fix.
2025/06/19 merge latest ddpmopt fix.
2025/06/25 readme.md moved into implant lieonn.hh as a comment.
2025/06/28 refactor and fix around lieonn, re-compat with gcc2953.
2025/06/29-30 merge latest ddpmopt result.
2025/07/01 merge latest ddpmopt result.
2025/07/02-03 merge latest ddpmopt resut, no logic change.
2025/07/04 merge latest ddpmopt resut, some speed remedy.
2025/07/06 merge latest lieonn.
2025/07/13 merge latest lieonn.
2025/07/14-16 merge latest lieonn.
2025/07/17-19 merge latest lieonn.
2025/07/20 merge latest lieonn.
2025/07/24 merge latest lieonn.
2025/07/25 merge latest lieonn.
2025/07/26-28 merge latest lieonn.
2025/08/01 merge latest lieonn, pseudo multiple variable pred not exact.
2025/08/01 re-enable step option.
2025/08/03 fix last up shift delay, merge latest lieonn, however something buggy.
2025/08/04-06 merge latest lieonn.
2025/08/11 merge latest lieonn.
2025/08/12-16 merge latest lieonn.
2025/08/17-23 merge latest lieonn.
2025/08/25 merge latest lieonn.
2025/09/01 merge latest lieonn.
2025/09/05 merge latest lieonn.
2025/09/25 merge latest lieonn.
2025/10/06 add step < 0 option for first half prediction with applying second half and return whole.
2025/10/15 fix step < 0 case last column. also eliminate step < 0 case with by replacing the predictor as a linear one.
2025/10/16 we need step < 0 actually, so revert and fixed last bug.
2025/10/30 revert eliminte step < 0 condition, we don't need them.

