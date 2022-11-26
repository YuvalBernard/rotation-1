<TeXmacs|2.1>

<style|<tuple|generic|SIUnits-simple>>

<\body>
  <doc-data|<doc-title|Numerical Solution of Bloch Equations for Liquid
  Phase>>

  <underline|<strong|Goal:> >

  Write a code which solves the Bloch equations for a vector of the form
  \ <rigid|<math|<wide|M|\<vect\>>=[M<rsub|x>,M<rsub|y>,M<rsub|z>]>>. It
  should act on an initial magnetization vector
  <math|<wide|M|\<vect\>>=M<rsub|0>[0,0,1]> and use a matrix operator which
  includes the static external field <math|H<rsub|0>>, the perturbing field
  <math|H<rsub|1>> rotating at frequency <math|\<omega\><rsub|z>>, the
  relaxation constants <math|T<rsub|1>,T<rsub|2>> and the offset from
  resonance.

  <section|Different representations of Bloch equations>

  Assuming <math|H<rsub|0>> is statically aligned on the <math|z> axis, that
  <math|H<rsub|1>> is initially aligned on the <math|x> axis and alternates
  at frequency <math|\<omega\><rsub|z>=-\<gamma\>H<rsub|0>+\<Delta\>> (i.e.
  off-resonance). We also define <math|\<omega\><rsub|1>\<equiv\>\<gamma\>H<rsub|1>>
  and <math|\<omega\><rsub|0>=\<gamma\>H<rsub|0>>.

  <subsection|At the laboratory frame>

  Within the laboratory frame, Bloch equations are:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|\<mathd\>M<rsub|x>|\<mathd\>t>>|<cell|=>|<cell|\<gamma\><around*|(|M<rsub|y>*B<rsub|z>-M<rsub|z>*B<rsub|y>|)>-<frac|M<rsub|x>|T<rsub|2>>>>|<row|<cell|<frac|\<mathd\>M<rsub|y>|\<mathd\>t>>|<cell|=>|<cell|\<gamma\><around*|(|M<rsub|z>*B<rsub|x>-M<rsub|x>*B<rsub|z>|)>-<frac|M<rsub|y>|T<rsub|2>>>>|<row|<cell|<frac|\<mathd\>M<rsub|z>|\<mathd\>t>>|<cell|=>|<cell|\<gamma\><around*|(|M<rsub|x>*B<rsub|y>-M<rsub|y>*B<rsub|x>|)>-<frac|M<rsub|z>-M<rsub|0>|T<rsub|1>>>>>>
  </eqnarray*>

  where <math|B<rsub|x>=H<rsub|1>*cos \<omega\><rsub|z
  >t,B<rsub|y>=-H<rsub|1>*sin \<omega\><rsub|z>t,B<rsub|z>=H<rsub|0>>.

  And in matrix form:

  <\equation*>
    <frac|\<mathd\>|\<mathd\>t><bmatrix|<tformat|<table|<row|<cell|M<rsub|x>>>|<row|<cell|M<rsub|y>>>|<row|<cell|M<rsub|z>>>>>>=<bmatrix|<tformat|<table|<row|<cell|-<frac|1|T<rsub|2>>>|<cell|\<gamma\>B<rsub|z>>|<cell|-\<gamma\>B<rsub|y>>>|<row|<cell|-\<gamma\>B<rsub|z>>|<cell|-<frac|1|T<rsub|2>>>|<cell|\<gamma\>B<rsub|x>>>|<row|<cell|\<gamma\>B<rsub|y>>|<cell|-\<gamma\>B<rsub|x>>|<cell|-<frac|1|T<rsub|1>>>>>>><bmatrix|<tformat|<table|<row|<cell|M<rsub|x>>>|<row|<cell|M<rsub|y>>>|<row|<cell|M<rsub|z>>>>>>+<bmatrix|<tformat|<table|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|<frac|M<rsub|0>|T<rsub|1>>>>>>>
  </equation*>

  \;

  That can be solved numerically.

  <subsection|At the rotating frame>

  Within the rotating frame, Bloch equations become:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|\<mathd\>M<rprime|'><rsub|x>|\<mathd\>t>>|<cell|=>|<cell|\<gamma\>M<rsub|y><rprime|'>*h<rsub|0>-<frac|M<rsub|x><rprime|'>|T<rsub|2>>>>|<row|<cell|<frac|\<mathd\>M<rprime|'><rsub|y>|\<mathd\>t>>|<cell|=>|<cell|\<gamma\><around*|(|M<rsub|z><rprime|'>*H<rsub|1>-M<rsub|x><rprime|'>*h<rsub|0>|)>-<frac|M<rsub|y><rprime|'>|T<rsub|2>>>>|<row|<cell|<frac|\<mathd\>M<rprime|'><rsub|z>|\<mathd\>t>>|<cell|=>|<cell|-\<gamma\>M<rsub|y><rprime|'>*H<rsub|1>+<frac|M<rsub|0>-M<rsub|z><rprime|'>|T<rsub|1>>>>>>
  </eqnarray*>

  where <math|h<rsub|0>\<equiv\>H<rsub|0>+\<omega\><rsub|z>/\<gamma\>>. And
  in matrix form:

  <\equation*>
    <frac|\<mathd\>|\<mathd\>t><bmatrix|<tformat|<table|<row|<cell|M<rsub|x><rprime|'>>>|<row|<cell|M<rsub|y><rprime|'>>>|<row|<cell|M<rsub|z><rprime|'>>>>>>=<bmatrix|<tformat|<table|<row|<cell|-<frac|1|T<rsub|2>>>|<cell|\<gamma\>h<rsub|0>>|<cell|-\<gamma\>H<rsub|1>>>|<row|<cell|-\<gamma\>h<rsub|0>>|<cell|-<frac|1|T<rsub|2>>>|<cell|\<gamma\>H<rsub|1>>>|<row|<cell|-\<gamma\>H<rsub|1>>|<cell|-\<gamma\>H<rsub|1>>|<cell|-<frac|1|T<rsub|1>>>>>>><bmatrix|<tformat|<table|<row|<cell|M<rsub|x><rprime|'>>>|<row|<cell|M<rsub|y><rprime|'>>>|<row|<cell|M<rsub|z><rprime|'>>>>>>+<bmatrix|<tformat|<table|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|<frac|M<rsub|0>|T<rsub|1>>>>>>>
  </equation*>

  That can also be solved numerically and analytically.

  <subsection|Effects of parameters>

  <underline|Alternating <math|T<rsub|2>>:>

  For <math|T<rsub|1>=<SI|1|s>>, <math|\<omega\><rsub|1>=<SI|150|Hz>>,
  <math|\<Delta\>\<omega\>=<SI|2500|Hz><rsub|>>:

  <\with|par-mode|center>
    <tabular|<tformat|<cwith|1|1|1|1|par-mode|center>|<cwith|1|8|1|3|cell-halign|c>|<cwith|1|8|1|3|cell-valign|c>|<cwith|1|8|1|3|cell-tborder|1ln>|<cwith|1|8|1|3|cell-bborder|1ln>|<cwith|1|8|1|3|cell-lborder|1ln>|<cwith|1|8|1|3|cell-rborder|1ln>|<cwith|1|1|1|3|cell-background|pastel
    yellow>|<cwith|7|7|1|-1|cell-background|#faa>|<cwith|6|6|1|-1|cell-background|#faa>|<table|<row|<cell|<math|T<rsub|2>>>|<cell|<math|M<rsub|z><rsup|ss>>>|<cell|<math|t<rsup|ss>>
    (s)>>|<row|<cell|<SI|20|ms>>|<cell|0.848>|<cell|6.35>>|<row|<cell|<SI|5|ms>>|<cell|0.583>|<cell|5.18>>|<row|<cell|<SI|1|ms>>|<cell|0.244>|<cell|2.53>>|<row|<cell|<SI|500|<math|\<mu\>>s>>|<cell|0.186>|<cell|1.98>>|<row|<cell|<SI|100|<math|\<mu\>>s>>|<cell|0.321>|<cell|3.19>>|<row|<cell|<SI|10|<math|\<mu\>>s>>|<cell|0.817>|<cell|6.30>>>>>
  </with>

  Untill a certain threshold, <em|lowering <math|T<rsub|2>> decreases
  <math|M<rsub|z><rsup|ss>> and <math|t<rsup|ss>>>.

  When reaching the range of <math|T<rsub|2>> that represents solids(?), the
  effect is switched.\ 

  \;

  <underline|Alternating <math|T<rsub|1>>:>

  For <math|T<rsub|2>=<SI|2|ms>>, <math|\<omega\><rsub|1>=<SI|150|Hz>>,
  <math|\<Delta\>\<omega\>=<SI|2500|Hz><rsub|>>:

  <\with|par-mode|center>
    <tabular|<tformat|<cwith|1|1|1|1|par-mode|center>|<cwith|1|7|1|3|cell-halign|c>|<cwith|1|7|1|3|cell-valign|c>|<cwith|1|7|1|3|cell-tborder|1ln>|<cwith|1|7|1|3|cell-bborder|1ln>|<cwith|1|7|1|3|cell-lborder|1ln>|<cwith|1|7|1|3|cell-rborder|1ln>|<cwith|1|1|1|3|cell-background|pastel
    yellow>|<table|<row|<cell|<math|T<rsub|1>>>|<cell|<math|M<rsub|z><rsup|ss>>>|<cell|<math|t<rsup|ss>>
    (s)>>|<row|<cell|<SI|0.1|s>>|<cell|0.852>|<cell|0.83>>|<row|<cell|<SI|0.5|s>>|<cell|0.536>|<cell|2.62>>|<row|<cell|<SI|1|s>>|<cell|0.366>|<cell|3.58>>|<row|<cell|<SI|2|s>>|<cell|0.224>|<cell|4.38>>|<row|<cell|<SI|10|s>>|<cell|0.055>|<cell|5.34>>>>>
  </with>

  <em|Increasing <math|T<rsub|1>> decreases <math|M<rsub|z><rsup|ss>> and
  increases <math|t<rsup|ss>>.>

  <\underline>
    Alternating <math|\<omega\><rsub|1>>:
  </underline>

  For <math|T<rsub|1>=<SI|1|s>>, <math|T<rsub|2>=<SI|5|ms>>,
  <math|\<Delta\>\<omega\>=<SI|2500|Hz>>:

  <\with|par-mode|center>
    <tabular|<tformat|<cwith|1|1|1|1|par-mode|center>|<cwith|1|7|1|3|cell-halign|c>|<cwith|1|7|1|3|cell-valign|c>|<cwith|1|7|1|3|cell-tborder|1ln>|<cwith|1|7|1|3|cell-bborder|1ln>|<cwith|1|7|1|3|cell-lborder|1ln>|<cwith|1|7|1|3|cell-rborder|1ln>|<cwith|1|1|1|3|cell-background|pastel
    yellow>|<table|<row|<cell|<math|\<omega\><rsub|1>>
    (Hz)>|<cell|<math|M<rsub|z><rsup|ss>>>|<cell|<math|t<rsup|ss>>
    (s)>>|<row|<cell|50>|<cell|0.926>|<cell|6.18>>|<row|<cell|150>|<cell|0.583>|<cell|5.18>>|<row|<cell|300>|<cell|0.259>|<cell|2.69>>|<row|<cell|750>|<cell|0.053>|<cell|0.69>>|<row|<cell|1500>|<cell|0.014>|<cell|0.24>>>>>
  </with>

  <em|Increasing <math|\<omega\><rsub|1>> decreases <math|M<rsub|z><rsup|ss>>
  and <math|t<rsup|ss>>>.

  <\underline>
    Alternating <math|\<Delta\>\<omega\>>:
  </underline>

  For <math|T<rsub|1>=<SI|1|s>>, <math|T<rsub|2>=<SI|5|ms>>,
  <math|\<omega\><rsub|1>=<SI|150|Hz>>:

  <\with|par-mode|center>
    <tabular|<tformat|<cwith|1|1|1|1|par-mode|center>|<cwith|1|7|1|3|cell-halign|c>|<cwith|1|7|1|3|cell-valign|c>|<cwith|1|7|1|3|cell-tborder|1ln>|<cwith|1|7|1|3|cell-bborder|1ln>|<cwith|1|7|1|3|cell-lborder|1ln>|<cwith|1|7|1|3|cell-rborder|1ln>|<cwith|1|1|1|3|cell-background|pastel
    yellow>|<table|<row|<cell|<math|\<Delta\>\<omega\>>
    (Hz)>|<cell|<math|M<rsub|z><rsup|ss>>>|<cell|<math|t<rsup|ss>>
    (s)>>|<row|<cell|0>|<cell|0.0088>|<cell|0.13>>|<row|<cell|250>|<cell|0.022>|<cell|0.31>>|<row|<cell|1000>|<cell|0.188>|<cell|2.03>>|<row|<cell|2500>|<cell|0.583>|<cell|5.18>>|<row|<cell|5000>|<cell|0.848>|<cell|6.35>>>>>
  </with>

  Increasing <math|\<Delta\>\<omega\>> increases <math|M<rsub|z><rsup|ss>>
  and <math|t<rsup|ss>>.

  <underline|Optimal saturation pulse>

  The saturation pulse is defined such that at <math|t<rsup|ss>>,
  <math|M<rsub|z><rsup|ss>\<rightarrow\>0>.

  We want to find the <with|font-series|bold|shortest> saturation pulse
  obtainable.

  Turns out it's always on resonance.
</body>

<\initial>
  <\collection>
    <associate|font-base-size|12>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|1.2|2>>
    <associate|auto-4|<tuple|1.3|2>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Different
      representations of Bloch equations>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>At the laboratory frame
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>At the rotating frame
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|1.3<space|2spc>Effects of parameters
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>
    </associate>
  </collection>
</auxiliary>