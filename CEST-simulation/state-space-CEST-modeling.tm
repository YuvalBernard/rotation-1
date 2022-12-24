<TeXmacs|2.1.2>

<style|generic>

<\body>
  <doc-data|<doc-title|State space modeling of Bloch-McConnell equations to
  calculate Z and <math|MTR<rsub|Assym>> spectra>>

  Given the set of differential equations:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|M|\<dot\>><rsub|x><rsup|a>>|<cell|=>|<cell|<around*|(|-R<rsub|2><rsup|a>+k<rsub|a>|)>*M<rsub|x><rsup|a>+\<Delta\>\<omega\><rsub|a>*M<rsub|y><rsup|a>+k<rsub|b>*M<rsub|x><rsup|b>>>|<row|<cell|<wide|M|\<dot\>><rsub|y><rsup|a>>|<cell|=>|<cell|-\<Delta\>\<omega\><rsub|a><rsub|><rsup|>*M<rsub|x><rsup|a>-<around*|(|R<rsub|2><rsup|a>+k<rsub|a>|)>*M<rsub|y><rsup|a>-\<omega\><rsub|1>*M<rsub|z><rsup|a>+k<rsub|b>*M<rsub|y><rsup|b>>>|<row|<cell|<wide|M|\<dot\>><rsub|z><rsup|a>>|<cell|=>|<cell|\<omega\><rsub|1>*M<rsub|y><rsup|a>-<around*|(|R<rsub|1><rsup|a>+k<rsub|a>|)>*M<rsub|z><rsup|a>+k<rsub|b>*M<rsub|z><rsup|b>+R<rsub|1><rsup|a>*M<rsub|z><rsup|a,0>>>|<row|<cell|<wide|M|\<dot\>><rsub|x><rsup|b>>|<cell|=>|<cell|k<rsub|a>*M<rsub|x><rsup|a>-<around*|(|R<rsub|2><rsup|b>+k<rsub|b>|)>*M<rsub|x><rsup|b>+\<Delta\>\<omega\><rsub|><rsup|><rsub|b>*M<rsub|y><rsup|b>>>|<row|<cell|<wide|M|\<dot\>><rsub|y><rsup|b>>|<cell|=>|<cell|k<rsub|a>*M<rsub|y><rsup|a>-\<Delta\>\<omega\><rsup|><rsub|b>*M<rsub|x><rsup|b>-<around*|(|R<rsub|2><rsup|b>+k<rsub|b>|)>*M<rsub|y><rsup|b>-\<omega\><rsub|1>*M<rsub|z><rsup|b>>>|<row|<cell|<wide|M|\<dot\>><rsub|z><rsup|b>>|<cell|=>|<cell|k<rsub|a>*M<rsub|z><rsup|a>+\<omega\><rsub|1>*M<rsub|y><rsup|b>-<around*|(|R<rsub|1><rsup|b>+k<rsub|b>|)>*M<rsub|z><rsup|b>+R<rsub|1><rsup|b>*M<rsub|z><rsup|b,0>>>>>
  </eqnarray*>

  The state variables are

  <\equation*>
    <wide|x|\<vect\>>=<choice|<tformat|<table|<row|<cell|x<rsub|1>=M<rsub|x><rsup|a>>>|<row|<cell|x<rsub|2>=M<rsub|y><rsup|a>>>|<row|<cell|x<rsub|3>=M<rsub|z><rsup|a>>>|<row|<cell|x<rsub|4>=M<rsub|x><rsup|b>>>|<row|<cell|x<rsub|5>=M<rsub|y><rsup|b>>>|<row|<cell|x<rsub|6>=M<rsub|z><rsup|b>>>>>>,<application-space|1em><wide|x|\<vect\>><around*|(|0|)>=<bmatrix|<tformat|<table|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|M<rsub|z><rsup|a,0>>>|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|M<rsub|z><rsup|b,0>>>>>>
  </equation*>

  The state matrix <math|A> and input matrix <math|B> are:

  <\equation*>
    A=<bmatrix|<tformat|<table|<row|<cell|-R<rsub|2><rsup|a>+k<rsub|a>>|<cell|\<Delta\>\<omega\><rsub|a>>|<cell|0>|<cell|k<rsub|b>>|<cell|0>|<cell|0>>|<row|<cell|-\<Delta\>\<omega\><rsub|a>>|<cell|-R<rsub|2><rsup|a>+k<rsub|a>>|<cell|\<omega\><rsub|1>>|<cell|0>|<cell|k<rsub|b>>|<cell|0>>|<row|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|-<around*|(|R<rsub|1><rsup|a>+k<rsub|a>|)>>|<cell|0>|<cell|0>|<cell|k<rsub|b>>>|<row|<cell|k<rsub|a>>|<cell|0>|<cell|0>|<cell|-<around*|(|R<rsub|2><rsup|b>+k<rsub|b>|)>>|<cell|\<Delta\>\<omega\><rsub|b>>|<cell|0>>|<row|<cell|0>|<cell|k<rsub|a>>|<cell|0>|<cell|-\<Delta\>\<omega\><rsub|b>>|<cell|-<around*|(|R<rsub|2><rsup|b>+k<rsub|b>|)>>|<cell|\<omega\><rsub|1>>>|<row|<cell|0>|<cell|0>|<cell|k<rsub|a>>|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|-<around*|(|R<rsub|1><rsup|b>+k<rsub|b>|)>>>>>>
  </equation*>

  <\equation*>
    <wide|b|\<vect\>>=<bmatrix|<tformat|<table|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|R<rsub|1><rsup|a>*<frac|M<rsub|z><rsup|a,0>|M<rsub|z><rsup|b,0>>>>|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|R<rsub|1><rsup|b>*>>>>>
  </equation*>

  and the input is <math|M<rsub|z><rsup|b,0>>. As such,

  <\equation*>
    <wide|<wide|x|\<dot\>>|\<vect\>>=A*<wide|x|\<vect\>>+<wide|b|\<vect\>>*u
  </equation*>

  As we are interested in <math|M<rsub|z><rsup|a>> only, we define the output
  matrix <math|C> and the direct transmittion matrix <math|D> such that:

  <\equation*>
    C=<bmatrix|<tformat|<table|<row|<cell|0>|<cell|0>|<cell|1>|<cell|0>|<cell|0>|<cell|0>>>>>,<application-space|1em>D=0
  </equation*>

  This makes it that <math|y=M<rsub|z><rsup|a>> is the output.

  <\equation*>
    y=C*<wide|x|\<vect\>>+D*u
  </equation*>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>