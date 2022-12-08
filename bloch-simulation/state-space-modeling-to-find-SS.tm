<TeXmacs|2.1>

<style|generic>

<\body>
  <doc-data|<doc-title|State space modeling of bloch equations to find steady
  state magnetization>>

  Given the set of differential equations:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|M|\<dot\>><rsub|x>>|<cell|=>|<cell|\<Delta\>\<omega\>*M<rsub|y>-R<rsub|2>*M<rsub|x>>>|<row|<cell|<wide|M|\<dot\>><rsub|y>>|<cell|=>|<cell|-\<Delta\>\<omega\>*M<rsub|x>-R<rsub|2>*M<rsub|y>+\<omega\><rsub|1>*M<rsub|z>>>|<row|<cell|<wide|M|\<dot\>><rsub|z>>|<cell|=>|<cell|-\<Delta\>\<omega\><rsub|1>*M<rsub|y>-R<rsub|1><around*|(|M<rsub|z>-M<rsub|z><rsup|0>|)>>>>>
  </eqnarray*>

  Our <strong|state variables> are\ 

  <\equation*>
    <choice|<tformat|<table|<row|<cell|x<rsub|1>=M<rsub|x>>>|<row|<cell|x<rsub|2>=M<rsub|y>>>|<row|<cell|x<rsub|3>=M<rsub|z>>>|<row|<cell|x<rsub|4>=1>>>>>
  </equation*>

  If we write the set of equation in matrix format, we get the following:

  <\equation*>
    <around*|{|<tabular|<tformat|<table|<row|<cell|<wide|x|\<dot\>><rsub|1>>>|<row|<cell|<wide|x|\<dot\>><rsub|2>>>|<row|<cell|<wide|x|\<dot\>><rsub|3>>>|<row|<cell|<wide|x|\<dot\>><rsub|4>>>>>>|}>=<bmatrix|<tformat|<table|<row|<cell|-R<rsub|2>>|<cell|\<Delta\>\<omega\>>|<cell|0>|<cell|0>>|<row|<cell|-\<Delta\>\<omega\>>|<cell|-R<rsub|1>>|<cell|\<omega\><rsub|1>>|<cell|0>>|<row|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|-R<rsub|1>>|<cell|R<rsub|1>M<rsub|z><rsup|0>>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>>>>>\<cdot\><around*|{|<tabular|<tformat|<table|<row|<cell|x<rsub|1>>>|<row|<cell|x<rsub|2>>>|<row|<cell|x<rsub|3>>>|<row|<cell|x<rsub|4>>>>>>|}>+<bmatrix|<tformat|<table|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|0>>>>>\<cdot\><around*|{|<tabular|<tformat|<table|<row|<cell|u>>>>>|}>
  </equation*>

  The <strong|state matrix> A and <strong|input matrix> B are

  <\equation*>
    A=<bmatrix|<tformat|<table|<row|<cell|-R<rsub|2>>|<cell|\<Delta\>\<omega\>>|<cell|0>|<cell|0>>|<row|<cell|-\<Delta\>\<omega\>>|<cell|-R<rsub|1>>|<cell|\<omega\><rsub|1>>|<cell|0>>|<row|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|-R<rsub|1>>|<cell|R<rsub|1>*M<rsub|z><rsup|0>>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>>>>>,B=<bmatrix|<tformat|<table|<row|<cell|0>>>>>
  </equation*>

  To identify the <strong|output> and <strong|direct transmission matrices>,
  we need to decide which is the measurable output. We want to measure
  <math|M<rsub|z>>. So we set

  <\equation*>
    y=x<rsub|3>
  </equation*>

  If we write the output equation in matrix format, we get

  <\equation*>
    <around*|{|<tabular|<tformat|<table|<row|<cell|y>>>>>|}>=<around*|[|<tabular|<tformat|<table|<row|<cell|0>|<cell|0>|<cell|1>|<cell|0>>>>>|]>\<cdot\><around*|{|<tabular|<tformat|<table|<row|<cell|x<rsub|1>>>|<row|<cell|x<rsub|2>>>|<row|<cell|x<rsub|3>>>|<row|<cell|x<rsub|4>>>>>>|}>+<around*|[|<tabular|<tformat|<table|<row|<cell|0>>>>>|]>\<cdot\><around*|{|<tabular|<tformat|<table|<row|<cell|u>>>>>|}>
  </equation*>

  Which gives the <strong|output matrix> <math|C> and the <strong|direct
  transmittion matrix> <math|D>:

  <\equation*>
    C=<around*|[|<tabular|<tformat|<table|<row|<cell|0>|<cell|0>|<cell|1>|<cell|0>>>>>|]>,<application-space|1em>D=<around*|[|<tabular|<tformat|<table|<row|<cell|0>>>>>|]>
  </equation*>

  Now, we want to use the <verbatim|dcgain> function to infer what the steady
  state value is.
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>