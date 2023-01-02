<TeXmacs|2.1.2>

<style|generic>

<\body>
  <doc-data|<doc-title|3 pool CEST model>>

  Consider a system where a large pool <math|A> can exchange protons with two
  solute pools: <math|B> and <math|C>.

  <\equation*>
    A<long-arrow|\<rubber-rightleftharpoons\>|k<rsub|ab>|k<rsub|ba>><rsub|\<nosymbol\>><rsup|\<nosymbol\>><rsub|\<nosymbol\>>B
  </equation*>

  <\equation*>
    A<long-arrow|\<rubber-rightleftharpoons\>|k<rsub|ac>|k<rsub|ca>>C
  </equation*>

  Then:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|\<mathd\><around*|[|A|]>|\<mathd\>t>>|<cell|=>|<cell|-k<rsub|ab>*<around*|[|A|]>+k<rsub|ba>*<around*|[|B|]>-k<rsub|ac>*<around*|[|A|]>+k<rsub|ca>*<around*|[|C|]>>>|<row|<cell|<frac|\<mathd\><around*|[|B|]>|\<mathd\>t>>|<cell|=>|<cell|k<rsub|ab>*<around*|[|A|]>-k<rsub|ba>*<around*|[|B|]>>>|<row|<cell|<frac|\<mathd\><around*|[|C|]>|\<mathd\>t>>|<cell|=>|<cell|k<rsub|ac>*<around*|[|A|]>-k<rsub|ca>*<around*|[|C|]>>>>>
  </eqnarray*>

  And in matrix form:

  <\equation*>
    <frac|\<mathd\>|dt><bmatrix|<tformat|<table|<row|<cell|<around*|[|A|]>>>|<row|<cell|<around*|[|B|]>>>|<row|<cell|<around*|[|C|]>>>>>>=<bmatrix|<tformat|<table|<row|<cell|-<around*|(|k<rsub|ab>+k<rsub|ac>|)>>|<cell|k<rsub|ba>>|<cell|k<rsub|ca>>>|<row|<cell|k<rsub|ab>>|<cell|-k<rsub|ba>>|<cell|0>>|<row|<cell|k<rsub|ac>>|<cell|0>|<cell|-k<rsub|ca>>>>>>
  </equation*>

  The magnetization space vector representing magnatization of the pools is
  given by:

  <\equation*>
    <bmatrix|<tformat|<table|<row|<cell|1>>|<row|<cell|1>>|<row|<cell|1>>>>>\<otimes\><bmatrix|<tformat|<table|<row|<cell|M<rsub|x>>>|<row|<cell|M<rsub|y>>>|<row|<cell|M<rsub|z>>>>>>=<bmatrix|<tformat|<table|<row|<cell|M<rsub|x><rsup|A>>>|<row|<cell|M<rsub|y><rsup|A>>>|<row|<cell|M<rsub|z><rsup|A>>>|<row|<cell|M<rsub|x><rsup|B>>>|<row|<cell|M<rsub|y><rsup|B>>>|<row|<cell|M<rsub|z><rsup|B>>>|<row|<cell|M<rsub|x><rsup|C>>>|<row|<cell|M<rsub|y><rsup|C>>>|<row|<cell|M<rsub|z><rsup|C>>>>>>
  </equation*>

  A new exchange matrix in the basis of the magnetization space is calculated
  by taking the direct product of the exchange matrix with the identity
  operator in the same dimension as the magnetization space:

  <\equation*>
    <bmatrix|<tformat|<table|<row|<cell|-<around*|(|k<rsub|ab>+k<rsub|ac>|)>>|<cell|k<rsub|ba>>|<cell|k<rsub|ca>>>|<row|<cell|k<rsub|ab>>|<cell|-k<rsub|ba>>|<cell|0>>|<row|<cell|k<rsub|ac>>|<cell|0>|<cell|-k<rsub|ca>>>>>>\<otimes\><bmatrix|<tformat|<table|<row|<cell|1>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|1>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|1>>>>>=\<cdots\>
  </equation*>

  <\equation>
    \<cdots\>=<bmatrix|<tformat|<table|<row|<cell|-<around*|(|k<rsub|ab>+k<rsub|ac>|)>>|<cell|k<rsub|ba>>|<cell|k<rsub|ca>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|0>>|<row|<cell|k<rsub|ab>>|<cell|-k<rsub|ba>>|<cell|0>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|k<rsub|ac>>|<cell|0>|<cell|-k<rsub|ca>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|-<around*|(|k<rsub|ab>+k<rsub|ac>|)>>|<cell|k<rsub|ba>>|<cell|k<rsub|ca>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|k<rsub|ab>>|<cell|-k<rsub|ba>>|<cell|0>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|k<rsub|ac>>|<cell|0>|<cell|-k<rsub|ca>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|-<around*|(|k<rsub|ab>+k<rsub|ac>|)>>|<cell|k<rsub|ba>>|<cell|k<rsub|ca>>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|k<rsub|ab>>|<cell|-k<rsub|ba>>|<cell|0>>|<row|<cell|0>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|k<rsub|ac>>|<cell|0>|<cell|-k<rsub|ca>>>>>>
  </equation>

  Likewise, the matrix describing coherent and incoherent magnetization
  interactions can be recast in a similar fashion to give

  <\equation*>
    <bmatrix|<tformat|<table|<row|<cell|1>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|1>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|1>>>>>\<otimes\><bmatrix|<tformat|<table|<row|<cell|R<rsub|2>>|<cell|\<Delta\>*\<omega\>>|<cell|0>>|<row|<cell|-\<Delta\>*\<omega\>>|<cell|R<rsub|2>>|<cell|\<omega\><rsub|1>>>|<row|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|R<rsub|1>>>>>>=\<cdots\>
  </equation*>

  <\equation>
    \<cdots\>=<bmatrix|<tformat|<table|<row|<cell|R<rsup|A><rsub|2>>|<cell|\<Delta\>*\<omega\><rsup|A>>|<cell|0>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|0>>|<row|<cell|-\<Delta\>*\<omega\><rsup|A>>|<cell|R<rsub|2><rsup|A>>|<cell|\<omega\><rsub|1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|R<rsub|1><rsup|A>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|R<rsup|B><rsub|2>>|<cell|\<Delta\>*\<omega\><rsup|B>>|<cell|0>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|-\<Delta\>*\<omega\><rsup|B>>|<cell|R<rsub|2><rsup|B>>|<cell|\<omega\><rsub|1>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|R<rsub|1><rsup|B>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|R<rsup|C><rsub|2>>|<cell|\<Delta\>*\<omega\><rsup|C>>|<cell|0>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|-\<Delta\>*\<omega\><rsup|C>>|<cell|R<rsub|2><rsup|C>>|<cell|\<omega\><rsub|1>>>|<row|<cell|0>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|0>|<cell|-\<omega\><rsub|1>>|<cell|R<rsub|1><rsup|C>>>>>>
  </equation>

  Add (1) and (2) to get:

  <\equation*>
    <frac|\<mathd\>|\<mathd\>t><bmatrix|<tformat|<table|<row|<cell|M<rsub|x><rsup|A>>>|<row|<cell|M<rsub|y><rsup|A>>>|<row|<cell|M<rsub|z><rsup|A>>>|<row|<cell|M<rsub|x><rsup|B>>>|<row|<cell|M<rsub|y><rsup|B>>>|<row|<cell|M<rsub|z><rsup|B>>>|<row|<cell|M<rsub|x><rsup|C>>>|<row|<cell|M<rsub|y><rsup|C>>>|<row|<cell|M<rsub|z><rsup|C>>>>>>=-<around*|[|<around*|(|1|)>+<around*|(|2|)>|]>*<bmatrix|<tformat|<table|<row|<cell|M<rsub|x><rsup|A>>>|<row|<cell|M<rsub|y><rsup|A>>>|<row|<cell|M<rsub|z><rsup|A>>>|<row|<cell|M<rsub|x><rsup|B>>>|<row|<cell|M<rsub|y><rsup|B>>>|<row|<cell|M<rsub|z><rsup|B>>>|<row|<cell|M<rsub|x><rsup|C>>>|<row|<cell|M<rsub|y><rsup|C>>>|<row|<cell|M<rsub|z><rsup|C>>>>>>+<bmatrix|<tformat|<table|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|R<rsub|1><rsup|A>*M<rsub|z,eq><rsup|A>>>|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|R<rsub|1><rsup|B>*M<rsub|z,eq><rsup|B>>>|<row|<cell|0>>|<row|<cell|0>>|<row|<cell|R<rsub|1><rsup|C>*M<rsub|z,eq><rsup|C>>>>>>
  </equation*>

  This system of equations can be easily solved<text-dots>
</body>

<\initial>
  <\collection>
    <associate|font-base-size|12>
    <associate|page-medium|paper>
  </collection>
</initial>