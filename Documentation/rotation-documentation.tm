<TeXmacs|2.1>

<style|<tuple|generic|SIUnits-simple>>

<\body>
  <doc-data|<doc-title|Rotation Documentation>>

  <section|Paradigm and Expectations><marginal-note|normal|c|27/10/22>

  <em|Final goal>: How does CEST (Chemical Exchange Saturation Transfer)
  behave when the subject material diffuses between solute/liquid and solid
  phases?

  <em|Initial goal>: At the liquid phase, given <math|T<rsub|1> > and
  <math|T<rsub|2>>, when does the dillute spin population reaches saturation?

  <em|Intermediate goal>: When saturation is reached (and the net
  magnetization of the dillute pool approaches zero), how are the spins
  transfered to another pool?

  <section|Initial Assignments><marginal-note|normal|c|31/10/22>

  Roadblock:

  <\itemize-dot>
    <item>Understand the Bloch equations. (How the magnetization is
    obtained.) <with|color|blue|\<checkmark\>>

    <item>Break down Assaf Tal's Matlab code (spintool) t<strong|>o
    understand how the magnetization is calculated.

    <item>Run the code and simulate saturation pulses to gain insight on how
    the saturation is affected by <math|T<rsub|1>,\<space\>T<rsub|2>> and the
    saturation pulse amplitude.
  </itemize-dot>

  Alternative approach:

  <\enumerate>
    <item>Understand the Bloch equations. <with|color|blue|\<checkmark\>>

    <item>Write a code which solves the Bloch equations for a vector of the
    form \ <rigid|<math|<wide|M|\<vect\>>=<around*|[|M<rsub|x>,\<space\>M<rsub|y>,\<space\>M<rsub|z>|]>>.>
    It should act on an initial magnetization vector
    <rigid|<math|<wide|M|\<vect\>>=M<rsub|0><around*|[|0,0,1|]>>> and use a
    matrix operator which includes the static external field
    <math|H<rsub|0>>, the pertubating field <math|H<rsub|1>> rotating at
    frequency <math|\<omega\><rsub|z>>, the relaxation constants
    <math|T<rsub|1>,\<space\>T<rsub|2>>, and the offset from resonance.
    <with|color|blue|\<checkmark\>>\ 

    <item>Verify that the result commutes with the analytical model
    prediction.

    <item>Check how <math|M<rsub|z><rsup|ss>> and <math|t<rsup|sat>> is
    affected by <math|\<omega\><rsub|1>> and <math|\<Delta\>\<omega\>> for
    different pairs of <math|T<rsub|1>,T<rsub|2>>.
  </enumerate>

  After which,

  <\itemize-dot>
    <item>Understand Bloch-McConnell equation and how CEST works.

    <item>Look up how CEST in pure liquid state is affected by the system
    parameters.
  </itemize-dot>

  \;

  \;
</body>

<\initial>
  <\collection>
    <associate|font-base-size|12>
    <associate|page-medium|paper>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Paradigm
      and Expectations> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Initial
      Assignments> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>