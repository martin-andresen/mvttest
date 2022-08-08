{smcl}
{cmd:help mvttest}{e})}
{hline}

{title:Title}

{p2colset 5 15 17 2}{...}
{p2col:{cmd:mvttest} {hline 2}}Tests for instrument validity with binarized multivalued treatments{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 11 2}
{cmd:mvttest} {depvar:} {varname:_iv} [{indepvars:}] {ifin}, {opt jstar:(#)} [{opt cmi:opts(string)} {opt graph:opts(string)} {opt noplot} {opt nocmi}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}

{syntab:Options}
{synopt:{opt jstar:{#}}}is required, and specifies the cutoff value for constructing the binarized treatment{p_end}
{synopt:{opt cmi:opts{string}}}estimation options that are passed directly to the {cmd:cmi_test} command{p_end}
{synopt:{opt graph:opts{string}}}options for the graph. Overrides default looks. See {helpb twoway_options}{p_end}
{synopt:{opt nocmi}}turns off condititional moment inequality tests.{p_end}
{synopt:{opt noplot}}Does not display graph.{p_end}
{synopt:{opt cluster(clustvar)}}Cluster standard errors on {it:clustvar}.{p_end}

{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:mvttest} performs F-tests and conditional moment inequality tests of instrument validity in settings where 
a multivalued treatment is binarized, see Andresen and Huber (2020).


{marker remarks}{...}
{title:Remarks}

{pstd}In a setting where Z is a valid binary instrument for D, a multivalued 
treatment (conditional on X), {cmd:mvttest} performs various tests of the validity 
of Z as an instrument for D*, a binarized treatment constructed so that 1{D>=j*}, 
where j* is a threshold value.{p_end}

{pstd}{cmd:mvvtest} first performs a multivariate regression of all possible binary 
treatment indicators on the instrument using {cmd:reghdfe}. If [{indepvars:}] is specified, these
coefficients are estimated separately in cells of the independent variables.{p_end}

{pstd}Unless "noplot" is specified results are presented graphically. If no controls are used, the 
coefficients are plotted across j. If controls are used, mvttest instead plots the maximum violation
of Assumption 3 across cells of X for each level j, which amounts to b_j-b_{j+1} for j<j* and 
b_{j+1}-b_j for j>=j*. For comparison, the violations in a model with no controls is also plotted, 
allowing the user to see whether violations in some subgroups of X are averaged over.{p_end}

{pstd}The paper provides testable conditions for the b_j coefficients, where b_j is the 
first stage effect of Z on D* using j as the threshold. The special case in Assumption 3*
require all b_j to be 0 except b_j*, while Assumption 4 requires all b_j to be the same.
Both these conditions can be tested using G-tests.{p_end}

{pstd}More generally, a necessary condition for instrument validity in this setting
(Assumption 3) is given by{p_end}

{p 10} b_{j+1}>=b_j for all j<j*{p_end}
{p 10} b_{j+1}<=b_j for all j>=j*{p_end}

{pstd}This is tested in the conditional moment equality/inequality framework of Andrews
et al. (2013), implemented using {helpb cmi_test}.{p_end}

{marker examples}{...}
{title:Examples}

{pstd}
Download data from NLSYM, used by Card (1995) {p_end}
{phang2}{stata "use http://www.stata.com/data/jwooldridge/eacsap/card.dta, clear": {cmd:. use "http//www.stata.com/data/jwooldridge/eacsap/card.dta", clear}} ({it:click to run}){p_end}

{pstd}
Perform tests of instrument validity for the nearc4 instrument, instrumenting for treatment defined as years of education >=16{p_end}
{phang2}{stata mvttest educ nearc4, jstar(16): {cmd:. mvttest educ nearc4, jstar(16)}} ({it:click to run}){p_end}

{pstd}
Looking for violations within cells of smsa and region:{p_end}
{phang2}{stata mvttest educ nearc4 smsa region, jstar(16): {cmd:. mvttest educ nearc4 smsa region, jstar(16)}} ({it:click to run}){p_end}


{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:mvttest} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(levels_X)}}Number of cells of indepvars{p_end}
{synopt:{cmd:e(F3star)}}F-stat for the test of Assumption 3*{p_end}
{synopt:{cmd:e(p_val3star)}}p-value for the F-test of Assumption 3*{p_end}
{synopt:{cmd:e(df3star)}}Constraints tested for the F-test of Assumption 3*{p_end}
{synopt:{cmd:e(df_r3star)}}Degrees of freedom for the F-test of Assumption 3*{p_end}
{synopt:{cmd:e(F4)}}F-stat for the test of Assumption 4{p_end}
{synopt:{cmd:e(p_val4)}}p-value for the F-test of Assumption 4{p_end}
{synopt:{cmd:e(df4)}}Constraints tested for the F-test of Assumption 4{p_end}
{synopt:{cmd:e(df_r4)}}Degrees of freedom for the F-test of Assumption 4{p_end}
{synopt:{cmd:e(N_ineq)}}number of inequalities tested (number of beta_j+1-beta_j pairs){p_end}
{synopt:{cmd:e(N_ineqcells)}}number of inequalities by cells tested{p_end}'
{synopt:{cmd:e(cmi_stat)}}The test statistic for the CMI-test of Assumption 2{p_end}
{synopt:{cmd:e(cmi_cv01)}}1% critical value for the CMI-test{p_end}
{synopt:{cmd:e(cmi_cv05)}}5% critical value for the CMI-test{p_end}
{synopt:{cmd:e(cmi_cv10)}}10% critical value for the CMI-test{p_end}
{synopt:{cmd:e(cmi_pval)}}p-value for the F-test of Assumption 2{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:mvttest}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}threshold{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}; {cmd:V}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector for b_j{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:b}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{p2colreset}{...}



{marker references}{...}
{title:References}

{phang}
Andresen, M. E. and M. Huber 2020. "Instrument-based estimation with binarized treatments: Issues and tests for the exclusion restriction", working paper

{phang}
Andrews, D., and X. Shi. 2013. "Inference based on conditional moment inequalities", Econometrica, Vol.81, No. 2, 609-666.

{phang}
Andrews, D., W. Kim and X. Shi. 2017. "Commands for testing conditional moment inequalities and equalities", Stata Journal, Vol.17, No. 1, 56-72.

{phang}
Card, D. 1995. "Using geographic variation in college proximity to estimate the return to schooling", In "Aspects of labour market behavior: Essays in honour of John Vanderkamp", pp. 201-222, University of Toronto Press, Toronto.


{marker Author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}Statistics Norway{p_end}
{pstd}Oslo, NO{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}

{marker also_see}{...}
{title:Also see}
{p 7 14 2}Help:  {helpb cmi_test} and {helpb reghdfe}.{p_end}
