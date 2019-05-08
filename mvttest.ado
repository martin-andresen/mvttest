*! mvttest v 1.2 31aug2018
*! Author: Martin E. Andresen
*! For "Testing Binarized Treatments: Issues and tests for the exclusion restriction, joint with Martin Huber

cap program drop mvttest
{
	program define mvttest, eclass
		version 13.0
		syntax varlist(min=2) [fweight pweight], jstar(numlist max=1 min=1) [ /*
				*/cmi_opts(string) 		/* Any other options for the cmi_test
				*/graph_opts(string)  	/* Options for the twoway graph. Overrides default
				*/noPlot				/* Do not plot graph 
			/*	*/ftest					/* Perform an additional F-test of the non-complying beta pairs */
				*/keepsingletons		/* Keeps observations in singleton groups
				*/nocmi					/* Do not perform CMI-tests
				*/]
	
	marksample touse
	preserve
	
	keep `touse' `varlist'
	qui {
		//Parse and check input
		tokenize `varlist'
		loc D "`1'"
		mac shift
		loc Z "`1'"
		mac shift
		loc X `*'
		
		if strpos("`cmi_opts'","cmv")>0 loc stattype "CvM"
		else loc stattype KS
		
		_fv_check_depvar `Z'
		tab `Z' if `touse'
		if `r(r)'==2 loc c 1.
		else loc c c.
		
		foreach prog in reghdfe cmi_test fastreshape gtools {
			capture which `prog'
			if _rc!=0 {
				noi di in red "mvttest uses `prog', which is not installed. Install using ssc install `prog'"
				exit 
				}
			}

		tempvar id regno
		gen `id'=_n
		
		foreach var in `X' {
			loc xabs `xabs'#`var'
			}
		
		keep if `touse'
		sort `D'
		
		if "`X'"!="" {
			tempvar N_g group
			egen `group'=group(`X')
			if "`keepsingletons'"!="" {
				loc totcount=0
				while `count'!=0 {
					bys `group': gen `N_g'=_N if `touse'
					count if `N_g'==1&`touse'
					loc count=r(N)
					loc totcount=`totcount'+`count'
					replace `touse'=`touse'&`N_g'>1
					}
				if `totcount'>0 noi di as text "Note: Dropped `totcount' observations in singleton groups."
				}
			}
			
		levelsof `D' if `touse', local(values)
		
		loc no=0
		foreach j in `values' {
			loc ++no
			if `no'>1 gen `D'`j'=`D'>=`j'
			}
			
		su `D' if `touse'
		loc J=r(max)
		if `J'<=2 {
			noi di in red "Treatment `D' not multivalued"
			exit
			}
		loc jzero=r(min)
		
		drop `D'
		
		reshape long `D' `d', i(`id') j(`regno')
	
		reghdfe `D' `c'`Z'#`regno' if `touse', absorb(`regno'`xabs') vce(cluster `id') `keepsingletons' nocons
		loc N=e(N)
		
		tempname b V
		mat `b'=e(b)
		mat `V'=e(V)
		
		foreach j in `values' {
			capture di _b[`c'`Z'#`j'.`regno']
			if _rc==0  {
				if `j'<`jstar' loc colnames `colnames' Below:`j'
				if `j'>=`jstar' loc colnames `colnames' Above:`j'
				loc values2 `values2' `j'
				}
			
			}
		loc values `values2'
		
		mat colnames `b'=`colnames'
		mat colnames `V'=`colnames'
		mat rownames `V'=`colnames'
		
		loc no=0
		tokenize `values'
		foreach j in `values' {
			loc ++no
			capture di _b[`c'`Z'#`j'.`regno']
			if _rc==0 {
				if `j'<`jstar' {
					loc teststring4 `teststring4' `c'`Z'#`j'.`regno'
					loc teststring5 `teststring5' `c'`Z'#`j'.`regno' =
					
					/* if "`ftest'"!="" {
						if _b[`c'`Z'#`j'.`regno']>_b[`c'`Z'#``=`no'+1''.`regno'] {
							loc teststring3 `teststring3' (`c'`Z'#`j'.`regno' = `c'`Z'#``=`no'+1''.`regno')
							loc test3desc `test3desc' `j'=``=`no'+1'',
							}
						}
					*/
					
					}
				if `j'>=`jstar' {
					if "`thresholdj'"=="" {
						loc thresholdj=`j'
						}
					else {
						loc teststring5 `teststring5' `c'`Z'#`j'.`regno' =	
						loc teststring4 `teststring4' `c'`Z'#`j'.`regno'
						}
					if "`ftest'"!="" {
						if `j'<`J' {
							if _b[`c'`Z'#`j'.`regno']<_b[`c'`Z'#``=`no'+1''.`regno'] {
								loc teststring2 `teststring2' (`c'`Z'#`j'.`regno' = `c'`Z'#``=`no'+1''.`regno')
								loc test2desc `test2desc' `j'=``=`no'+1'',
								}
							}
						}
					
					}
				}
			}
		
		/*if "`ftest'"!="" {
			if "`teststring3'"!="" {
				noi test `teststring3', df(`=`J'-1')
				loc F3=r(F)
				loc p3=r(p)
				}
			}
		*/

		if "`teststring4'"!="" {
			test `teststring4'
			loc F4=r(F)
			loc p4=r(p)
			}
		
		if "`teststring5'"!="" {
			test `teststring5' `c'`Z'#`thresholdj'.`regno'
			loc F5=r(F)
			loc p5=r(p)
			}
		
		parmest, norestore
		loc no=0
		tempvar jno
		gen `jno'=_n
		foreach j in `values' {
			loc ++no
			loc xlabel `xlabel' `no' "`j'"
			if `j'>=`jstar'&"`xline'"=="" loc xline=`=`no'-0.5'
			}
		
		if `J'>20 loc small ,labsize(vsmall) angle(90)
		if "`plot'"!="noplot" twoway (bar estimate `jno', color(navy) lcolor(white) lwidth(medium)) (rcap min95 max95 `jno', color(navy)) ///
			,  scheme(s1color) graphregion(color(white)) plotregion(lcolor(black)) ///
			xtitle("Treatment at least") legend(off) title("First stage effect using various thresholds") ///
			xlabel(`xlabel' `small') xline(`xline', lcolor(black) lpattern(dash)) `graph_opts'
		
		restore
		
		
		
		
		//CMI-test
		if "`cmi'"!="nocmi" {
				
			if "`X'"!="" {
				tempvar zbar
				egen `group'=group(`X')
				bys `group': gen `N_g'=_N
				bys `group': egen double `zbar'=mean(`Z')
				count if `N_g'==1&`touse'
				replace `touse'=`touse'&`N_g'>1
				}
			else {
				tempvar constant
				gen `constant'=1
				su `Z' if `touse', meanonly
				loc zbar=r(mean)	
				}
			
			count if `touse'
			

			tempvar d
			gen `d'=.
			loc num=0
			loc N_ineq=0 
			tokenize `values'
			foreach j in `values' {
				if `j'==`J' continue
				loc ++num
				loc ++N_ineq
				tempvar moment`j' moment`j'
				loc moments `moments' `moment`j''
				replace `d'=`D'==`j' if `touse'
				
				if "`c'"=="c." {
					if "`X'"=="" {
						su `d', meanonly
						loc dbar=r(mean)
						}
					else {
						tempvar dbar
						bys `group': egen double `dbar'=mean(`d')
						}
					if `j'<`jstar' 	gen double `moment`j''=-(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'
					else gen double `moment`j''=(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'
					}
				else {
					if `j'<`jstar' 	gen double `moment`j''=`d'*(`zbar'-`Z')/(`zbar'*(1-`zbar')) if `touse'
					else gen double `moment`j''=`d'*(`Z'-`zbar')/(`zbar'*(1-`zbar')) if `touse'
					}
				}
			
			cmi_test (`moments') () `constant' `X' if `touse', `cmi_opts'
			}
		
		//Post results
		ereturn post `b' `V', depname("`D'>=j") esample(`touse') obs(`N')
		ereturn scalar F_stat4=`F4'
		ereturn scalar p_val4=`p4'
		ereturn scalar F_stat5=`F5'
		ereturn scalar p_val5=`p5'
		
		if "`cmi'"!="nocmi" {
			foreach stat in stat cv01 cv05 cv10 pval {
				loc cmi_`stat'=r(`stat')
				ereturn scalar cmi_`stat'=r(`stat')
				}
		
			ereturn scalar N_ineq=`N_ineq'
			}
		
		ereturn local title 	"Tests of instrument validity with multivalued binarized treatment"
		ereturn local cmdline 	"mvttest `0'"
		ereturn local cmd 		"mvttest"
		
		noi {
			di _newline
			di "`e(title)'"
			di "{hline 78}"
			di as text "Levels of `D': `values'"
			di as text "Threshold value: `jstar'"
			di _newline
			di "First stage effect of `Z' for various thresholds"
			eret di, allbaselevels
			if "`cmi'"!="nocmi" {
				if "`ftest'"!="" {
					if "`test2desc'"!="" {
						di "Test of Assumption 3: {col 45}F-value {col 67}`: di %12.4f `F2''"
						di "{col 5}`test3desc' {col 45}p-value {col 67}`: di %12.4f `p2''"
						di "{hline 78}"
						}
					}
				di "CMI-test of Assumption 3: {col 45}`stattype' test statistic {col 67}`: di %12.4f `cmi_stat''"
				di "{col 45}Critical value, 1% {col 67}`: di %12.4f `cmi_cv01''"
				di "{col 5}b_{j+1}>=b_j> for j<`thresholdj'{col 45}Critical value, 5% {col 67}`: di %12.4f `cmi_cv05''"
				di "{col 5}b_{j+1}<=b_j for j>=`thresholdj' {col 45}Critical value, 10% {col 67}`: di %12.4f `cmi_cv10''"
				di "{col 45}p-value {col 67}`: di %12.4f `cmi_pval''"
				di "{hline 78}"
				}
			di "Test of Assumption 4: {col 45}F-value {col 67}`: di %12.4f `F4''"
			di "{col 5}all b_j=0 except b_`thresholdj' {col 45}p-value {col 67}`: di %12.4f `p4''"
			di "{hline 78}"
			di "Test of Assumption 5: {col 45}F-value {col 67}`: di %12.4f `F5''"
			di "{col 5}all b_j are the same {col 45}p-value {col 67}`: di %12.4f `p5''"
			di "{hline 78}"

			if "`X'"!="" di as text "Conditioning variables: `X'"
			}
		
	}
	
end
}
