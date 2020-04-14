*! mvttest v 0.9 14apr2020
*! Author: Martin E. Andresen
*! For "Testing Binarized Treatments: Issues and tests for the exclusion restriction, joint with Martin Huber

cap program drop mvttest
{
	program define mvttest, eclass
		version 13.0
		syntax varlist(min=2) [fweight pweight] [if] [in], jstar(numlist max=1 min=1) [ /*
				*/cmi_opts(string) 		/* Any other options for the cmi_test
				*/graph_opts(string)  	/* Options for the twoway graph. Overrides default
				*/noPlot				/* Do not plot graph 
				*/keepsingletons		/* Keeps observations in singleton groups
				*/nocmi					/* Do not perform CMI-tests
				*/constrained			/* Constrained coefficients on Z to be the same across cells of X
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

		tempvar id regno group
		gen `id'=_n
		
		keep if `touse'
		sort `D'
		
		if "`X'"==""|"`constrained'"!="" loc plottype normal
		else loc plottype minmax
		
		if "`X'"!="" {
			tempvar N_g 
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
			local numX: word count `X'
			if `numX'==1 loc Xspand i.`X'
			loc no=0
			else foreach var in `X' {
				loc ++no
				if `no'==1 loc Xspand `var'
				else loc Xspand `Xspand'#`var'
				}
			fvexpand `Xspand' if `touse'
			loc Xlevs `r(varlist)'
			if "`constrained'"!="" loc Xlevs 1 
			}
		else {
			gen `group'=1
			loc Xlevs 1
			}	
		levelsof `D' if `touse', local(values)
		local numvals: word count `values'
		if `numvals'<=2 {
			noi di in red "Treatment `D' not multivalued"
			exit
			}
		
		
		loc no=0
		foreach j in `values' {
			loc ++no
			if `no'>1 gen `D'`j'=`D'>=`j'
			}
			
		su `D' if `touse'
		loc J=r(max)
		loc jzero=r(min)
		
		loc no=0
		foreach j in `values' {
			if `j'<`jstar' loc valuesbelow `valuesbelow' `j'
			if `j'>=`jstar' {
				loc ++no
				if `no'==1 loc thresholdj=`j'
				else loc valuesabove `valuesabove' `j'
				}
			}
		
		foreach var in `X' {
			loc xabs `xabs'#`var'
			}

		
		drop `D'
		
		reshape long `D' `d', i(`id') j(`regno')
		
		if "`constrained'"=="" reghdfe `D' `c'`Z'`xabs'#`regno' if `touse', absorb(`regno'`xabs') vce(cluster `id') `keepsingletons' nocons
		else reghdfe `D' `c'`Z'#`regno' if `touse', absorb(`regno'`xabs') vce(cluster `id') `keepsingletons' nocons
		
		loc N=e(N)
		
		tempname b V
		mat `b'=e(b)
		mat `V'=e(V)
		
		local colnames: colnames `b'
		local colnames `=subinstr("`colnames'","`c'`Z'#","",.)'
		if "`c'"=="1." local colnames `=subinstr("`colnames'","1o.`Z'#","",.)'
		local colnames `=subinstr("`colnames'","`regno'","j",.)'
		
		mat colnames `b'=`colnames'
		mat colnames `V'=`colnames'
		mat rownames `V'=`colnames'
		

		if "`plot'"!="noplot" {
			if `J'>20 loc small ,labsize(vsmall) angle(90)
			if "`plottype'"=="minmax" {
				loc i=0
				foreach val in `values' {
					loc ++i
					if `i'==2 loc prevval=`val'
					if `i'<=2 continue
					loc maxvio
					loc maxparm
					foreach xgroup in `Xlevs' {
						if `prevval'<`jstar' loc testval=_b[`c'`Z'#`xgroup'#`prevval'.`regno']-_b[`c'`Z'#`xgroup'#`val'.`regno']
						else loc testval=_b[`c'`Z'#`xgroup'#`val'.`regno']-_b[`c'`Z'#`xgroup'#`prevval'.`regno']
						if "`maxvio'"=="" loc maxvio=`testval'
						else if `testval'>`maxvio' loc maxvio=`testval'
						if `maxvio'==`testval' loc maxparm `testval'
						}
					loc maxparms `maxparms' `maxparm'
					loc prevval=`val'
					}
				
				reghdfe `D' `c'`Z'#`regno' if `touse', absorb(`regno'`xabs') vce(cluster `id') `keepsingletons' nocons	
				}

				
			parmest, norestore
			loc no=0
			tempvar jno
			gen `jno'=_n
			foreach j in `values' {
				if (`j'==`jzero') continue
				loc ++no
				if "`plottype'"=="normal"|`no'>1 loc xlabel `xlabel' `no' "`j'"
				if `j'>=`jstar'&"`xline'"=="" {
					if "`plottype'"=="normal" loc xline=`=`no'-0.5'
					else loc xline=`=`no'+0.5'
					}
				}
			
			if "`plottype'"=="minmax" {
				gen vio=estimate[_n-1]-estimate[_n] if _n<`xline'
				replace vio=estimate[_n]-estimate[_n-1] if _n>=`xline'
				drop if vio==.
				gen double maxvio=.
				
				loc i=0
				foreach val in `maxparms' {
					loc ++i
					replace maxvio=`val' in `i'
					}
				gen jnomax=`jno'+0.2
				replace `jno'=`jno'-0.2
				
				
				twoway (bar vio `jno', color(navy) lcolor(white) lwidth(medium) barwidth(0.4)) (bar maxvio jnomax, color(maroon) lcolor(white) lwidth(medium) barwidth(0.4)) ///
					,  scheme(s1color) graphregion(color(white)) plotregion(lcolor(black)) ///
					xtitle("cutoff value j") title("Maximum violation across cells of X") ///
					legend(label(1 "violation without X") label(2 "maximum violation across cells of X") ring()) ///
					xlabel(`xlabel' `small') xline(`xline', lcolor(black) lpattern(dash)) `graph_opts'
				}
			else {
				twoway (bar estimate `jno', color(navy) lcolor(white) lwidth(medium)) (rcap min95 max95 `jno', color(navy)) ///
					,  scheme(s1color) graphregion(color(white)) plotregion(lcolor(black)) ///
					xtitle("`D' at least") legend(off) title("First stage effect using various thresholds") ///
					xlabel(`xlabel' `small') xline(`xline', lcolor(black) lpattern(dash)) `graph_opts'
				}
				
			}
			
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
				gen `group'=1
				su `Z' if `touse', meanonly
				loc zbar=r(mean)	
				}
			
			tempvar d
			gen `d'=.
			loc N_ineq=0 
			levelsof `group' if `touse', local(grouplevs)
			if "`constrained'"!="" loc grouplevs 1
			foreach j in `values' {
				if `j'==`J'|`j'==`jzero' continue
				replace `d'=`D'==`j' if `touse'
				if "`X'"=="" {
					su `d', meanonly
					loc dbar=r(mean)
						}
					else {
						tempvar dbar
						bys `group': egen double `dbar'=mean(`d')
						}
				foreach groupval in `grouplevs' { 
					loc ++N_ineq
					tempvar moment`N_ineq' 
					loc moments `moments' `moment`N_ineq''
					if "`constrained'"=="" {
						if "`c'"=="c." {					
							if `j'<`jstar' 	gen double `moment`N_ineq''=-(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'&`group'==`groupval'
							else gen double `moment`N_ineq''=(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'&`group'==`groupval'
							}
						else {
							if `j'<`jstar' 	gen double `moment`N_ineq''=`d'*(`zbar'-`Z')/(`zbar'*(1-`zbar')) if `touse'&`group'==`groupval'
							else gen double `moment`N_ineq''=`d'*(`Z'-`zbar')/(`zbar'*(1-`zbar')) if `touse'&`group'==`groupval'
							}
						}
					else {
						if "`c'"=="c." {					
							if `j'<`jstar' 	gen double `moment`N_ineq''=-(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'
							else gen double `moment`N_ineq''=(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'
							}
						else {
							if `j'<`jstar' 	gen double `moment`N_ineq''=`d'*(`zbar'-`Z')/(`zbar'*(1-`zbar')) if `touse'
							else gen double `moment`N_ineq''=`d'*(`Z'-`zbar')/(`zbar'*(1-`zbar')) if `touse'
							}
						}
					replace `moment`N_ineq''=0 if `moment`N_ineq''==.
					}
				
				}
			
			cmi_test (`moments') () `group' if `touse', `cmi_opts'
			
			if "`cmi'"!="nocmi" {
				foreach stat in stat cv01 cv05 cv10 pval {
					loc cmi_`stat'=r(`stat')
					}
				}
				
			}
		
		//Post results
		count if `touse'
		ereturn post `b' `V', depname("`D'>=j") esample(`touse') obs(`r(N)')
		
		//perform Chi2-tests of A4 and A5
		foreach Xgroup in `Xlevs' {
			loc test4
			loc test5
			if "`X'"==""|"`constrained'"!="" loc loc pre
			else loc pre `Xgroup'#
			foreach j in `values' {
				if `j'==`jzero' continue
				capture di _b[`pre'`j'.j]
				if _rc==0 {
					if `j'<`jstar' {
						loc test4 `test4' `pre'`j'.j
						loc test5 `test5' `pre'`j'.j =				
						}

					else if `j'>`thresholdj' {
						loc test4 `test4' `pre'`j'.j
						loc test5 `test5' `pre'`j'.j =
						}					
					}
				}
			loc teststring4 `teststring4' (`test4')
			loc teststring5 `teststring5' (`test5' `pre'`thresholdj'.j )
			}
			

		if "`teststring4'"!="" {
			test `teststring4'
			loc chi2_4=r(chi2)
			loc p4=r(p)
			}
		
		if "`teststring5'"!="" {
			test `teststring5'
			loc chi2_5=r(chi2)
			loc p5=r(p)
			}
		
		ereturn scalar chi2_4=`chi2_4'
		ereturn scalar p_val4=`p4'
		ereturn scalar chi2_5=`chi2_5'
		ereturn scalar p_val5=`p5'
		
		if "`cmi'"!="nocmi" {
			foreach stat in stat cv01 cv05 cv10 pval {
				ereturn scalar cmi_`stat'=`cmi_`stat''
				}
		
			ereturn scalar N_ineq=`N_ineq'
			}
		
		ereturn local title 	"Tests of instrument validity with multivalued binarized treatment"
		ereturn local cmdline 	"mvttest `0'"
		ereturn local cmd 		"mvttest"
		
		//Display results
		noi {
			di _newline
			di "`e(title)'"
			di "{hline 78}"
			di as text "Levels of `D': `valuesbelow' || {ul:{bf:`jstar'}} `valuesabove'"
			if "`X'"!="" di as text "Conditioning variables: `X'"
			if "`cmi'"!="nocmi" di as text "Number of inequalities tested: `=e(N_ineq)'"
			
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
			di "Test of Assumption 4: {col 45}Chi2-value {col 67}`: di %12.4f `chi2_4''"
			di "{col 5}all b_j=0 except b_`thresholdj' {col 45}p-value {col 67}`: di %12.4f `p4''"
			di "{hline 78}"
			di "Test of Assumption 5: {col 45}Chi2-value {col 67}`: di %12.4f `chi2_5''"
			di "{col 5}all b_j are the same {col 45}p-value {col 67}`: di %12.4f `p5''"
			di "{hline 78}"

			}
		
	}
	
end
}
