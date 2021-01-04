*! mvttest v 1.01 04jan2021
*! Author: Martin E. Andresen
*! For "Instrument-based estimation with binarized treatments: Issues and tests for the exclusion restriction", joint with Martin Huber

cap program drop mvttest
{
	program define mvttest, eclass
		version 13.0
		syntax varlist(min=2) [if] [in], jstar(numlist max=1 min=1) [ /*
				*/cmi_opts(string) 		/* Any other options for the cmi_test
				*/graph_opts(string)  	/* Options for the twoway graph. Overrides default
				*/noPlot				/* Do not plot graph 
				*/keepsingletons		/* Keeps observations in singleton groups
				*/nocmi					/* Do not perform CMI-tests
				*/]
	
	
	qui {
		marksample touse
		preserve
		
		keep `touse' `varlist'
		keep if `touse'
		
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
		
		//chech user written programs installed
		foreach prog in reghdfe cmi_test {
			capture which `prog'
			if _rc!=0 {
				noi di in red "mvttest uses `prog', which is not installed. Install using ssc install `prog'"
				exit 
				}
			}
		
		tempvar id regno group d expand include
				
		//drop singleton groups, determine cells of X
		if "`X'"!="" {
			tempvar N_g
			egen `group'=group(`X')
			if "`keepsingletons'"!="" {
				loc totcount=0
				while `count'!=0 {
					bys `group': gen `N_g'=_N 
					count if `N_g'==1
					drop if `N_g'==1
					loc count=r(N)
					loc totcount=`totcount'+`count'
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
			loc levelsX: word count `r(varlist)'
			
			foreach var in `X' {
				loc xabs `xabs'#`var'
				}
			}
		else {
			gen `group'=1
			loc Xlevs 1
			}
		levelsof `group' if `touse', local(grouplevs)
		
		//determine levels of D
		levelsof `D', local(values)
		local numvals: word count `values'
		if `numvals'<=2 {
			noi di in red "Treatment `D' not multivalued"
			exit
			}
		
		su `D'
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
		
		if "`cmi'"!="nocmi" {
			tempfile tmpdata
			save `tmpdata', replace
			}
		
		//expand and estimate b_j
		gen `id'=_n
		expand `=`numvals'-1'
		loc no=0
		gen `regno'=.
		foreach j in `values' {
			if `j'==`jzero' continue
			loc ++no
			bys `id': replace `regno'=`j' if _n==`no'
			}
		
		gen byte `include'=0
		if "`X'"!="" foreach xgroup in `grouplevs' {
			levelsof `D' if `group'==`xgroup', local(xvals)
			loc no=0
			foreach j in `xvals' {
				loc ++no
				if `no'==1 continue
				replace `include'=1 if `group'==`xgroup'&`regno'==`j'
				}
			}
		else replace `include'=1
	
		
		replace `D'=`D'>=`regno'

		reghdfe `D' `c'`Z'`xabs'#`regno' if `include', absorb(`regno'`xabs') vce(cluster `id') `keepsingletons' nocons
				
		//perform F-tests of A3* and A4
		foreach Xgroup in `Xlevs' {
			loc test3star
			loc test4
			if "`X'"=="" loc loc pre
			else loc pre `Xgroup'#
			foreach j in `values' {
				if `j'==`jzero' continue
				capture di _b[`c'`Z'#`pre'`j'.`regno']
				if _rc==0 {
					if `j'<`jstar' {
						loc test3star `test3star' `c'`Z'#`pre'`j'.`regno'
						loc test4 `test4' `c'`Z'#`pre'`j'.`regno' =				
						}
					else if `j'>`thresholdj' {
						loc test3star `test3star' `c'`Z'#`pre'`j'.`regno'
						loc test4 `test4' `c'`Z'#`pre'`j'.`regno' =
						}					
					}
				}
			loc teststring3star `teststring3star' (`test3star')
			loc teststring4 `teststring4' (`test4' `c'`Z'#`pre'`thresholdj'.`regno' )
			}
			

		foreach t in 3star 4 {
			if "`teststring`t''"!="" {
				test `teststring`t''
				loc F_`t'=r(F)
				loc p`t'=r(p)
				loc df`t'=r(df)
				loc df_r`t'=r(df_r)
				}
			}
			
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
		
		//plot figure
		tempname r
		if "`plot'"!="noplot" {
			if `J'>20 loc small ,labsize(vsmall) angle(90)
			if "`X'"!="" {
				tempname orig noX maxvios
				est sto `orig'
				reghdfe `D' `c'`Z'#`regno', absorb(`regno') nosample `keepsingletons' nocons vce(cluster `id')
				clear
				est sto `noX'
				est restore `orig'
				eret di, level(95)
				mat `r'=r(table)'
				loc parmnames: rownames `r'
				svmat `r', names(col)
				loc i=0
				gen parm=""
				foreach name in `parmnames' {
					loc ++i
					replace parm="`name'" in `i'
					}
				split parm, parse(#)
				loc nvars=r(nvars)
				split parm`nvars', parse(.)
				destring parm`nvars'1, ignore(o b) replace
				rename parm`nvars'1 j
				forvalues n=2/`=`nvars'-1' {
					if `n'==2 loc gen parm`n'
					else loc gen `gen'+parm`n'
					}
				gen str group=`gen'
				sort group j
				by group: gen maxvio=b-b[_n+1] if j<`jstar'
				by group: replace maxvio=-(b-b[_n+1]) if j>=`jstar'
				drop if maxvio==.
				sort j maxvio
				bys j: drop if _n!=_N
				keep maxvio j
				save `maxvios', replace
				est restore `noX'
				}
			
			clear
			eret di, level(95)
			mat `r'=r(table)'
			svmat `r', names(col)
			loc parmnames: rownames `r'
			loc i=0
			gen parm=""
			foreach name in `parmnames' {
				loc ++i
				replace parm="`name'" in `i'
				}
			split parm, parse(#)
			split parm2, parse(.)
			destring parm21, replace ignore(o b)
			rename parm21 j
			
			if "`X'"=="" {
				levelsof j, local(levj)
				loc no=0
				foreach lev in `levj' {
					loc ++no
					loc xlabels `xlabels' `no' "`lev'"
					if `lev'>=`jstar'&"`xline'"=="" loc xline=`no'-0.5
					}
				gen jno=_n
				twoway (bar b jno, color(navy) lcolor(white) lwidth(medium)) (rcap ll ul jno, color(navy)) ///
					,  scheme(s1color) graphregion(color(white)) plotregion(lcolor(black)) ///
					xtitle("`D' at least") legend(off) title("First stage effect using various thresholds") ///
					xlabel(`xlabels' `small') xline(`xline', lcolor(black) lpattern(dash)) `graph_opts'
				

				}
				
			else {
				gen vio=b-b[_n+1] if j<`jstar'
				replace vio=-(b-b[_n+1]) if j>=`jstar'
				drop if vio==.
				merge 1:1 j using `maxvios', nogen keep(1 3)
				levelsof j, local(levj)
				loc no=0
				foreach lev in `levj' {
					loc ++no
					loc xlabels `xlabels' `no' "`lev'"
					if `lev'>=`jstar'&"`xline'"=="" loc xline=`no'-0.5
					}
				gen jno=_n-0.2
				gen jmax=jno+0.4

				twoway (bar vio jno, color(navy) lcolor(white) lwidth(medium) barwidth(0.4)) (bar maxvio jmax, color(maroon) lcolor(white) lwidth(medium) barwidth(0.4)) ///
					,  scheme(s1color) graphregion(color(white)) plotregion(lcolor(black)) ///
					xtitle("`D' at least") title("Maximum violation across cells of X") ///
					legend(label(1 "violation without X") label(2 "maximum violation across cells of X") ring()) ///
					xlabel(`xlabels' `small') xline(`xline', lcolor(black) lpattern(dash)) `graph_opts'
					
				est restore `orig'
				est drop `orig' `noX'
				}
			
		
				
			}
			
		
		//CMI-test
		
		if "`cmi'"!="nocmi" {
					
			//determine coef comparisons  - drop j if beta_j+1 and beta_j cannot both be estimated (within a cell of  X)
			clear
			eret di, level(95)
			mat `r'=r(table)'
			loc parmnames: rownames `r'
			svmat `r', names(col)
			loc i=0
			gen parm=""
			foreach name in `parmnames' {
				loc ++i
				replace parm="`name'" in `i'
				}
			drop if se==0|se==.
			split parm, parse(#)
			loc nvars=r(nvars)
			split parm`nvars', parse(.)
			destring parm`nvars'1, ignore(o b) replace
			rename parm`nvars'1 j
			if "`X'"!="" {
				forvalues n=2/`=`nvars'-1' {
					if `n'==2 loc gen parm`n'
					else loc gen `gen'+parm`n'
					}
				gen str group=`gen'
				}
			else  rename parm1 group
			sort group j
			bys group: drop if _n==_N
			levelsof j, local(levj)
			count
			loc ineq_cells=r(N)

			
			//construct moment inqualities	
			u `tmpdata', clear
			if "`X'"!="" {
				tempvar zbar
				bys `group': egen double `zbar'=mean(`Z') if `touse'
				if "`c'"=="c." tempvar dbar
				}
			else {
				su `Z' if `touse', meanonly
				loc zbar=r(mean)	
				}
			
			tempvar d
			gen `d'=.
			
			
			loc N_ineq=0 
			foreach j in `levj' {
				if inlist(`j',`jzero',`J') continue
				loc ++N_ineq
				replace `d'=`D'==`j' if `touse'		
				
				tempvar moment`N_ineq' 
				if "`c'"=="c." {	
					if "`X'"=="" {
						su `d' if `touse', meanonly
						loc dbar=r(mean)
						}
					else {
						cap drop `dbar'
						bys `group': egen `dbar'=mean(`d') if `touse'
						}
					if `j'<`jstar' 	gen double `moment`N_ineq''=-(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'
					else gen double `moment`N_ineq''=(`d'*(`Z'-`zbar')-`dbar'*(`Z'-`zbar'))/(`Z'^2-2*`Z'*`zbar'+`zbar'^2) if `touse'
					}
				else {
					if `j'<`jstar' 	gen double `moment`N_ineq''=`d'*(`zbar'-`Z')/(`zbar'*(1-`zbar')) if `touse'
					else gen double `moment`N_ineq''=`d'*(`Z'-`zbar')/(`zbar'*(1-`zbar')) if `touse'
					}
				
				loc moments `moments' `moment`N_ineq''
				
				}
			cmi_test (`moments') () `group' if `touse', `cmi_opts'
			foreach stat in stat cv01 cv05 cv10 pval {
				loc cmi_`stat'=r(`stat')
				}		
			}
		
		//Post results
		restore
		count if `touse'
		ereturn post `b' `V', depname("`D'>=j") esample(`touse') obs(`r(N)')
		
		foreach  t in 3star 4 {
			ereturn scalar F`t'=`F_`t''
			ereturn scalar p_val`t'=`p`t''
			ereturn scalar df`t'=`df`t''
			ereturn scalar df_r`t'=`df_r`t''
			}
		if "`levelsX'"!="" ereturn scalar levels_X=`levelsX'
		
		if "`cmi'"!="nocmi" {
			foreach stat in stat cv01 cv05 cv10 pval {
				ereturn scalar cmi_`stat'=`cmi_`stat''
				}
		
			ereturn scalar N_ineq=`N_ineq'
			ereturn scalar N_ineqcells=`ineq_cells'
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
			eret di, noemptycells
			if "`cmi'"!="nocmi" {
				di "CMI-test of Assumption 3: {col 45}`stattype' test statistic {col 67}`: di %12.4f `cmi_stat''"
				di "{col 45}Critical value, 1% {col 67}`: di %12.4f `cmi_cv01''"
				di "{col 5}b_{j+1}>=b_j for `thresholdj'>j>0 {col 45}Critical value, 5% {col 67}`: di %12.4f `cmi_cv05''"
				di "{col 5}b_j>=b_{j+1} for J>j>=`thresholdj' {col 45}Critical value, 10% {col 67}`: di %12.4f `cmi_cv10''"
				di "{col 45}p-value {col 67}`: di %12.4f `cmi_pval''"
				di "{hline 78}"
				}
			di "Test of Assumption 3*: {col 45}F(`df3star',`df_r3star') {col 67}`: di %12.4f `F_3star''"
			di "{col 5}all b_j=0 except b_`thresholdj' {col 45}p-value {col 67}`: di %12.4f `p3star''"
			di "{hline 78}"
			di "Test of Assumption 4: {col 45}F(`df4',`df_r4') {col 67}`: di %12.4f `F_4''"
			if "`X'"=="" di "{col 5}all b_j are the same {col 45}p-value {col 67}`: di %12.4f `p4''"
			else di "{col 5}all b_j are the same within cells of X {col 45}p-value {col 67}`: di %12.4f `p4''"
			di "{hline 78}"

			}
	
	
	}
	
	
end
}
