*! mvttest v 1.1 09aug2022
*! Author: Martin E. Andresen
*! For "Instrument-based estimation with binarized treatments: Issues and tests for the exclusion restriction", joint with Martin Huber

cap program drop mvttest bayesboot
{
	program define mvttest, eclass
		version 13.0
		syntax varlist(min=2) [if] [in], jstar(numlist max=1 min=1) [ /*
				*/cmi_opts(string) 		/* Any other options for the cmi_test
				*/graph_opts(string)  	/* Options for the twoway graph. Overrides default
				*/noPlot				/* Do not plot graph 
				*/keepsingletons		/* Keeps observations in singleton groups
				*/nocmi					/* Do not perform CMI-tests
				*/cluster(varlist)		/* Cluster standard errors
				*/bootreps(integer 0)	/* Perform sup-t-tests in addition/instead of CMI-tests
				*/]
	
	
	qui {
		marksample touse
		preserve
		
		keep `touse' `varlist' `cluster'
		keep if `touse'
		
		//Parse and check input
		tokenize `varlist'
		loc D "`1'"
		mac shift
		loc Z "`1'"
		mac shift
		loc X `*'
		
		if strpos("`cmi_opts'","ks")==0 loc stattype "CvM"
		else loc stattype KS
		
		_fv_check_depvar `Z'
		tab `Z' if `touse'
		if `r(r)'==2 loc c 1. //binary instrument
		else loc c c. //continuous instrument
		
		//chech user written programs installed
		foreach prog in reghdfe cmi_test {
			capture which `prog'
			if _rc!=0 {
				noi di in red "mvttest uses `prog', which is not installed. Install using ssc install `prog'"
				exit 
				}
			}
		
		if `bootreps'<0 {
			noi di in red "Bootstrap reps specified in bootreps() must be nonnegative."
			exit
		}
		
		tempvar id regno group d expand include
		tempname b0 boot orig
		
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
	
		tempvar d
		gen `d'=`D'!=`regno'-1 if `regno'-1<`jstar'
		replace `d'=`D'==`regno'-1 if `regno'-1>=`jstar'
		replace `D'=`D'>=`regno'
		
		if "`cluster'"=="" loc cluster `id'
		reghdfe `D' `c'`Z'`xabs'#`regno' if `include', absorb(`regno'`xabs') vce(cluster `cluster') `keepsingletons' nocons
		est sto `orig'
		
		//perform F-tests of A3* and A4
		foreach Xgroup in `Xlevs' {
			loc test3star
			if "`X'"=="" loc loc pre
			else loc pre `Xgroup'#
			loc i=0
			foreach j in `values' {
				if `j'==`jzero' continue
				capture di _b[`c'`Z'#`pre'`j'.`regno']
				if _rc==0 {
					loc ++i
					if `j'<`jstar' {
						loc test3star `test3star' `c'`Z'#`pre'`j'.`regno'
						if `i'==1 loc test4 `c'`Z'#`pre'`j'.`regno'	
						else loc test4 `test4'=`c'`Z'#`pre'`j'.`regno'	
						}
					else if `j'>`thresholdj' {
						loc test3star `test3star' `c'`Z'#`pre'`j'.`regno'
						if `i'==1 loc test4  `c'`Z'#`pre'`j'.`regno'
						else loc test4 `test4'=`c'`Z'#`pre'`j'.`regno'	
						}					
					}
				}
				loc teststring3star `teststring3star' (`test3star')
				capture di _b[`c'`Z'#`pre'`thresholdj'.`regno']
				if _rc==0 loc teststring4 `teststring4' (`test4' = `c'`Z'#`pre'`thresholdj'.`regno' )
				else loc teststring4 `teststring4' (`test4')
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
		else {
			local colnames `=subinstr("`colnames'","#co.`Z'","",.)'
			local colnames `=subinstr("`colnames'","#c.`Z'","",.)'
		}
		local colnames `=subinstr("`colnames'","`regno'","j",.)'
		
		mat colnames `b'=`colnames'
		mat colnames `V'=`colnames'
		mat rownames `V'=`colnames'
		
		//plot figure
		tempname r
		if "`plot'"!="noplot" {
			if `J'>20 loc small ,labsize(vsmall) angle(90)
			if "`X'"!="" {
				tempname noX maxvios

				reghdfe `D' `c'`Z'#`regno', absorb(`regno') nosample `keepsingletons' nocons vce(cluster `cluster')
				drop _all
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
				if "`c'"=="c." {
					split parm`=`nvars'-1', parse(.)
					destring parm`=`nvars'-1'1, ignore(o b) replace
					rename parm`=`nvars'-1'1 j
					drop parm parm`=`nvars'-1'2 parm`=`nvars'-1'
				}
				
				else {
					split parm`nvars', parse(.)
					destring parm`nvars'1, ignore(o b) replace
					rename parm`nvars'1 j
					drop  parm parm`nvars'2 parm`nvars'
				}
				
				egen group=group(parm*)
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
			
			drop _all
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
			if "`c'"=="1." {
				split parm2, parse(.)
				destring parm21, replace ignore(o b)
				rename parm21 j
			}
			else {
				split parm1, parse(.)
				destring parm11, replace ignore(o b)
				rename parm11 j
			}
			
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
					xtitle("`D' at least") legend(off) title("First stage effect of `Z' using various thresholds") ///
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
		
	
		//Sup-t test
		if `bootreps'>0 {
			tempname b0 t0 Vboot se0 sort pmat
			reghdfe `d' `c'`Z'`xabs'#`regno' if `include', absorb(`regno'`xabs') vce(cluster `cluster') `keepsingletons' nocons
			mat `b0'=e(b)
			loc names: colnames `b0'
			
			mat `Vboot'=e(V)
			mata: st_matrix("`se0'",sqrt(diagonal(st_matrix("`Vboot'")))')
			mata: st_matrix("`t0'",st_matrix("`b0'"):/st_matrix("`se0'"))
			mata: st_matrix("`sort'",sort(st_matrix("`t0'")',1)')
			
			noi bayesboot `d' `c'`Z'`xabs'#`regno', absorb(`regno'`xabs') cluster(`cluster') opts(`keepsingletons' nocons)
			noi simulate  _b _se, reps(`bootreps'): bayesboot `d' `c'`Z'`xabs'#`regno', absorb(`regno'`xabs') cluster(`cluster') opts(`keepsingletons' nocons)
			
			loc numcells=colsof(`b0')
			forvalues i=1/`numcells' {
				gen t`i'=(_sim_`i'-`b0'[1,`i'])/_sim_`=`i'+`numcells''
			}

			loc i=0
			foreach name in `names' {
				loc ++i
				if abs(`t0'[1,`i']-`sort'[1,1])<1e-12 continue, break
			}

			egen mint=rowmin(t*)
			loc tval_boot=`t0'[1,`i']
			count if mint<.
			loc validreps=r(N)
			count if mint<=`t0'[1,`i']
			loc pboot =  `=(r(N)+1)/(`validreps'+1)'
			noi su mint, d
			noi _pctile mint, percentiles(1 5 10)
			loc boot_cv01=r(r1)
			loc boot_cv05=r(r2)
			loc boot_cv10=r(r3)
		}
		
		
		//CMI-test
		
		if "`cmi'"!="nocmi" {
					
			//determine coef comparisons  - drop j if beta_j+1 and beta_j cannot both be estimated (within a cell of  X)
			drop _all
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
			drop if se==0|se==.
			split parm, parse(#)
			loc nvars=r(nvars)
			
			loc nvars=r(nvars)
							if "`c'"=="c." {
					split parm`=`nvars'-1', parse(.)
					destring parm`=`nvars'-1'1, ignore(o b) replace
					rename parm`=`nvars'-1'1 j
					drop parm parm`=`nvars'-1'2 parm`=`nvars'-1'
				}
				
				else {
					split parm`nvars', parse(.)
					destring parm`nvars'1, ignore(o b) replace
					rename parm`nvars'1 j
					drop  parm parm`nvars'2 parm`nvars'
				}
				
			
			egen group=group(parm*)
			sort group j
			bys group: drop if _n==_N
			levelsof j, local(levj)
			count
			loc ineq_cells=r(N)

			
			//construct moment inqualities	
			drop _all
			u `tmpdata'
			if "`X'"!="" {
				tempvar zbar
				bys `group': egen double `zbar'=mean(`Z') if `touse'
				if "`c'"=="c." tempvar dbar
				}
			else {
				su `Z' if `touse', meanonly
				loc zbar=r(mean)	
				}
			
			tempvar d denom
			gen `d'=.
			
			
			loc N_ineq=0
			if "`c'"=="c." {
				gen double `denom'=(`Z'-`zbar')^2
				noi su `denom'
				loc denom=r(sum)
				loc N=r(N)
			}
			foreach j in `levj' {
				if inlist(`j',`J') continue
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
					
					if `j'<`jstar' 	gen double `moment`N_ineq''=-`N'*((`d'-`dbar')*(`Z'-`zbar'))/`denom' if `touse'
					else gen double `moment`N_ineq''=`N'*((`d'-`dbar')*(`Z'-`zbar'))/`denom' if `touse'
					}
				else {
					if `j'<`jstar' 	gen double `moment`N_ineq''=`d'*(`zbar'-`Z')/(`zbar'*(1-`zbar')) if `touse'
					else gen double `moment`N_ineq''=`d'*(`Z'-`zbar')/(`zbar'*(1-`zbar')) if `touse'
					}
				
				loc moments `moments' `moment`N_ineq''
			
				}
			cmi_test (`moments') () `X' if `touse', `cmi_opts'
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
				di "{col 5}b_{j+1}>=b_j for `thresholdj'>j>=`jzero' {col 45}Critical value, 5% {col 67}`: di %12.4f `cmi_cv05''"
				di "{col 5}b_j>=b_{j+1} for `J'>j>=`thresholdj' {col 45}Critical value, 10% {col 67}`: di %12.4f `cmi_cv10''"
				di "{col 45}p-value {col 67}`: di %12.4f `cmi_pval''"
				di "{hline 78}"
				}
			if `bootreps'>0 {
				di "Bootstrap test of Assumption 3: {col 45}Test statistic {col 67}`: di %12.4f `tval_boot''"
				di "{col 45}Critical value, 1% {col 67}`: di %12.4f `boot_cv01''"
				di "{col 5}b_{j+1}>=b_j for `thresholdj'>j>=`jzero' {col 45}Critical value, 5% {col 67}`: di %12.4f `boot_cv05''"
				di "{col 5}b_j>=b_{j+1} for `J'>j>=`thresholdj' {col 45}Critical value, 10% {col 67}`: di %12.4f `boot_cv10''"
				di "{col 45}p-value {col 67}`: di %12.4f `pboot''"
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