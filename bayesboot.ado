program bayesboot, eclass
	syntax anything, [cluster(varlist) absorb(string) opts(string)]
	
	tempname r
	gen `r'=rgamma(2,2)
	if "`cluster'"!="" bys `cluster': replace `r'=`r'[1]
	reghdfe `anything' [aw=`r'], absorb(`absorb') `opts'
	drop `r'
end
