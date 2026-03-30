* Create combined clock variables (using whichever array the participant has)
* Age 9
gen z_grimAge_9 = z_grimAge_me9
replace z_grimAge_9 = z_grimAge_mk9 if missing(z_grimAge_9)

gen z_phenoAge_9 = z_phenoAge_me9
replace z_phenoAge_9 = z_phenoAge_mk9 if missing(z_phenoAge_9)

gen z_poam45_9 = z_poam45_me9
replace z_poam45_9 = z_poam45_mk9 if missing(z_poam45_9)

* Age 15
gen z_grimAge_15 = z_grimAge_me15
replace z_grimAge_15 = z_grimAge_mk15 if missing(z_grimAge_15)

gen z_phenoAge_15 = z_phenoAge_me15
replace z_phenoAge_15 = z_phenoAge_mk15 if missing(z_phenoAge_15)

gen z_poam45_15 = z_poam45_me15
replace z_poam45_15 = z_poam45_mk15 if missing(z_poam45_15)

* Create array indicator
gen array_type_9 = .
replace array_type_9 = 1 if !missing(z_grimAge_me9)  // EPIC
replace array_type_9 = 0 if !missing(z_grimAge_mk9)  // 450K

gen array_type_15 = .
replace array_type_15 = 1 if !missing(z_grimAge_me15)  // EPIC
replace array_type_15 = 0 if !missing(z_grimAge_mk15)  // 450K

label define array_lbl 0 "450K" 1 "EPIC"
label values array_type_9 array_type_15 array_lbl

* Create combined cell composition variables
* Age 9
gen cell_epi_9 = k5me_epi
replace cell_epi_9 = k5mk_epi if missing(cell_epi_9)

gen cell_fib_9 = k5me_fib
replace cell_fib_9 = k5mk_fib if missing(cell_fib_9)

gen cell_ic_9 = k5me_ic
replace cell_ic_9 = k5mk_ic if missing(cell_ic_9)

* Age 15
gen cell_epi_15 = k6me_epi
replace cell_epi_15 = k6mk_epi if missing(cell_epi_15)

gen cell_fib_15 = k6me_fib
replace cell_fib_15 = k6mk_fib if missing(cell_fib_15)

gen cell_ic_15 = k6me_ic
replace cell_ic_15 = k6mk_ic if missing(cell_ic_15)

* Verify counts
tab array_type_9, m
tab array_type_15, m


* Example: Main model with array type as covariate
regress z_grimAge_9 ln_incbirth_imp female i.race age_9 i.array_type_9, vce(cluster idnum)

* Sensitivity: Add cell type proportions
regress z_grimAge_9 ln_incbirth_imp female i.race age_9 i.array_type_9 cell_epi_9 cell_ic_9, vce(cluster idnum)

