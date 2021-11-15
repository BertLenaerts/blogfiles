
clear all
cls clears screen

import delimited "df.csv", clear 
qui tab id, gen(k)


//// GROUP-FIXED EFFECTS

xtset id year

// NO CORRECTION
reg y u x k*
xtreg y u x, fe

// HC CORRECTION
reg y u x k*, robust

// GROUP CLUSTER
reg y u x k*, vce(cluster id)
xtreg y u x, fe cluster(id) dfadj

// TIME CLUSTER
reg y u x k*, vce(cluster year)
* no xtreg as panels are not nested within clusters

// HIGHER-LEVEL CLUSTERING
reg y u x k*, vce(cluster gid)
xtreg y u x, fe cluster(gid) dfadj
cgmreg y u x k*, cluster(gid) nocons


//// TIME-FIXED EFFECTS

xtset year id

// NO CORRECTION
reg y u x i.year
xtreg y u x, fe

// HC CORRECTION
reg y u x i.year, robust

// GROUP CLUSTER
reg y u x i.year, vce(cluster id)
* no xtreg as panels are not nested within clusters

// TIME CLUSTER
reg y u x i.year, vce(cluster year)
xtreg y u x, fe cluster(year) dfadj

// HIGHER-LEVEL CLUSTERING
reg y u x i.year, vce(cluster gid)
* no xtreg as panels are not nested within clusters


//// TWOWAY FIXED EFFECTS

xtset id year

// NO CORRECTION
reg y u x i.year k*
xtreg y u x i.year, fe

// HC CORRECTION
reg y u x i.year k*, robust

// GROUP CLUSTER
reg y u x k* i.year, vce(cluster id)
xtreg y u x i.year, fe cluster(id) dfadj

// TIME CLUSTER
reg y u x k* i.year, vce(cluster year)
* no xtreg as panels are not nested within clusters

// HIGHER-LEVEL CLUSTERING
reg y u x k* i.year, vce(cluster gid)
xtreg y u x i.year, fe cluster(gid) dfadj
