ddg_output <- read.delim("~/Desktop/fhcrc/stability/augur/data_sets/ddg_output_1_15_16.txt", header=FALSE)
trunk = ddg_output$V3
tip = ddg_output$V4
date = ddg_output$V5
date = as.Date(date)
ddg_outgroup_1HA0 = ddg_output$V6
ddg_outgroup_2YP7 = ddg_output$V7
ddg_outgroup_average = (ddg_outgroup_1HA0+ddg_outgroup_2YP7)/2
ddg_parent_1HA0 = ddg_output$V8
ddg_parent_2YP7 = ddg_output$V9
ddg_parent_1HA0 = as.numeric(levels(ddg_parent_1HA0))[ddg_parent_1HA0]
ddg_parent_2YP7 = as.numeric(levels(ddg_parent_2YP7))[ddg_parent_2YP7]
ddg_parent_average = (ddg_parent_1HA0+ddg_parent_2YP7)/2
old.par <- par(mfrow=c(1, 1))
plot(date, ddg_outgroup_1HA0, main="Change in stability from the outgroup over time for 1968 structure\n(decreasing stability as moving away from 1968)", yaxt = "n")
axis(2, at=c(-60,-40,-20,0,20,40,60,80,100,120,140,160))
plot(date, ddg_outgroup_2YP7, main="Change in stability from the outgroup over time for 2011 structure\n(increasing stability as getting closer to 2011)")
plot(date, ddg_outgroup_average, main="Change in stability from the outgroup over time for average of two structures")
plot(date, ddg_parent_1HA0, main="Change in stability from the parent node over time for 1968 structure")
plot(date, ddg_parent_2YP7, main="Change in stability from the parent node over time for 2011 structure")
plot(date, ddg_parent_average, main="Change in stability from the parent node over time for average of two structures")
plot(ddg_outgroup_1HA0,ddg_outgroup_2YP7, main="ddG from outgroup for 1968(1HA0) and 2011(2YP7) structure\n-0.629 correlation")
cor(ddg_outgroup_1HA0,ddg_outgroup_2YP7)
par(old.par)