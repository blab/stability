ddg_output <- read.delim("~/Desktop/fhcrc/stability/augur/data_sets/ddg_output_1_21_16.txt", header=FALSE)
data = data.frame(trunk=as.logical(ddg_output$V3), tip=as.logical(ddg_output$V4), 
                  date=as.Date(ddg_output$V5), ddg_outgroup_1HA0=ddg_output$V6, 
                  ddg_outgroup_2YP7=ddg_output$V7, mutations_from_1968=ddg_output$V12)

data$classification = "branch"
data$classification[data$trunk==TRUE] = "trunk"
data$classification[data$tip==TRUE] = "twig"
data$classification = ordered(data$classification)
data$year=as.numeric(format(data$date, "%Y"))


quadratic_fit = lm(year~poly(mutations_from_1968,2), data=data)
data$predict_year = predict(quadratic_fit, data.frame(data))

p_predict = ggplot(data=data, aes(x =mutations_from_1968, year))
p_predict + geom_point() + stat_smooth(method="lm", formula = y~poly(x,2)) + labs(title="Quadratic fit to get estimated year from the number of mutations from 1968 sequence", x="Mutations from 1968 structure sequence", y="Estimated Year")

p_pr = ggplot(data=data, aes(x =year,predict_year))
p_pr + geom_point() + stat_smooth(method="lm") + labs(title="Result of Quadratic Fit to get Estimated Years", x="Actual Year", y="Estimated Year")

# Residuals for each structure from Loess fit
loess_fit_1HA0 = loess(ddg_outgroup_1HA0~predict_year, data=data)
data$predict_ddg_1HA0 = predict(loess_fit_1HA0)
data$resid_1HA0 = data$ddg_outgroup_1HA0 - data$predict_ddg_1HA0
p_resid_1HA0 = ggplot(data=data, aes(x=predict_year,resid_1HA0))
p_resid_1HA0 + geom_point() + geom_hline(yintercept = 0, colour = 'red') + labs(title="Residual Plot for 1968 structure", x="Predicted Year", y="Residual")

loess_fit_2YP7 = loess(ddg_outgroup_2YP7~predict_year, data=data)
data$predict_ddg_2YP7 = predict(loess_fit_2YP7)
data$resid_2YP7 = data$ddg_outgroup_2YP7 - data$predict_ddg_2YP7
p_resid_2YP7 = ggplot(data=data, aes(x=predict_year,resid_2YP7))
p_resid_2YP7 + geom_point() + geom_hline(yintercept = 0, colour = 'red') + labs(title="Residual Plot for 2005 structure", x="Predicted Year", y="Residual")

p_resid = ggplot(data=data, aes(x=resid_1HA0, resid_2YP7, colour = predict_year))
p_resid + geom_point() + stat_smooth(method='lm',formula=y~x) + labs(title="Detrended ddG for 1968 and 2005 structures", x="1968 Residuals", y="2005 Residuals")

p_resid = ggplot(data=data, aes(x=ddg_outgroup_1HA0, ddg_outgroup_2YP7, colour = predict_year))
p_resid + geom_point() + stat_smooth(method='lm',formula=y~x) + labs(title="Detrended ddG for 1968 and 2005 structures", x="1968 Residuals", y="2005 Residuals")


# Loess curves drawn for each classification
p_classify_ddg_1HA0 = ggplot(data=data, aes(x=predict_year, ddg_outgroup_1HA0, colour = classification))
p_classify_ddg_1HA0 + geom_point() + geom_smooth(method="loess", fill=NA) + labs(title="Change in Stability over time for 1968 structure for trunk, branch and twig", x="Predicted Year", y="Change in Stability from Outgroup on 1HA0 (ddg)")
p_classify_ddg_2YP7 = ggplot(data=data, aes(x=predict_year, ddg_outgroup_2YP7, colour = classification))
p_classify_ddg_2YP7 + geom_point() + geom_smooth(method="loess", fill=NA) + labs(title="Change in Stability over time for 2005 structure for trunk, branch and twig", x="Predicted Year", y="Change in Stability from Outgroup on 2YP7 (ddg)")

# Jitter plots for each classification showing ddg
p_jitter_1HA0 = ggplot(data=data, aes(x =classification, ddg_outgroup_1HA0, colour=predict_year))
p_jitter_1HA0 + geom_jitter() + labs(title="Change in Stability over time for 1968 structure for trunk, branch and twig colored by year", y="Change in Stability from Outgroup on 1HA0 (ddg)")
p_jitter_2YP7 = ggplot(data=data, aes(x =classification, ddg_outgroup_2YP7, colour=predict_year))
p_jitter_2YP7 + geom_jitter() + labs(title="Change in Stability over time for 2005 structure for trunk, branch and twig colored by year", y="Change in Stability from Outgroup on 2YP7 (ddg)")

