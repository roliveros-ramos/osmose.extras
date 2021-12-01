
configDir4  = "../osmose-peru/osmose-hum_v4"
configFile4a = file.path(configDir4, "osmose4.3.3_hum.R")

xx = .readConfiguration(file=configFile4a)


inifile = "maite.R"
initialize_osmose(input=configFile4a, output=inifile)


sim = init$init$anchoveta
par(mfrow=c(2,3), mar=c(3,3,1,1), oma=c(1,3,3,1))
plot(sim$biomass, type="l", ylim=c(0, 1.2*max(sim$biomass, sim$observed$biomass, na.rm=TRUE)), las=1)
points(sim$observed$biomass, pch=19, col="blue")
plot(sim$yield, type="h", ylim=c(0, 1.2*max(sim$yield)), las=1)
points(sim$observed$yield, pch=19, col="blue")
plot(sim$F, type="h", las=1)
plot(sim$size, sim$distB, type="h", las=1)
plot(sim$size, sim$selectivity, type="l")
title(main=.getPar(this, "species.name"), outer=TRUE)

