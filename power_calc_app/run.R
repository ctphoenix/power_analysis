library(rsconnect)
rsconnect::setAccountInfo(name='ctphoenix',
                          token='879507F9586444A0C72448FD664E216C',
                          secret='YuvnUIDoImWiSX+/lE2Gbls0FuGlA3FUu9npCIXR')
rsconnect::deployApp("/Users/patrick.staples/repos/power/power_calc_app")

apps <- rsconnect::applications()
print(apps)

# Terminate a specific app by name
# rsconnect::terminateApp("r_app")
