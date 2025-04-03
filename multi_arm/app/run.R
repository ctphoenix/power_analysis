# run.R for multi_arm/app

# Load the rsconnect library
library(rsconnect)

# Deploy the application in the current directory ('multi_arm/app')
# Assumes you have already configured your rsconnect account
# using setAccountInfo() or connectUser()
rsconnect::deployApp(
  appDir = ".", # Deploy the current directory
  appName = "multi_arm_power_sim", # Optional: specify an app name
  forceUpdate = TRUE
)

# Optional: List applications after deployment
# print(rsconnect::applications()) 