
Comprehensive list of output variables. Variable types include instantaneous (I) values at the archival time step, means between archival times (M), and miNimum (N) or maXimum (X) values between archival times. The status column shows the availability of the variables via THREDDS: available at time of writing (T), available later (L), or not available (blank). 
When a variable is shared on THREDDS at a coarser frequency than what was originally saved, the THREDDS frequency is indicated in parenthesis.

|Variable    |Description                                                |Frequency  |Realm   |Type|Status|
|------------|-----------------------------------------------------------|-----------|--------|----|------|
|baresoilFrac|Bare Soil Fraction                                         |3-hr       |Land 2D |I   |      |
|cape        |Convective Available Potential Energy                      |3-hr       |Atmo 2D |I   |L     |
|clt         |Total Cloud Fraction                                       |1-hr       |Atmo 2D |M   |L     |
|clwvi       |Condensed Water Path                                       |3-hr       |Atmo 2D |I   |T     |
|ct          |Air Pressure at Cloud Top                                  |3-hr       |Atmo 2D |I   |      |
|dds         |Near-Surface Dewpoint Depression                           |3-hr       |Atmo 2D |I   |      |
|dedt        |Tendency of Integrated Cloud Water and Ice                 |monthly    |Atmo 2D |M   |      |
|dfq         |Horizontal Divergence of Water Vapor Flux                  |monthly    |Atmo 2D |M   |      |
|drdt        |Tendency of Integrated Water Vapor                         |monthly    |Atmo 2D |M   |      |
|ebq         |Residual Term of Atmosphere Water Budget                   |monthly    |Atmo 2D |M   |      |
|evspsbl     |Evaporation                                                |1-hr       |Atmo 2D |M   |T     |
|evspsblland |Water Evaporation from Land                                |3-hr       |Land 2D |M   |      |
|hfls        |Surface Upward Latent Heat Flux                            |3-hr       |Atmo 2D |M   |T     |
|hfss        |Surface Upward Sensible Heat Flux                          |3-hr       |Atmo 2D |M   |L     |
|hrmax       |Maximum Near-Surface Relative Humidity                     |3-hr       |Atmo 2D |X   |      |
|hrmin       |Minimum Near-Surface Relative Humidity                     |3-hr       |Atmo 2D |N   |      |
|hurs        |Near-Surface Relative Humidity                             |1-hr       |Atmo 2D |I   |T     |
|hus         |Specific Humidity                                          |3-hr (6-hr)|Atmo 3D |I   |L     |
|huss        |Near-Surface Specific Humidity                             |1-hr       |Atmo 2D |I   |T     |
|lflrt       |Lake Floor Temperature                                     |3-hr       |Lake 2D |I   |      |
|lif         |Lake Ice Fraction                                          |3-hr       |Lake 2D |I   |      |
|lit         |Lake Ice Thickness                                         |3-hr       |Lake 2D |I   |      |
|lmlt        |Lake Mixed-Layer Temperature                               |3-hr       |Lake 2D |I   |      |
|lmlthick    |Lake Mixed-Layer Thickness                                 |3-hr       |Lake 2D |I   |      |
|mrfsl       |Soil Layer Frozen Water Content                            |monthly    |Land 3D |I   |      |
|mrfso       |Soil Frozen Water Content                                  |3-hr       |Land 2D |I   |L     |
|mrlsl       |Water Content of Soil Layer                                |monthly    |Land 3D |I   |      |
|mrro        |Total Runoff                                               |3-hr       |Land 2D |M   |T     |
|mrros       |Surface Runoff                                             |3-hr       |Land 2D |M   |T     |
|mrso        |Total Soil Moisture Content                                |3-hr       |Land 2D |I   |      |
|mrsos       |Moisture in Upper Portion of Soil Column                   |3-hr       |Land 2D |I   |L     |
|pr          |Precipitation Flux                                         |1-hr, daily|Atmo 2D |M   |T     |
|prc         |Convective Precipitation                                   |1-hr       |Atmo 2D |M   |L     |
|prdc        |Deep Convective Precipitation                              |monthly    |Atmo 2D |M   |      |
|prfr        |Freezing Rain                                              |1-hr       |Atmo 2D |M   |T     |
|prlp/prra   |Liquid Precipitation/Rainfall Flux                         |3-hr       |Atmo 2D |M   |T     |
|prrp        |Refrozen Rain                                              |3-hr       |Atmo 2D |M   |      |
|prsn        |Snowfall Flux                                              |3-hr       |Atmo 2D |M   |T     |
|prw         |Water Vapor Path                                           |3-hr       |Atmo 2D |I   |T     |
|ps          |Surface Air Pressure                                       |1-hr       |Atmo 2D |I   |T     |
|psl         |Sea Level Pressure                                         |1-hr       |Atmo 2D |I   |T     |
|pw          |Precipitable Water                                         |3-hr       |Atmo 2D |I   |      |
|rhsmax      |Surface Daily Maximum Relative Humidity                    |daily      |Atmo 2D |X   |      |
|rhsmin      |Surface Daily Minimum Relative Humidity                    |daily      |Atmo 2D |N   |      |
|rlds        |Surface Downwelling Longwave Radiation                     |1-hr       |Atmo 2D |M   |T     |
|rls         |Net LW Surface Radiation                                   |1-hr       |Atmo 2D |M   |      |
|rlus        |Surface Upwelling Longwave Radiation                       |1-hr       |Atmo 2D |M   |L     |
|rlut        |TOA Outgoing Longwave Radiation                            |1-hr       |Atmo 2D |M   |L     |
|rsaa        |Shortwave Radiation Absorbed by Atmosphere                 |1-hr       |Atmo 2D |M   |      |
|rsds        |Surface Downwelling Shortwave Radiation                    |1-hr       |Atmo 2D |M   |T     |
|rsdt        |TOA Incident Shortwave Radiation                           |1-hr       |Atmo 2D |M   |L     |
|rss         |Net SW Surface Radiation                                   |1-hr       |Atmo 2D |M   |      |
|rsus        |Surface Upwelling Shortwave Radiation                      |1-hr       |Atmo 2D |M   |L     |
|rsut        |TOA Outgoing Shortwave Radiation                           |1-hr       |Atmo 2D |M   |L     |
|sfcWindmax  |Daily Maximum Near-Surface Wind Speed                      |daily      |Atmo 2D |X   |L     |
|sic         |Sea Ice Area Fraction                                      |3-hr       |Ocean 2D|I   |      |
|sit         |Sea Ice Thickness                                          |3-hr       |Ocean 2D|I   |      |
|snc         |Snow Area Fraction                                         |3-hr       |Land 2D |I   |      |
|snd         |Snow Depth                                                 |3-hr       |Land 2D |I   |L     |
|snm         |Surface Snow Melt                                          |3-hr       |Land 2D |M   |L     |
|snw         |Surface Snow Amount                                        |3-hr       |Land 2D |I   |T     |
|ta          |Air Temperature                                            |3-hr (6-hr)|Atmo 3D |I   |L     |
|tas         |Near-Surface Air Temperature                               |1-hr, daily|Atmo 2D |I   |T     |
|tasmax      |Daily Maximum Near-Surface Temperature                     |daily      |Atmo 2D |X   |T     |
|tasmin      |Daily Minimum Near-Surface Temperature                     |daily      |Atmo 2D |N   |T     |
|tke         |Turbulent Kinetic Energy                                   |monthly    |Atmo 2D |I   |      |
|ts          |Surface Temperature                                        |3-hr       |Atmo 2D |I   |L     |
|tsl         |Temperature of Soil                                        |3-hr (6-hr)|Land 3D |I   |L     |
|tsmax       |Maximum Near-Surface Temperature                           |3-hr       |Atmo 2D |X   |      |
|tsmin       |Minimum Near-Surface Temperature                           |3-hr       |Atmo 2D |N   |      |
|tso         |Sea Surface Temperature                                    |3-hr       |Ocean 2D|I   |      |
|ttop        |Air Temperature at Cloud Top                               |3-hr       |Atmo 2D |I   |      |
|ua          |Eastward Wind                                              |3-hr (6-hr)|Atmo 3D |I   |L     |
|uas         |Eastward Near-Surface Wind                                 |1-hr       |Atmo 2D |I   |T     |
|uqvc        |Vertically Integrated East. Comp. of Specific Humidity Flux|3-hr       |Atmo 2D |M   |      |
|uvmax       |Maximum Near-Surface Wind Speed                            |3-hr       |Atmo 2D |X   |      |
|va          |Northward Wind                                             |3-hr (6-hr)|Atmo 3D |I   |L     |
|vas         |Northward Near-Surface Wind                                |1-hr       |Atmo 2D |I   |T     |
|volmrfsl    |Soil Layer Volumetric Frozen Water Content                 |3-hr       |Land 3D |I   |      |
|volmrliqsl  |Soil Layer Volumetric Liquid Water Content                 |3-hr       |Land 3D |I   |      |
|vqvc        |Vertically Integrated North Comp. of Specific Humidity Flux|3-hr       |Atmo 2D |M   |      |
|wa          |Upward air velocity                                        |3-hr (6-hr)|Atmo 3D |I   |L     |
|zg          |Geopotential Height                                        |3-hr (6-hr)|Atmo 3D |I   |L     |
|zmla        |Height of Boundary Layer                                   |3-hr       |Atmo 2D |I   |L     |
