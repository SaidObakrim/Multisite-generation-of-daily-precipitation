## Multisite generation of daily, monthly, and annual precipitation

This project presents an implementation of the familiar two-part model (Wilks 1998) for daily simulation of precipitation
at multiple locations. And the nested model (Srikanthan 2009) for monthly and annual simulations.
As its name indicates, the two-part model is consisting of two models. The first model, the occurrence process, simulates binary
sequences which represent the precipitation occurrence (wet or dry) at each day on each location.
The second model, the amounts process, simulates non zero amounts using the Gamma distribution.
