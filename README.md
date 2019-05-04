## Multisite generation of daily precipitation

In this project I presente an implementation of the familiar two-part model (Wilks 1998) for daily simulation of precipitation
at multiple locations.
As its name indicates, the two-part model is consisting of two models. The first model, the occurrence process, simulates binary
sequences which represent the precipitation occurrence (wet or dry) at each day on each location.
The second model, the amounts process, simulates precipitation amounts using the mixed exponential distribution.
