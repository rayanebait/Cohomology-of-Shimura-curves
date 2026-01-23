The packages fdompres.gp and
fdomprestests.gp both depend on 

J. Rickards's Fundamental-domains-for-Shimura-Curves PARI/GP package which can be found at :
    https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves

# HOW TO RUN THE PACKAGES : 
(Can be done once for all)
In fdompres folder :
gp < configure.gp

## Running rgraph and slp packages
In fdompres folder :
gp rgraph.gp slp.gp

## Running fdompres 
In fdom folder :
gp fdom
\r runfdompres.gp


# HOW TO RUN TESTS :
After gp < configure.gp, tests are located
in the fdompres subfolder "tests/". Running
gp rgraphtests.gp slptests.gp will run a 
series of tests. 

To run fdomprestests.gp, 

In fdom folder :
gp fdom
\r runfdomprestests.gp
