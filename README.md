The packages fdompres.gp and
fdomprestests.gp both depend on J. Rickards's

    Fundamental-domains-for-Shimura-Curves

PARI/GP package which can be found at :

    https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves

In the following denote the folder in which it is
downloaded by fdom. Also denote the current package's folder
by fdompres.

# Recommendation
For both the fdompres and fdom packages, it is recommended to
run the line
    
    default(parisizemax, "1G")

or to put 

    parisizemax="1G"

in your gprc which can be created at $HOME/.gprc 

# Dependency tree
The packages are organized as follows :

    slp.gp ->slptests.gp
        |-> rgraph.gp -> rgraphtests.gp
        |          |            
    rgraphutils.gp -> fdompres.gp -> fdomprestests.gp

# The timings folder
It contains the code ran for the timings in the article.

# HOW TO RUN THE PACKAGES : 
(Can be done once for all)
In fdompres folder :

    gp < configure.gp

the command

    gp loadpackages.gp

then loads every fdompres package.

## Running fdompres 
As fdompres is dependant on fdom it has to be run
as follows. In fdom folder :

    gp fdom
    \r runfdompres.gp


# HOW TO RUN TESTS :
After gp < configure.gp, tests are located
in the fdompres subfolder "tests/". Running

    \r rgraphtests.gp slptests.gp

will run a series of tests. The packages rgraph.gp and
slp.gp need to be loaded for those to work.

To run fdomprestests.gp, in fdom folder :

    gp fdom
    \r runfdomprestests.gp 

(Sometimes need to be run twice for obscure reasons.)
