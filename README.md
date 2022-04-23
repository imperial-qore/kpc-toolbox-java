# KPC-Toolbox Java
Java version of the KPC-Toolbox (MATLAB)

https://github.com/kpctoolboxteam/kpc-toolbox

Used to model network traces with Markovian Arrival Processes (MAPs).



## Compile and Run
Compile the project

    mvn compile

And execute 
        
    mvn exec:java -Dexec.mainClass=KPC.DemoRun


## Adding the JAR

Build the project JAR

    mvn package -DskipTests

Then, if using the JAR for another application built with maven

    mvn install -DskipTests

To add it to your local maven repository. Check pom.xml for dependency information.

See [DemoRun](src/main/java/KPC/DemoRun.java) for an example of how to execute a MAP fitting to a custom trace file.
## License

See LICENSE.TXT
