# KPC-Toolbox Java
Java version of the KPC-Toolbox 

https://github.com/kpctoolboxteam/kpc-toolbox

Used to model network traces with Markovian Arrival Processes (MAPs).

## Creating the JAR

Build the project JAR

    mvn package -DskipTests

Then, if using the JAR for another application built with maven

    mvn install -DskipTests

To add it to your local maven repository. Check pom.xml for dependency information.

## License

See LICENSE.TXT
