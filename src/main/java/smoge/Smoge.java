package smoge;

import org.apache.log4j.Logger;
import smoge.managers.AlgorithmManager;
import smoge.managers.PropertiesManager;
import smoge.species.Gene;
import smoge.species.RNApolymerase;
import smoge.species.Ribosome;
import smoge.species.Spliceosome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

public class Smoge {
    public static final Logger log = Logger.getLogger(Smoge.class);

    private static Gene gene = null;
    private static RNApolymerase[] polimerases = null;
    private static Spliceosome[] spliceosomes = null;
    private static Ribosome[] ribosomes = null;


    public static void initializeEnvironment(PropertiesManager properties) {
        log.info("Initializing objects.");

        gene = new Gene(properties.getGeneLength());
        polimerases = new RNApolymerase[properties.getRnaPolNumber()];
        spliceosomes = new Spliceosome[properties.getSpliceosomeNumber()];
        ribosomes = new Ribosome[properties.getRibosomeNumber()];
        //args[6] = Number of splicing sites

        for (int i = 0; i < polimerases.length; i++) {
            polimerases[i] = new RNApolymerase(
                    gene,
                    properties.getPolymeraseKc(),
                    properties.getPolymeraseKp(),
                    properties.getPolymeraseKd(),
                    properties.getPolymeraseKdg()
            );
        }

        for (int i = 0; i < spliceosomes.length; i++) {
            spliceosomes[i] = new Spliceosome(
                    polimerases,
                    properties.getSpliceosomeKc(),
                    properties.getSpliceosomeKs(),
                    properties.getSpliceosomeKt(),
                    properties.getSpliceosomeKdg()
            );
        }

        for (int i = 0; i < ribosomes.length; i++) {
            ribosomes[i] = new Ribosome(
                    spliceosomes,
                    properties.getRibosomeKc(),
                    properties.getRibosomeKp(),
                    properties.getRibosomeKd(),
                    properties.getRibosomeKdg()
            );
        }

        log.info("All objects were initialized successfully.");
    }

    public static void startSimulation(PropertiesManager properties) {

        try (BufferedWriter out = new BufferedWriter(new FileWriter(new File("output.txt")))) {

            long i = 0;
            long time = 0L;
            while(time < properties.getSimulationTime()) {

                double[] ProbArray = AlgorithmManager.calculateProbabilities(polimerases, spliceosomes, ribosomes); // Probabilities calculus

                double t = AlgorithmManager.firstRaffle(ProbArray, gene, polimerases, spliceosomes, ribosomes, properties.getSpliceSitesNumber());

                if(i % properties.getOutputStep() == 0) { // Print output
                    int cp = 0;
                    for(int j = 0; j < polimerases.length; j++){
                        cp += polimerases[j].getPmRNA().size();
                    }
                    int cm = 0;
                    for(int j = 0; j < spliceosomes.length; j++){
                        cm += spliceosomes[j].getArrayM().size();
                    }
                    int c = 0;
                    for(int j = 0; j < ribosomes.length; j++){
                        c += ribosomes[j].getProtArray().size();
                    }
                    // Writing to file
                    out.write(time + "," + cp + "," + cm + "," + c);
                    out.newLine();
                    out.flush();
                }
                time += t;
                i++;
            }
        }
        catch (FileNotFoundException e){
            System.err.println("File not found!");
            System.exit(1);
        }
        catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}