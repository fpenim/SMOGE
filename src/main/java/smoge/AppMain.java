package main.java.smoge;

import main.java.smoge.main.AlgorithmManager;
import main.java.smoge.species.Gene;
import main.java.smoge.species.RNApol;
import main.java.smoge.species.Ribosome;
import main.java.smoge.species.Spliceosome;
import main.java.smoge.utils.Timer;

import java.io.*;
import java.util.logging.Logger;

public class AppMain {
    private static final Logger log = Logger.getLogger(AppMain.class.getName()); //FIXME - fix logging

    public static void main ( String[] args ) {
        Timer timer = new Timer();
        log.info("Starting application...");

        Gene gene = null;
        RNApol[] polimerases = null;
        Spliceosome[] spliceosomes = null;
        Ribosome[] ribosomes = null;

        try {
            gene = new Gene(Integer.parseInt(args[2]));                  // Creating Gene
            polimerases = new RNApol[Integer.parseInt(args[3])];         // Creating RNApol array
            spliceosomes = new Spliceosome[Integer.parseInt(args[4])];   // Creating Splicesomes array
            ribosomes = new Ribosome[Integer.parseInt(args[5])];         // Creating Ribosomes array
            //args[6] = Number of splicing sites

            // Fill in arrays
            for (int i = 0; i < polimerases.length; i++) {
                polimerases[i] = new RNApol(gene, Double.parseDouble(args[7]), Double.parseDouble(args[8]), Double.parseDouble(args[9]), Double.parseDouble(args[10]));
            }

            for (int i = 0; i < spliceosomes.length; i++) {
                spliceosomes[i] = new Spliceosome(polimerases, Double.parseDouble(args[11]), Double.parseDouble(args[12]), Double.parseDouble(args[13]), Double.parseDouble(args[14]));
            }

            for (int i = 0; i < ribosomes.length; i++) {
                ribosomes[i] = new Ribosome(spliceosomes, Double.parseDouble(args[15]), Double.parseDouble(args[16]), Double.parseDouble(args[17]), Double.parseDouble(args[18]));
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            e.printStackTrace();
            log.info("An error occurred while initializing arrays, elapsed time: " + timer.getElapsedTime());
            System.exit(1);
        }

        log.info("All arrays were initialized successfully.");

        File file = new File("teste.txt");
        try (BufferedWriter out = new BufferedWriter(new FileWriter(file))) {

            int i = 0;
            double tempo = 0.0;
            while(tempo < Double.parseDouble(args[0])){

                double[] ProbArray = AlgorithmManager.calculateProbabilities(polimerases, spliceosomes, ribosomes);                                   // Probabilities calculus
                double t = AlgorithmManager.firstRaffle(ProbArray, gene, polimerases, spliceosomes, ribosomes, Integer.parseInt(args[6])); // Algorithm execution

                if(i % (Integer.parseInt(args[1])) == 0) { // Values output
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
                    out.write(tempo+","+cp+","+cm+","+c);
                    out.newLine();
                    out.flush();
                }
                tempo += t;
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

        log.info("Application terminated, elapsed time: " + timer.getElapsedTime());
    }
}