import org.apache.log4j.Logger;
import smoge.managers.PropertiesManager;
import smoge.species.Gene;
import smoge.species.RNApolymerase;
import smoge.species.Ribosome;
import smoge.species.Spliceosome;


public class SmogeApplication {
    public static final Logger log = Logger.getLogger(SmogeApplication.class);

    private static Gene gene = null;
    private static RNApolymerase[] polimerases = null;
    private static Spliceosome[] spliceosomes = null;
    private static Ribosome[] ribosomes = null;


    public static void main ( String[] args ) {
        log.info("Starting application...");

        PropertiesManager propertiesManager = PropertiesManager.getInstance();

        initializeEnvironment(propertiesManager);

        //startSimulation();

        log.info("Application terminated.");
    }

    private static void initializeEnvironment(PropertiesManager properties) {
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
/*
    private static void startSimulation() {
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
    }
*/
}