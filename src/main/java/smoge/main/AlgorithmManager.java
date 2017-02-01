package main.java.smoge.main;

import main.java.smoge.species.*;

import java.util.Random;

public class AlgorithmManager {
    private int time;
    private int pmRNA;
    private int mRNA;
    private int prot;

    public AlgorithmManager(){
        time = 0;
        pmRNA = 0;
        mRNA = 0;
        prot = 0;
    }

    public static double[] calculateProbabilities(RNApol[] polimerases, Spliceosome[] spliceosomes, Ribosome[] ribosomes) {
        double[] probabilitiesArray = new double[3];

        // sum of k - RNApolimerases
        double sr = 0.0;
        for(int i = 0; i < polimerases.length; i++){
            sr += polimerases[i].sum();
        }

        // sum of k - spliceosomes
        double ss = 0.0;
        for(int i = 0; i < spliceosomes.length; i++){
            ss += spliceosomes[i].sum();
        }

        // sum of k - ribosomes
        double sri = 0.0;
        for(int i = 0; i < ribosomes.length; i++){
            sri += ribosomes[i].sum();}

        double t = sr + ss + sri;
        probabilitiesArray[0] = sr / t;
        probabilitiesArray[1] = (sr + ss) / t;
        probabilitiesArray[2] = (sr + ss + sri) / t;

        return probabilitiesArray;
    }

    public static double sum(RNApol[] rna, Spliceosome[] spl, Ribosome[] rib){
        double sm = 0.0;

        for(int i = 0; i < rna.length; i++) {
            sm += rna[i].sum();
        }
        for(int i = 0; i < spl.length; i++) {
            sm += spl[i].sum();
        }
        for(int i = 0; i < rib.length; i++) {
            sm += rib[i].sum();
        }
        return sm;
    }

    public static double firstRaffle(double[] prob, Gene gene, RNApol[] rna, Spliceosome[] spl, Ribosome[] rib, int s) { // Raffle between: RNApolimerase, Spliceosome or Ribosome
        Random pseudoRandom = new Random();
        double randomDouble = pseudoRandom.nextDouble();

        double sum = sum(rna, spl, rib);
        double time = Math.log(1 - randomDouble) / (-sum); // Exponential distribution - time

        if(randomDouble <= prob[0]) { // -> RNApolimerase
            secondRaffle(rna, gene, rna, spl, s, gene.getGeneLength());

        } else if(randomDouble <= prob[1]) { // -> Spliceosome
            secondRaffle(spl, gene, rna, spl, s, gene.getGeneLength());

        } else { // -> Ribosome
            secondRaffle(rib, gene, rna, spl, s, gene.getGeneLength());
        }
        return time;
    }

    public static void secondRaffle(Element[] ele, Gene gene, RNApol[] rna, Spliceosome[] spl, int s, int dg) { // Given the result of the firstRaffle, raffle the specific element
        double[] prob = new double[ele.length]; // Probabilities array

        for(int i = 0; i < ele.length; i++) { // Individual probability
            prob[i] = ele[i].sum();
        }

        for(int i = 1; i < prob.length;i++) { // Cumulative probability
            prob[i] += prob[i-1];
        }

        for(int i = 0; i < prob.length; i++) { // Division by the total K
            prob[i] /= prob[prob.length - 1];
        }

        // Raffle
        Random pseudoRandom = new Random();
        double randomDouble = pseudoRandom.nextDouble();
        int i = 0;
        while(randomDouble > prob[i]) {
            i++;
        }
        thirdRaffle(ele[i], gene, rna, spl, s, dg);
    }

    public static void thirdRaffle(Element element, Gene gene, RNApol[] rna, Spliceosome[] spl, int s, int dg) { // Raffle the action
        double[] probabilities = new double[4]; // Probabilities array
        Random pseudoRandom = new Random();

        if(element instanceof RNApol) {
            polimeraseActionRaffle((RNApol) element, gene, s, probabilities, pseudoRandom);
        }
        else if(element instanceof Spliceosome) {
            spliceosomeActionRaffle((Spliceosome) element, rna, dg, probabilities, pseudoRandom);
        }
        else {
            ribosomeActionRaffle((Ribosome) element, spl, probabilities, pseudoRandom);
        }
    }

    private static void polimeraseActionRaffle(RNApol element, Gene gene, int s, double[] probabilities, Random pseudoRandom) {
        double aux = element.sum();
        probabilities[0] = element.getKc() / aux;                                                           // Connect
        probabilities[1] = (element.getKc() + element.getKp()) / aux;                                       // Slide
        probabilities[2] = (element.getKc() + element.getKp() + element.getKd()) / aux;                     // Disconnect
        probabilities[3] = (element.getKc() + element.getKp() + element.getKd() + element.getKdg()) / aux;  // Degrade

        double randomDouble = pseudoRandom.nextDouble();
        int i = 0;
        while(randomDouble > probabilities[i]) { i++; }
        if(i == 0) {
            element.connect(gene, s);
        } else if(i == 1) {
            element.progress(s);
        } else if(i == 2) {
            element.Disconnect();
        } else {
            element.degradePmRNA();
        }
    }

    private static void spliceosomeActionRaffle(Spliceosome element, RNApol[] rna, int dg, double[] probabilities, Random pseudoRandom) {
        double aux = element.sum();
        probabilities[0]= element.getKc() / aux;                                                          // Connect
        probabilities[1]= (element.getKc() + element.getKs()) / aux;                                      // Splice
        probabilities[2]= (element.getKc() + element.getKs() + element.getKt()) / aux;                    // Transport
        probabilities[3]= (element.getKc() + element.getKs() + element.getKt() + element.getKdg()) / aux; // Degrade

        double randomDouble = pseudoRandom.nextDouble();
        int i = 0;
        while(randomDouble > probabilities[i]){ i++; }
        if(i == 0) {
            for(int a = 0; a < rna.length; a++) {
                for(int b = rna[a].getPmRNA().size() - 1; b >- 1; b--){
                    element.connect(rna[a].getPmRNA().get(b));
                }
            }
        } else if(i == 1){
            element.splice();
        } else if(i==2){
            element.transport(dg);
        } else {
            element.degrademRNA();
        }
    }

    private static void ribosomeActionRaffle(Ribosome element, Spliceosome[] spl, double[] probabilities, Random pseudoRandom) {
        double aux = element.sum();
        probabilities[0] = element.getKc() / aux;                                                           // Connect
        probabilities[1] = (element.getKc() + element.getKp()) / aux;                                       // Slide
        probabilities[2] = (element.getKc() + element.getKp() + element.getKd()) / aux;                     // Disconnect
        probabilities[3] = (element.getKc() + element.getKp() + element.getKd() + element.getKdg()) / aux;  // Degrade protein (!!! - for code simplicity)

        double randomDouble = pseudoRandom.nextDouble();
        int i = 0;
        while(randomDouble > probabilities[i]){ i++; }
        if(i == 0){
            for(int a = 0;a < spl.length;a++){
                for(int b = spl[a].getArrayM().size() - 1; b >- 1; b--){
                    element.connect(spl[a].getArrayM().get(b));
                }
            }
        } else if(i == 1){
            element.progress();
        } else if(i == 2){
            element.disconnect();
        } else {
            element.degradeProtein();
        }
    }
}