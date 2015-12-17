package smoge.main;

import smoge.species.Gene;
import smoge.species.RNApol;
import smoge.species.Ribosome;
import smoge.species.Spliceosome;

import java.io.*;
import java.util.Calendar;

/**
 * Created by fpenim on 10/12/2015.
 */
public class main {

    public static void main ( String[] args ) {

        Gene gene1 = new Gene(Integer.parseInt(args[2])); //Criar gene
        RNApol[] RNApolArray = new RNApol[Integer.parseInt(args[3])]; //Criar Array de RNApol
        Spliceosome[] SplArray = new Spliceosome[Integer.parseInt(args[4])]; //Criar Array de Splicesomas
        Ribosome[] RibArray = new Ribosome[Integer.parseInt(args[5])]; //Criar Array de Ribossomas
        //args[6] = nr locais de splicing

        //Preencher Arrays
        for (int i=0; i < RNApolArray.length; i++){
            RNApolArray[i] = new RNApol(gene1, Double.parseDouble(args[7]), Double.parseDouble(args[8]), Double.parseDouble(args[9]), Double.parseDouble(args[10]));
        }

        for (int i=0; i < SplArray.length; i++){
            SplArray[i] = new Spliceosome(RNApolArray, Double.parseDouble(args[11]), Double.parseDouble(args[12]), Double.parseDouble(args[13]), Double.parseDouble(args[14]));
        }

        for (int i=0; i < RibArray.length; i++){
            RibArray[i] = new Ribosome(SplArray, Double.parseDouble(args[15]), Double.parseDouble(args[16]), Double.parseDouble(args[17]), Double.parseDouble(args[18]));
        }

        Calendar ci = Calendar.getInstance();
        String ti = ci.getTime().toString();

        System.out.println("");
        System.out.println("Start: "+ti);
        System.out.println("------------------------------------");

        File file = new File("teste.txt");
        try{
            BufferedWriter out = new BufferedWriter(new FileWriter(file));

            int i=0;
            double tempo = 0.0;
            while(tempo<Double.parseDouble(args[0])){

                double[] ProbArray = AlgorithmManager.CalcProb(RNApolArray, SplArray, RibArray); //Calcular as probabilidades
                //Execucao do algoritmo
                double t = AlgorithmManager.SorteioT(ProbArray, gene1, RNApolArray, SplArray, RibArray, Integer.parseInt(args[6]));

                if(i%(Integer.parseInt(args[1]))==0){
                    //Recolha dos valores
                    int cp = 0;
                    for(int j=0; j<RNApolArray.length; j++){
                        cp+=RNApolArray[j].getPmRNA().size();
                    }
                    int cm = 0;
                    for(int j=0; j<SplArray.length; j++){
                        cm+=SplArray[j].getArrayM().size();
                    }
                    int c = 0;
                    for(int j=0; j<RibArray.length; j++){
                        c+=RibArray[j].getProtArray().size();
                    }

                    //Escrita no ficheiro
                    out.write(tempo+","+cp+","+cm+","+c);
                    out.newLine();
                    out.flush();
                }

                //Incrementos
                tempo+=t;
                i++;
            }
            out.close();
        }
        catch (FileNotFoundException e){
            System.err.println("Ficheiro nao encontrado!");
            System.exit(0);
        }
        catch (IOException e) {
            e.printStackTrace();
            System.exit(0);
        }

        Calendar cf = Calendar.getInstance();
        String tf = cf.getTime().toString();
        System.out.println("End: "+tf);
        long et=cf.getTimeInMillis()-ci.getTimeInMillis();
        System.out.println("------------------------------------");
        if(et<60000){
            System.out.println("Execution time: "+et/1000.0+"s");
        }
        if(et<3600000){
            System.out.println("Execution time: "+et/60000+"min"+(et%60000)/1000+"s");
        }
        else{
            System.out.println("Execution time: "+et/3600000+"h"+(et%3600000)/60000+"m"+((et%3600000)%60000)/1000+"s");
        }
    }
}

