package smoge.main;

import smoge.species.Gene;
import smoge.species.RNApol;
import smoge.species.Ribossoma;
import smoge.species.Spliceosoma;

/**
 * Created by fpenim on 10/12/2015.
 */
public class main {

    public static void main ( String[] args ) {

        Gene gene1 = new Gene(Integer.parseInt(args[2])); //Criar gene
        RNApol[] RNApolArray = new RNApol[Integer.parseInt(args[3])]; //Criar Array de RNApol
        Spliceosoma[] SplArray = new Spliceosoma[Integer.parseInt(args[4])]; //Criar Array de Splicesomas
        Ribossoma[] RibArray = new Ribossoma[Integer.parseInt(args[5])]; //Criar Array de Ribossomas
        //args[6] = nr locais de splicing

        //Preencher Arrays
        for (int i=0; i < RNApolArray.length; i++){
            RNApolArray[i] = new RNApol(gene1, Double.parseDouble(args[7]), Double.parseDouble(args[8]), Double.parseDouble(args[9]), Double.parseDouble(args[10]));
        }

        for (int i=0; i < SplArray.length; i++){
            SplArray[i] = new Spliceosoma(RNApolArray, Double.parseDouble(args[11]), Double.parseDouble(args[12]), Double.parseDouble(args[13]), Double.parseDouble(args[14]));
        }

        for (int i=0; i < RibArray.length; i++){
            RibArray[i] = new Ribossoma(SplArray, Double.parseDouble(args[15]), Double.parseDouble(args[16]), Double.parseDouble(args[17]), Double.parseDouble(args[18]));
        }
    }
}
