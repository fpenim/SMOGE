package smoge.main;

import smoge.species.Gene;
import smoge.species.RNApol;
import smoge.species.Ribossoma;
import smoge.species.Spliceosoma;

import java.util.Random;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Relogio {
    private int tempo;
    private int pmRNA;
    private int mRNA;
    private int prot;

    /**
     * Construtor da Classe
     */
    public Relogio(){
        tempo = 0;
        pmRNA = 0;
        mRNA = 0;
        prot = 0;
    }

    //------------------------------------------------------------------------------------------------//

    //C�lculo das probabilidades Comulativas (RNApol/Spliceosoma/Ribossoma)
    public static double[] CalcProb(RNApol[] rna, Spliceosoma[] spl, Ribossoma[] rib){
        double[] ProbArray = new double[3]; //Array das probabilidades

        //Soma dos k das RNApol
        double sr = 0.0;
        for(int i=0; i<rna.length; i++){
            sr += rna[i].sum();}
        //System.out.println();
        //Soma dos k dos spliceosomas
        double ss = 0.0;
        for(int i=0; i<spl.length; i++){
            ss += spl[i].sum();}
        //System.out.println();
        //Soma dos k dos ribossomas
        double sri = 0.0;
        for(int i=0; i<rib.length; i++){
            sri += rib[i].sum();}
        //System.out.println();
        //calculo das probabilidades
        double t = sr+ss+sri;
        ProbArray[0] = sr/t;
        ProbArray[1] = (sr+ss)/t;
        ProbArray[2] = (sr+ss+sri)/t;

        return ProbArray;
    }

    public static double Somatorio(RNApol[] rna, Spliceosoma[] spl, Ribossoma[] rib){
        double sm = 0.0;
        for(int i=0; i<rna.length; i++)
            sm+=rna[i].sum();
        for(int i=0; i<spl.length; i++)
            sm+=spl[i].sum();
        for(int i=0; i<rib.length; i++)
            sm+=rib[i].sum();
        return sm;
    }
    //-----------------------------------------Sorteio T|Q|A -----------------------------------------//

    //Sorteio Tipo (1)
    /**
     *
     * @param prob Probabilidades
     * @param gen  Gene
     * @param rna  Array de RNApol
     * @param spl  Array de Spliceosomas
     * @param rib  Array de Ribossomas
     * @param s    Nr de locais de splicing
     * @return
     */
    public static double SorteioT(double[] prob, Gene gen, RNApol[] rna, Spliceosoma[] spl, Ribossoma[] rib, int s){

        Random r = new Random();
        double rn = r.nextDouble();

        double sm = Somatorio(rna, spl, rib);
        double t = Math.log(1-rn)/(-sm); //Distribui��o Exponencial - tempo

        if(rn<=prob[0]){
            //Sorteio RNApol (2)
            SorteioQ(rna, gen, rna, spl, s, gen.getDimGene());
            //System.out.println("SorteioT (1)->RNApol \n");
        }else if(rn<=prob[1]){
            //Sorteio Spliceosoma (2)
            SorteioQ(spl, gen, rna, spl, s, gen.getDimGene());
            //System.out.println("SorteioT (1)->Spliceosoma \n");
        }else{
            //Sorteio Ribossoma (2)
            SorteioQ(rib, gen, rna, spl, s, gen.getDimGene());
            //System.out.println("SorteioT (1)->Ribossoma \n");
        }
        return t;
    }

    //------------ ---//
    //Sorteio Qual (2)
    //----------------//
    public static void SorteioQ(Elemento[] ele, Gene gen, RNApol[] rna, Spliceosoma[] spl, int s, int dg){
        double[] prob = new double[ele.length]; //Array das probabilidades
        Random r = new Random();

        //C�lculo Individual
        for(int i=0; i<ele.length; i++){
            prob[i]=ele[i].soma();
        }
        //C�lculo Cumulativo - Soma cumulativa
        for(int i=1; i<prob.length;i++){
            prob[i]+=prob[i-1];
        }
        //C�lculo da Probabilidades - Divis�o pelo K total
        for(int i=0; i<prob.length; i++){
            prob[i]/=prob[prob.length-1];
            //System.out.println(prob[i]);
        }
        //Sorteio
        double rn = r.nextDouble();
        int i =0;
        while(rn>prob[i]){
            i++;}
        //System.out.println(rn);
        //System.out.println("SorteioQ (2)->Element i= "+i);
        //Sorteio A (3)
        SorteioA(ele[i], gen, rna, spl, s, dg);
    }
    //----------------//
    //Sorteio Açao (3)
    //----------------//
    public static void SorteioA(Elemento ele, Gene gen, RNApol[] rna, Spliceosoma[] spl, int s, int dg){
        double[] p = new double[4]; //Array das probabildades
        Random r = new Random(); //Objeto aleat�rio

        //Preenchimento do array das probabilidades
        if(ele instanceof RNApol){ //Se sorteada RNApol
            double aux = ((RNApol) ele).sum();
            p[0]=((RNApol) ele).getKl()/aux; //Ligar
            p[1]=(((RNApol) ele).getKl()+((RNApol) ele).getKa())/aux; //Avan�ar
            p[2]=(((RNApol) ele).getKl()+((RNApol) ele).getKa()+((RNApol) ele).getKd())/aux; //Desligar
            p[3]=(((RNApol) ele).getKl()+((RNApol) ele).getKa()+((RNApol) ele).getKd()+((RNApol) ele).getKdg())/aux; //Degradar
            //Sorteio (3.1)
            double rn = r.nextDouble();
			/*System.out.println("Random: "+rn);
			for(int j=0; j<p.length; j++){
				System.out.println(p[j]);
			}*/
            int i = 0;
            while(rn>p[i]){
                i++;}
            if(i==0){
                //System.out.println("SorteioA (3.1)->Ligar");
                ((RNApol) ele).Ligar(gen, s);
            }else if(i==1){
                //System.out.println("SorteioA (3.1)->Avan�ar");
                ((RNApol) ele).Avan(s);
            }else if(i==2){
                //System.out.println("SorteioA (3.1)->Desligar");
                ((RNApol) ele).Desligar();
            }else{
                //System.out.println("SorteioA (3.1)->Degradar pre-mRNA");
                ((RNApol) ele).degradarPmRNA();
            }
        }
        else if(ele instanceof Spliceosoma){ //Se sorteado Spliceosoma
            double aux = ((Spliceosoma) ele).sum();
            p[0]=((Spliceosoma) ele).getKl()/aux; //Ligar
            p[1]=(((Spliceosoma) ele).getKl()+((Spliceosoma) ele).getKs())/aux; //Splicing
            p[2]=(((Spliceosoma) ele).getKl()+((Spliceosoma) ele).getKs()+((Spliceosoma) ele).getKt())/aux; //Transporte
            p[3]=(((Spliceosoma) ele).getKl()+((Spliceosoma) ele).getKs()+((Spliceosoma) ele).getKt()+((Spliceosoma) ele).getKd())/aux;
            //Sorteio (3.2)
            double rn = r.nextDouble();
			/*System.out.println("Random: "+rn);
			for(int j=0; j<p.length; j++){
				System.out.println(p[j]);
			}*/
            int i = 0;
            while(rn>p[i]){
                i++;}
            if(i==0){ //Ligar
                //System.out.println("SorteioA (3.2)->Ligar");
                for(int a=0; a<rna.length; a++){
                    for(int b=rna[a].getPmRNA().size()-1; b>-1; b--){
                        ((Spliceosoma)ele).ligar(rna[a].getPmRNA().get(b));
                    }
                }
            }else if(i==1){ //Splicing
                //System.out.println("SorteioA (3.2)->Splicing");
                ((Spliceosoma)ele).Splicing();
            }else if(i==2){ //Transporte
                //System.out.println("SorteioA (3.2)->Transporte");
                ((Spliceosoma) ele).transportar(dg);
            }else{
                //System.out.println("SorteioA (3.2)->Degradar mRNA");
                ((Spliceosoma) ele).degradarmRNA();
            }
        }
        else{ //Se sorteado Ribossoma
            double aux = ((Ribossoma) ele).sum();
            p[0]=((Ribossoma)ele).getKl()/aux; //Ligar
            p[1]=(((Ribossoma)ele).getKl()+((Ribossoma)ele).getKa())/aux; //Avan�ar
            p[2]=(((Ribossoma)ele).getKl()+((Ribossoma)ele).getKa()+((Ribossoma)ele).getKd())/aux; //Desligar
            p[3]=(((Ribossoma)ele).getKl()+((Ribossoma)ele).getKa()+((Ribossoma)ele).getKd()+((Ribossoma)ele).getKdg())/aux; //Degradar Prote�na
            //Sorteio (3.3)
            double rr = r.nextDouble();
			/*System.out.println("Random: "+rr);
			for(int j=0; j<p.length; j++){
				System.out.println(p[j]);
			}*/
            int i = 0;
            while(rr>p[i]){
                i++;}
            if(i==0){ //Ligar
                //System.out.println("SorteioA (3.3)->Ligar");
                for(int a=0;a<spl.length;a++){
                    for(int b=spl[a].getArrayM().size()-1; b>-1; b--){
                        ((Ribossoma)ele).Ligar(spl[a].getArrayM().get(b));
                    }
                }
            }else if(i==1){ //Avan�ar
                //System.out.println("SorteioA (3.3)->Avan�ar");
                ((Ribossoma)ele).Avan();
            }else if(i==2){
                //System.out.println("SorteioA (3.3)->Desligar");
                ((Ribossoma)ele).Desligar();
            }else{
                //System.out.println("SorteioA (3.3)->Degradar Prote�na");
                ((Ribossoma)ele).degradarProt();
            }
        }
    }

    //-----------------------------------------Getters & Setters-----------------------------------------//

    public int getTempo() {
        return tempo;
    }

    public void setTempo(int tempo) {
        this.tempo = tempo;
    }

    public int getPmRNA() {
        return pmRNA;
    }

    public void setPmRNA(int pmRNA) {
        this.pmRNA = pmRNA;
    }

    public int getmRNA() {
        return mRNA;
    }

    public void setmRNA(int mRNA) {
        this.mRNA = mRNA;
    }

    public int getProt() {
        return prot;
    }

    public void setProt(int prot) {
        this.prot = prot;
    }
}
