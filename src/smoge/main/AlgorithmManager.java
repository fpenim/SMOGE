package smoge.main;

import smoge.species.*;

import java.util.Random;

/**
 * Created by fpenim on 12/12/2015.
 */
public class AlgorithmManager {
    private int time;
    private int pmRNA;
    private int mRNA;
    private int prot;

    /**
     * Constructor
     */
    public AlgorithmManager(){
        time = 0;
        pmRNA = 0;
        mRNA = 0;
        prot = 0;
    }

    //------------------------------------------------------------------------------------------------//

    //C�lculo das probabilidades Comulativas (RNApol/Spliceosome/Ribosome)
    public static double[] CalcProb(RNApol[] rna, Spliceosome[] spl, Ribosome[] rib){
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

    public static double Somatorio(RNApol[] rna, Spliceosome[] spl, Ribosome[] rib){
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
    public static double SorteioT(double[] prob, Gene gen, RNApol[] rna, Spliceosome[] spl, Ribosome[] rib, int s){

        Random r = new Random();
        double rn = r.nextDouble();

        double sm = Somatorio(rna, spl, rib);
        double t = Math.log(1-rn)/(-sm); //Distribui��o Exponencial - time

        if(rn<=prob[0]){
            //Sorteio RNApol (2)
            SorteioQ(rna, gen, rna, spl, s, gen.getGeneLength());
            //System.out.println("SorteioT (1)->RNApol \n");
        }else if(rn<=prob[1]){
            //Sorteio Spliceosome (2)
            SorteioQ(spl, gen, rna, spl, s, gen.getGeneLength());
            //System.out.println("SorteioT (1)->Spliceosome \n");
        }else{
            //Sorteio Ribosome (2)
            SorteioQ(rib, gen, rna, spl, s, gen.getGeneLength());
            //System.out.println("SorteioT (1)->Ribosome \n");
        }
        return t;
    }

    //------------ ---//
    //Sorteio Qual (2)
    //----------------//
    public static void SorteioQ(Element[] ele, Gene gen, RNApol[] rna, Spliceosome[] spl, int s, int dg){
        double[] prob = new double[ele.length]; //Array das probabilidades
        Random r = new Random();

        //C�lculo Individual
        for(int i=0; i<ele.length; i++){
            prob[i]=ele[i].sum();
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
    public static void SorteioA(Element ele, Gene gen, RNApol[] rna, Spliceosome[] spl, int s, int dg){
        double[] p = new double[4]; //Array das probabildades
        Random r = new Random(); //Objeto aleat�rio

        //Preenchimento do array das probabilidades
        if(ele instanceof RNApol){ //Se sorteada RNApol
            double aux = ((RNApol) ele).sum();
            p[0]=((RNApol) ele).getKc()/aux; //connect
            p[1]=(((RNApol) ele).getKc()+((RNApol) ele).getKp())/aux; //progress�ar
            p[2]=(((RNApol) ele).getKc()+((RNApol) ele).getKp()+((RNApol) ele).getKd())/aux; //disconnect
            p[3]=(((RNApol) ele).getKc()+((RNApol) ele).getKp()+((RNApol) ele).getKd()+((RNApol) ele).getKdg())/aux; //Degradar
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
                //System.out.println("SorteioA (3.1)->connect");
                ((RNApol) ele).connect(gen, s);
            }else if(i==1){
                //System.out.println("SorteioA (3.1)->progress�ar");
                ((RNApol) ele).progress(s);
            }else if(i==2){
                //System.out.println("SorteioA (3.1)->disconnect");
                ((RNApol) ele).Disconnect();
            }else{
                //System.out.println("SorteioA (3.1)->Degradar pre-mRNA");
                ((RNApol) ele).degradePmRNA();
            }
        }
        else if(ele instanceof Spliceosome){ //Se sorteado Spliceosome
            double aux = ((Spliceosome) ele).sum();
            p[0]=((Spliceosome) ele).getKc()/aux; //connect
            p[1]=(((Spliceosome) ele).getKc()+((Spliceosome) ele).getKs())/aux; //splice
            p[2]=(((Spliceosome) ele).getKc()+((Spliceosome) ele).getKs()+((Spliceosome) ele).getKt())/aux; //Transporte
            p[3]=(((Spliceosome) ele).getKc()+((Spliceosome) ele).getKs()+((Spliceosome) ele).getKt()+((Spliceosome) ele).getKdg())/aux;
            //Sorteio (3.2)
            double rn = r.nextDouble();
			/*System.out.println("Random: "+rn);
			for(int j=0; j<p.length; j++){
				System.out.println(p[j]);
			}*/
            int i = 0;
            while(rn>p[i]){
                i++;}
            if(i==0){ //connect
                //System.out.println("SorteioA (3.2)->connect");
                for(int a=0; a<rna.length; a++){
                    for(int b=rna[a].getPmRNA().size()-1; b>-1; b--){
                        ((Spliceosome)ele).connect(rna[a].getPmRNA().get(b));
                    }
                }
            }else if(i==1){ //splice
                //System.out.println("SorteioA (3.2)->splice");
                ((Spliceosome)ele).splice();
            }else if(i==2){ //Transporte
                //System.out.println("SorteioA (3.2)->Transporte");
                ((Spliceosome) ele).transport(dg);
            }else{
                //System.out.println("SorteioA (3.2)->Degradar mRNA");
                ((Spliceosome) ele).degrademRNA();
            }
        }
        else{ //Se sorteado Ribosome
            double aux = ((Ribosome) ele).sum();
            p[0]=((Ribosome)ele).getKc()/aux; //connect
            p[1]=(((Ribosome)ele).getKc()+((Ribosome)ele).getKp())/aux; //progress�ar
            p[2]=(((Ribosome)ele).getKc()+((Ribosome)ele).getKp()+((Ribosome)ele).getKd())/aux; //disconnect
            p[3]=(((Ribosome)ele).getKc()+((Ribosome)ele).getKp()+((Ribosome)ele).getKd()+((Ribosome)ele).getKdg())/aux; //Degradar Prote�na
            //Sorteio (3.3)
            double rr = r.nextDouble();
			/*System.out.println("Random: "+rr);
			for(int j=0; j<p.length; j++){
				System.out.println(p[j]);
			}*/
            int i = 0;
            while(rr>p[i]){
                i++;}
            if(i==0){ //connect
                //System.out.println("SorteioA (3.3)->connect");
                for(int a=0;a<spl.length;a++){
                    for(int b=spl[a].getArrayM().size()-1; b>-1; b--){
                        ((Ribosome)ele).connect(spl[a].getArrayM().get(b));
                    }
                }
            }else if(i==1){ //progress�ar
                //System.out.println("SorteioA (3.3)->progress�ar");
                ((Ribosome)ele).progress();
            }else if(i==2){
                //System.out.println("SorteioA (3.3)->disconnect");
                ((Ribosome)ele).disconnect();
            }else{
                //System.out.println("SorteioA (3.3)->Degradar Prote�na");
                ((Ribosome)ele).degradeProtein();
            }
        }
    }

    //-----------------------------------------Getters & Setters-----------------------------------------//

    public int getTime() {
        return time;
    }

    public void setTime(int time) {
        this.time = time;
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
