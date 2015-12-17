package smoge.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Spliceosoma extends Element {

    private boolean lig; //Se esta ligada a um pre-mRNA (sim=true ; n�o=false)
    private premRNA currentP; //pr�-mRNA a q esta ligado
    private ArrayList<mRNA> arrayM;
    private RNApol[] polRef; //Refer�ncia �s RNApol

    //Constantes cineticas
    private double kl; //Ligar
    private double ks; //Evento splicing
    private double kt; //Transporte -> nuc-cit
    private double kd; //Degrada��o de mRNA

    /**
     * Construtor da Classe
     *
     */
    public Spliceosoma(RNApol[] rna, double kl, double ks, double kt, double kdg){
        this.kl = kl;
        this.ks = ks;
        this.kt = kt;
        this.kd = kdg;
        lig = false;
        currentP = null;
        arrayM = new ArrayList<mRNA>();
        polRef = rna;
    }

    //-----------------------------------------A��es do Spliceosoma -----------------------------------------//
    //Ligar
    public boolean ligar(premRNA pmR){
        if(!lig && !pmR.isSpliceosoma() && !pmR.isCitosol()){ //Se o Spliceosoma n�o estiver ligado; pr�-mRNA n�o estiver ligado
            pmR.ligarSpl(this); //Altera argumento do pmRNA - tem um spliceosoma ligado
            lig = true; //Passa a estar ligado
            currentP = pmR; //pr�-mRNA a que est� ligado
            //kl = 0.0; // N�o se pode ligar (ja esta)
            //ks = 125.5; //Pode ocorrer splicing
            return true;
        }
        return false;
    }

    //Splicing
    public void Splicing(){
        boolean a =currentP.splicing();//Ocorre splicing
		/*if(a){//Se todos os locais de splicing ja tiverem ocorrido
			//ks = 0.0; //Nao pode ocorrer mais splicing
			//kt = 533; //Pode correr transporte
		}*/
    }

    //Transportar n-c
    public void transportar(int dg){
        if(lig && currentP.getSplSitesTT()==currentP.getSplSitesD() && !currentP.isRnaP()){
            currentP.transporte();
            desligarS();

            mRNA m = new mRNA(dg);
            arrayM.add(0, m);
            //System.out.println("-----------------Transporte");
        }
    }

    //Desligar
    public void desligarS(){
        currentP = null; //N�o est� ligado a nenhum pr�-mRNA
        lig = false; //N�o est� ligado
        //kl = 56.2; //Pode ligar-se
        //ks = 0.0; //N�o pode 'fazer' splicing
        //kt = 0.0; //N�o pode transportar
    }

    //--------------------------------------------Degrada��o--------------------------------------------//
    public void degradarmRNA(){
        if(!arrayM.isEmpty()){
            mRNA m;
            m = arrayM.get(arrayM.size()-1); //mRNA mais antigo
            if(m.getRibossomasLig().isEmpty()){//Se n�o tiver nenhum ribossoma ligado
                arrayM.remove(m);
            }
            else{
                for(int i=0;i<m.getRibossomasLig().size();i++){ // percorre a lista de todos os ribossomas ligados
                    m.getRibossomasLig().get(i).Desligar(); //desliga os ribossomas
                }
                arrayM.remove(m);
            }
        }
    }

    //-----------------------------------------Soma - Algor�tmo-----------------------------------------//

    @Override
    public double sum(){
        double sm = getKl()+ getKs() + getKt() + getKd();
        return sm;
    }

    //---------------------------------------Getters and Setters---------------------------------------//
    public double getKl() {
        int a=0;
        for(int i=0; i<polRef.length; i++){
            if(!polRef[i].getPmRNA().isEmpty()){
                for(int j=0; j<polRef[i].getPmRNA().size(); j++){
                    if(!polRef[i].getPmRNA().get(j).isSpliceosoma() && !polRef[i].getPmRNA().get(j).isCitosol()){
                        a+=1;
                    }
                }
            }
        }
        double sm = kl*((double)a);
        return sm;
    }

    public double getKs() {
        if(lig && currentP.getSplSitesA()>0){
            return ks*((double)currentP.getSplSitesA());
        }
        return ks*0;
    }

    public double getKt() {
        if(lig && currentP.getSplSitesTT()==currentP.getSplSitesD() && !currentP.isRnaP())
            return kt;
        return 0.0;
    }

    public double getKd() {
        if(arrayM.isEmpty())
            return kd*(double)0;
        return kd*(arrayM.size());
    }

    public void setKl(double kl) {
        this.kl = kl;
    }

    public void setKs(double ks) {
        this.ks = ks;
    }

    public void setKt(double kt) {
        this.kt = kt;
    }

    public void setKd(double kd) {
        this.kd = kd;
    }

    public boolean isLig() {
        return lig;
    }

    public void setLig(boolean lig) {
        this.lig = lig;
    }

    public premRNA getCurrentP() {
        return currentP;
    }

    public void setCurrentP(premRNA currentP) {
        this.currentP = currentP;
    }

    public ArrayList<mRNA> getArrayM() {
        return arrayM;
    }

    public void setArrayM(ArrayList<mRNA> arrayM) {
        this.arrayM = arrayM;
    }

    //----------------------------------------- Override -----------------------------------------//
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = super.hashCode();
        result = prime * result + ((arrayM == null) ? 0 : arrayM.hashCode());
        result = prime * result
                + ((currentP == null) ? 0 : currentP.hashCode());
        long temp;
        temp = Double.doubleToLongBits(kd);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kl);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(ks);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kt);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        result = prime * result + (lig ? 1231 : 1237);
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (!super.equals(obj))
            return false;
        if (getClass() != obj.getClass())
            return false;
        Spliceosoma other = (Spliceosoma) obj;
        if (arrayM == null) {
            if (other.arrayM != null)
                return false;
        } else if (!arrayM.equals(other.arrayM))
            return false;
        if (currentP == null) {
            if (other.currentP != null)
                return false;
        } else if (!currentP.equals(other.currentP))
            return false;
        if (Double.doubleToLongBits(kd) != Double.doubleToLongBits(other.kd))
            return false;
        if (Double.doubleToLongBits(kl) != Double.doubleToLongBits(other.kl))
            return false;
        if (Double.doubleToLongBits(ks) != Double.doubleToLongBits(other.ks))
            return false;
        if (Double.doubleToLongBits(kt) != Double.doubleToLongBits(other.kt))
            return false;
        if (lig != other.lig)
            return false;
        return true;
    }
}
