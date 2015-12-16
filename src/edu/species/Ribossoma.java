package edu.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Ribossoma extends Element {
    private int velTrad; //Velocidade de Tradu��o
    private int posRib; //Posi��o do Ribossoma
    private int dimRib; //Dimens�o do Ribossoma
    private boolean ligado = false; //Se est� ligado
    private mRNA currentM; //A qual est� ligado
    private Proteina currentProt; //Prote�na em s�ntese
    private ArrayList<Proteina> ArrayPr;
    private Spliceosoma[] splRef;

    //Constantes cineticas
    private double kl; //Ligar
    private double kd; //Desligar
    private double ka; //Avan�ar
    private double kdg; //Degradar prote�na

    /**
     * Construtor da Classe
     */
    public Ribossoma(Spliceosoma[] spl, double kl, double ka, double kd, double kdg){
        velTrad = 15;
        posRib = -1;
        dimRib = 37;
        ligado = false;
        currentM = null;
        currentProt = null;
        ArrayPr = new ArrayList<Proteina>();
        splRef=spl;
        //Constantes cin�ticas
        this.kl = kl;
        this.ka = ka;
        this.kd = kd;
        this.kdg = kdg;
    }

    //Ligar Ribossoma-mRNA
    public void Ligar(mRNA mR){
        if(!ligado){ //N�o estando j� ligado
            mR.LigarRib(this);

            posRib = 0;
            ligado = true;
            currentM = mR;
            currentProt = new Proteina(); //Nova prote�na � criada
            ArrayPr.add(0, currentProt); //Adicionada � lista de prote�nas

			/*kl = 0.0;
			kd = 0.000114;
			ka = 1000;*/

        }
    }

    //Desligar Ribossoma-/-mRNA //Verificar aborto da transcri��o
    public void Desligar(){
        if(ligado){
            if(posRib<currentM.getDim()){ //Aborto da tradu��o
                ArrayPr.remove(currentProt);
                if(currentM.getRibossomasLig().get(0).equals(this)){
                    currentM.setDisp(true);
                }
            }
            currentM.DesligarRib(this);
            posRib = -1;
            ligado = false;
            currentM = null;
            currentProt=null;
			/*kl = 0.4;
			kd = 0.0;
			ka = 0.0;*/
        }
    }

    //Avan�ar
    public void Avan(){
        if(ligado){
            posRib += velTrad; //Incremento de tempo igual � velocidade de transcri��o?

            if(posRib>dimRib){
                currentM.setDisp(true);
            }

            if(posRib>currentM.getDim()){
                currentProt.Fim();
                Desligar();
            }
        }
    }

    //--------------------------------------------Degrada��o--------------------------------------------//
    public void degradarProt(){
        if(!ArrayPr.isEmpty()){
            Proteina prot;
            prot = ArrayPr.get(ArrayPr.size()-1);
            if(prot.isCompleta()){
                ArrayPr.remove(prot);
            }else{
                Desligar();
            }
        }
    }

    //-----------------------------------------Soma - Algor�tmo-----------------------------------------//
    @Override
    public double sum(){
        double sm = getKl() + getKd() + getKa() + getKdg();
        return sm;
    }

    //-----------------------------------------Getters & Setters-----------------------------------------//
    public double getKl() {
        double a=0;
        for(int i=0; i<splRef.length; i++){
            if(!splRef[i].getArrayM().isEmpty()){
                for(int j=0; j<splRef[i].getArrayM().size(); j++){
                    if(splRef[i].getArrayM().get(j).isDisp()){
                        a+=1;
                    }
                }
            }
        }
        return kl*a;
    }

    public double getKa() {

        if(ligado){
            int i = currentM.getRibossomasLig().indexOf(this);
            if(currentM.getRibossomasLig().indexOf(this)==currentM.getRibossomasLig().size()-1 ||
                    posRib+velTrad<currentM.getRibossomasLig().get(i+1).getPosRib()-dimRib){
                return ka;
            }
        }
        return ka*0;
    }

    public double getKd() {
        if(ligado)
            return kd;
        return 0.0;
    }

    public double getKdg() {
        if(ArrayPr.isEmpty())
            return kdg*0;
        return kdg*(ArrayPr.size());
    }

    public int getVelTrad() {
        return velTrad;
    }

    public void setVelTrad(int velTrad) {
        this.velTrad = velTrad;
    }

    public int getPosRib() {
        return posRib;
    }

    public void setPosRib(int posRib) {
        this.posRib = posRib;
    }

    public int getDimRib() {
        return dimRib;
    }

    public void setDimRib(int dimRib) {
        this.dimRib = dimRib;
    }

    public boolean isLigado() {
        return ligado;
    }

    public void setLigado(boolean ligado) {
        this.ligado = ligado;
    }

    public mRNA getCurrent() {
        return currentM;
    }

    public void setCurrent(mRNA current) {
        this.currentM = current;
    }

    public mRNA getCurrentM() {
        return currentM;
    }

    public void setCurrentM(mRNA currentM) {
        this.currentM = currentM;
    }

    public Proteina getCurrentProt() {
        return currentProt;
    }

    public void setCurrentProt(Proteina currentProt) {
        this.currentProt = currentProt;
    }

    public void setKl(double kl) {

        this.kl = kl;
    }

    public void setKd(double kd) {
        this.kd = kd;
    }

    public void setKa(double ka) {
        this.ka = ka;
    }

    public void setKdg(double kdg) {
        this.kdg = kdg;
    }

    public ArrayList<Proteina> getArrayPr() {
        return ArrayPr;
    }

    public void setArrayPr(ArrayList<Proteina> arrayPr) {
        ArrayPr = arrayPr;
    }

    //-----------------------------------------hashCode & equals-----------------------------------------//
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((ArrayPr == null) ? 0 : ArrayPr.hashCode());
        result = prime * result
                + ((currentM == null) ? 0 : currentM.hashCode());
        result = prime * result
                + ((currentProt == null) ? 0 : currentProt.hashCode());
        result = prime * result + dimRib;
        long temp;
        temp = Double.doubleToLongBits(ka);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kd);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kl);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        result = prime * result + (ligado ? 1231 : 1237);
        result = prime * result + posRib;
        result = prime * result + velTrad;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Ribossoma other = (Ribossoma) obj;
        if (ArrayPr == null) {
            if (other.ArrayPr != null)
                return false;
        } else if (!ArrayPr.equals(other.ArrayPr))
            return false;
        if (currentM == null) {
            if (other.currentM != null)
                return false;
        } else if (!currentM.equals(other.currentM))
            return false;
        if (currentProt == null) {
            if (other.currentProt != null)
                return false;
        } else if (!currentProt.equals(other.currentProt))
            return false;
        if (dimRib != other.dimRib)
            return false;
        if (Double.doubleToLongBits(ka) != Double.doubleToLongBits(other.ka))
            return false;
        if (Double.doubleToLongBits(kd) != Double.doubleToLongBits(other.kd))
            return false;
        if (Double.doubleToLongBits(kl) != Double.doubleToLongBits(other.kl))
            return false;
        if (ligado != other.ligado)
            return false;
        if (posRib != other.posRib)
            return false;
        if (velTrad != other.velTrad)
            return false;
        return true;
    }
}
