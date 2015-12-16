package edu.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class RNApol extends Element {
    private int velTrans;
    private int posRNApol;
    private int dimRNApol;
    private boolean ligada = false;
    private Gene g;
    private Gene currentG;
    private premRNA currentP;
    private ArrayList<premRNA> pmRNA;

    //Constantes cineticas
    private double kl; //Ligar
    private double kd; //Desligar
    private double ka; //Avan�ar
    private double kdg; //Degradar pre-mRNA's

    /**
     * Construtor da Classe
     * @param gt  Gene
     * @param kl  k ligar
     * @param ka  k avancar
     * @param kd  k desligar
     * @param kdg k degradar pre-mRNA
     */
    public RNApol(Gene gt, double kl, double ka, double kd, double kdg){
        velTrans = 39; //Velocidade de transcricao
        dimRNApol = 25; //Nr de nuc ocupados por uma RNApol
        posRNApol =-1; //Valor a dar qd se cria a RNApol e nao esta ainda ligada
        ligada = false; //Ligada a um gene
        currentG = null; //Gene a que esta ligada
        currentP = null; //pre-mRNA em s�ntese
        pmRNA = new ArrayList<premRNA>(); //pre-mRNA em sintese/sintetizados
        g = gt;

        this.kl = kl;
        this.ka = ka;
        this.kd = kd;
        this.kdg = kdg;
    }
    //---------------------------------------------- A��es ----------------------------------------//
    //Ligar
    public boolean Ligar(Gene gen, int sp) { //Recebe o gene e a lista geral dos pre-mRNA
        this.limpar();

        gen.LigarRNA(this);
        currentG = gen; //Associa o gene a que a RNApol esta ligada
        posRNApol = 0; //Passa a posicao a 0(promotor)
        ligada = true; //ligada passa a verdadeiro

        //kl = 0.0;
        //kd = 0.00019;
        //ka = 65;

        premRNA p = new premRNA(sp,0,0); //� criado um novo pr�-mRNA
        currentP = p;
        pmRNA.add(0, p); //Este � adicionado � lista de pre-mRNA's

        return true;
    }

    //Desligar
    public void Desligar() {
        this.limpar();
        if(ligada){
            currentP.setRnaP(false);
            //Verificar aborto da transcri��o
            if(posRNApol<currentG.getDimGene()){ //Se a RNApol n�o tiver terminado a transcri��o
                if(currentP.isSpliceosoma()){ //Se j� estiver algum spliceosoma ligado ao pre-mRNA
                    currentP.getSpl().desligarS();//Desligar spliceosoma
                }
                pmRNA.remove(currentP); //Eliminar o pre-mRNA incompleto
                if(currentG.getRNApolLigados().get(0).equals(this)){ //se for o 1� RNApol da lista (ultimo a ligar-se)
                    currentG.setPromotorDisp(true);
                }
            }
            currentG.setPromotorDisp(true);
            currentG.DesligarRNA(this); //Desligar do gene

            ligada=false; //ligada passa a falso
            currentG=null; //Nao esta associada a nenhum gene
            currentP=null; //Nao esta a sintetizar nenhum pre-mRNA
            posRNApol=-1; //Passa a posicao a -1(desligada)

            //kl = 0.0245;
            //kd = 0.0;
            //ka = 0.0;
        }

    }

    //Avan�ar (sem ultrapassar)
    public void Avan(int s){
        //Incremento da posicao
        posRNApol += velTrans;

        int gsize =currentG.getDimGene();
        int lim = gsize/(s+1);

        if(posRNApol>lim && currentP.getSplSitesA()+currentP.getSplSitesD()<1 && s>=1){
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }else
        if(posRNApol>(lim*2) && currentP.getSplSitesA()+currentP.getSplSitesD()<2 && s>=2){
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }else
        if(posRNApol>(lim*3) && currentP.getSplSitesA()+currentP.getSplSitesD()<3 && s>=3){
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }else
        if(posRNApol>(lim*4) && currentP.getSplSitesA()+currentP.getSplSitesD()<4 && s>=4){
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }else
        if(posRNApol>(lim*5) && currentP.getSplSitesA()+currentP.getSplSitesD()<5 && s>=5){
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }else
        if(posRNApol>(lim*6) && currentP.getSplSitesA()+currentP.getSplSitesD()<6 && s>=6){
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }else
        if(posRNApol>(lim*7) && currentP.getSplSitesA()+currentP.getSplSitesD()<7 && s>=7){
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }

		/*Locais de Splicing
		if(posRNApol>350 && currentP.getSplSitesA()+currentP.getSplSitesD()<1){
			currentP.setSplSitesA(currentP.getSplSitesA()+1);
		}else
			if(posRNApol>650 && currentP.getSplSitesA()+currentP.getSplSitesD()<2){
				currentP.setSplSitesA(currentP.getSplSitesA()+1);
			}else
				if(posRNApol>950 && currentP.getSplSitesA()+currentP.getSplSitesD()<3){
					currentP.setSplSitesA(currentP.getSplSitesA()+1);
				}else
					if(posRNApol>1250 && currentP.getSplSitesA()+currentP.getSplSitesD()<4){
						currentP.setSplSitesA(currentP.getSplSitesA()+1);
					}else
						if(posRNApol>1550 && currentP.getSplSitesA()+currentP.getSplSitesD()<5){
							currentP.setSplSitesA(currentP.getSplSitesA()+1);
						}else
							if(posRNApol>1850 && currentP.getSplSitesA()+currentP.getSplSitesD()<6){
								currentP.setSplSitesA(currentP.getSplSitesA()+1);
						}else
							if(posRNApol>2150 && currentP.getSplSitesA()+currentP.getSplSitesD()<7){
								currentP.setSplSitesA(currentP.getSplSitesA()+1);
						}*/

        //Quando RNApol avanca o suficiente
        if(posRNApol>dimRNApol){
            currentG.setPromotorDisp(true);
        }
        //Qd RNApol chega ao final do gene
        if(posRNApol>=currentG.getDimGene()){
            Desligar();
        }
    }

    //Eliminar pre-mRNA's j� processados em mRNA's
    public void limpar(){
        if(!pmRNA.isEmpty()){
            for(int i=0; i<pmRNA.size(); i++){
                if(pmRNA.get(i).isCitosol()){
                    pmRNA.remove(i);
                }
            }
        }
    }

    //--------------------------------------------Degrada��o--------------------------------------------//
    public void degradarPmRNA(){
        this.limpar();
        if(!pmRNA.isEmpty()){
            premRNA pm;
            pm = pmRNA.get(pmRNA.size()-1);
            if(pm.isRnaP() && pm.isSpliceosoma()){ //Se estiver ligado a um spliceosoma e a uma RNApol
                pm.getSpl().desligarS();
                this.Desligar();
            }else if(pm.isRnaP() && !pm.isSpliceosoma()){ //Se estiver ligado a uma RNApol
                this.Desligar();
            }else if(!pm.isRnaP() && pm.isSpliceosoma()){ //Se estiver ligado a um spliceosoma
                pm.getSpl().desligarS();
            }
            pmRNA.remove(pm);
        }
    }

    //---------------------------------------------- Soma ----------------------------------------//
    //Soma - Algor�tmo
    @Override
    public double sum(){

        double sm = getKl() + getKa() + getKd() + getKdg();
        return sm;
    }
    //-------------------------------------- Setters and Getters --------------------------------//
    public double getKl() {
        if(!g.getPromotorDisp()){
            return kl*0;
        }else{
            return kl;
        }
    }

    public double getKa() {
        if(ligada){
            int lim = currentG.getRNApolLigados().indexOf(this); //Posi��o do elemento na lista de RNApol ligados
            if(currentG.getRNApolLigados().indexOf(this)==currentG.getRNApolLigados().size()-1 ||
                    posRNApol<(currentG.getRNApolLigados().get(lim+1).posRNApol)-dimRNApol){
                return ka;
            }
        }
        return ka*0;
    }

    public double getKd(){
        if(ligada)
            return kd;
        return 0.0;
    }

    public double getKdg() {
        if(pmRNA.isEmpty())
            return kdg*0;
        return kdg*(pmRNA.size());
    }

    public int getDimRNApol() {
        return dimRNApol;
    }

    public ArrayList<premRNA> getPmRNA() {
        return pmRNA;
    }

    public void setPmRNA(ArrayList<premRNA> pmRNA) {
        this.pmRNA = pmRNA;
    }

    public int getVelTrans() {
        return velTrans;
    }

    public void setVelTrans(int velTrans) {
        this.velTrans = velTrans;
    }

    public int getPosRNApol() {
        return posRNApol;
    }

    public void setPosRNApol(int posRNApol) {
        this.posRNApol = posRNApol;
    }

    public boolean isLigada() {
        return ligada;
    }

    public void setLigada(boolean ligada) {
        this.ligada = ligada;
    }

    public Gene getCurrentG(){
        return currentG;
    }

    public void setCurrent(Gene gen) {
        this.currentG = gen;
    }

    public premRNA getCurrentP() {
        return currentP;
    }

    public void setCurrentP(premRNA currentP) {
        this.currentP = currentP;
    }

    //---------------------------------------------- Override ----------------------------------------//
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((currentG == null) ? 0 : currentG.hashCode());
        result = prime * result
                + ((currentP == null) ? 0 : currentP.hashCode());
        result = prime * result + dimRNApol;
        long temp;
        temp = Double.doubleToLongBits(ka);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kd);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kl);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        result = prime * result + (ligada ? 1231 : 1237);
        result = prime * result + ((pmRNA == null) ? 0 : pmRNA.hashCode());
        result = prime * result + posRNApol;
        result = prime * result + velTrans;
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
        RNApol other = (RNApol) obj;
        if (currentG == null) {
            if (other.currentG != null)
                return false;
        } else if (!currentG.equals(other.currentG))
            return false;
        if (currentP == null) {
            if (other.currentP != null)
                return false;
        } else if (!currentP.equals(other.currentP))
            return false;
        if (dimRNApol != other.dimRNApol)
            return false;
        if (Double.doubleToLongBits(ka) != Double.doubleToLongBits(other.ka))
            return false;
        if (Double.doubleToLongBits(kd) != Double.doubleToLongBits(other.kd))
            return false;
        if (Double.doubleToLongBits(kl) != Double.doubleToLongBits(other.kl))
            return false;
        if (ligada != other.ligada)
            return false;
        if (pmRNA == null) {
            if (other.pmRNA != null)
                return false;
        } else if (!pmRNA.equals(other.pmRNA))
            return false;
        if (posRNApol != other.posRNApol)
            return false;
        if (velTrans != other.velTrans)
            return false;
        return true;
    }


}
