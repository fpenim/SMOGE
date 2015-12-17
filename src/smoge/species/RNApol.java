package smoge.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class RNApol extends Element {
    private int transcriptionVelocity;
    private int RNApolPosition;
    private int RNApolDimension;
    private boolean connected = false;
    private Gene gene;
    private Gene currentGene;
    private premRNA currentP;
    private ArrayList<premRNA> pmRNA;

    // Kinetic constants
    private double kc; // connect
    private double kd; // disconnect
    private double kp; // progress
    private double kdg; // pre-mRNA degradation

    /**
     * Constructor
     * @param gene  Gene
     * @param kc  k connect
     * @param kp  k avancar
     * @param kd  k desligar
     * @param kdg k degradar pre-mRNA
     */
    public RNApol (Gene gene, double kc, double kp, double kd, double kdg) {
        transcriptionVelocity = 39;
        RNApolDimension = 25;
        RNApolPosition = -1;
        connected = false;
        currentGene = null;
        currentP = null;
        pmRNA = new ArrayList<premRNA>();
        this.gene = gene;

        this.kc = kc;
        this.kp = kp;
        this.kd = kd;
        this.kdg = kdg;
    }

    // Connect
    public boolean connect (Gene gen, int sp) { //Recebe o gene e a lista geral dos pre-mRNA
        this.clean();

        gen.LigarRNA(this);
        currentGene = gen; //Associa o gene a que a RNApol esta connected
        RNApolPosition = 0; //Passa a posicao a 0(promotor)
        connected = true; //connected passa a verdadeiro

        //kc = 0.0;
        //kd = 0.00019;
        //kp = 65;

        premRNA p = new premRNA(sp,0,0); //� criado um novo pr�-mRNA
        currentP = p;
        pmRNA.add(0, p); //Este � adicionado � lista de pre-mRNA's

        return true;
    }

    // Disconnect
    public void Disconnect () {
        this.clean();
        if (connected) {
            currentP.setRnaP(false);

            if (RNApolPosition < currentGene.getGeneLength()) { // Check transcription abortion
                if (currentP.isSpliceosome()) { // Any connected spliceosome to pre-mRNA
                    currentP.getSpl().disconnect(); //disconnect spliceosome
                }
                pmRNA.remove(currentP); //Eliminar o pre-mRNA incompleto
                if (currentGene.getConnectedRNApol().get(0).equals(this)) { //se for o 1� RNApol da lista (ultimo a connect-se)
                    currentGene.setAvailablePromotor(true);
                }
            }

            currentGene.setAvailablePromotor(true);
            currentGene.Disconnect(this); //disconnect gene

            connected = false;
            currentGene = null;
            currentP = null;
            RNApolPosition = -1;

            //kc = 0.0245;
            //kd = 0.0;
            //kp = 0.0;
        }

    }

    // Progress
    public void progress (int s) {
        RNApolPosition += transcriptionVelocity;

        int gsize = currentGene.getGeneLength();
        int lim = gsize/(s+1);

        if (RNApolPosition > lim && currentP.getSplSitesA() + currentP.getSplSitesD() < 1 && s >= 1) {
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }
        else if (RNApolPosition > (lim*2) && currentP.getSplSitesA() + currentP.getSplSitesD() < 2 && s >= 2) {
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }
        else if (RNApolPosition > (lim*3) && currentP.getSplSitesA() + currentP.getSplSitesD() < 3 && s >= 3) {
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }
        else if (RNApolPosition > (lim*4) && currentP.getSplSitesA() + currentP.getSplSitesD() < 4 && s >= 4) {
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }
        else if (RNApolPosition > (lim*5) && currentP.getSplSitesA() + currentP.getSplSitesD() < 5 && s >= 5) {
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }
        else if (RNApolPosition > (lim*6) && currentP.getSplSitesA() + currentP.getSplSitesD() < 6 && s >= 6) {
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }
        else if (RNApolPosition > (lim*7) && currentP.getSplSitesA() + currentP.getSplSitesD() < 7 && s >= 7) {
            currentP.setSplSitesA(currentP.getSplSitesA()+1);
        }

		/*Locais de splice
		if(RNApolPosition>350 && currentP.getSplSitesA()+currentP.getSplSitesD()<1){
			currentP.setSplSitesA(currentP.getSplSitesA()+1);
		}else
			if(RNApolPosition>650 && currentP.getSplSitesA()+currentP.getSplSitesD()<2){
				currentP.setSplSitesA(currentP.getSplSitesA()+1);
			}else
				if(RNApolPosition>950 && currentP.getSplSitesA()+currentP.getSplSitesD()<3){
					currentP.setSplSitesA(currentP.getSplSitesA()+1);
				}else
					if(RNApolPosition>1250 && currentP.getSplSitesA()+currentP.getSplSitesD()<4){
						currentP.setSplSitesA(currentP.getSplSitesA()+1);
					}else
						if(RNApolPosition>1550 && currentP.getSplSitesA()+currentP.getSplSitesD()<5){
							currentP.setSplSitesA(currentP.getSplSitesA()+1);
						}else
							if(RNApolPosition>1850 && currentP.getSplSitesA()+currentP.getSplSitesD()<6){
								currentP.setSplSitesA(currentP.getSplSitesA()+1);
						}else
							if(RNApolPosition>2150 && currentP.getSplSitesA()+currentP.getSplSitesD()<7){
								currentP.setSplSitesA(currentP.getSplSitesA()+1);
						}*/

        //Quando RNApol avanca o suficiente
        if (RNApolPosition > RNApolDimension) {
            currentGene.setAvailablePromotor(true);
        }
        //Qd RNApol chega ao final do gene
        if (RNApolPosition >= currentGene.getGeneLength()) {
            Disconnect();
        }
    }

    //Eliminar pre-mRNA's j� processados em mRNA's
    public void clean () {
        if (!pmRNA.isEmpty()) {
            for (int i=0; i<pmRNA.size(); i++) {
                if (pmRNA.get(i).isCitosol()) {
                    pmRNA.remove(i);
                }
            }
        }
    }

    // Degradation
    public void degradePmRNA () {
        this.clean();
        if (!pmRNA.isEmpty()) {
            premRNA pm;
            pm = pmRNA.get(pmRNA.size()-1);
            if (pm.isRnaP() && pm.isSpliceosome()) { //Se estiver ligado a um spliceosoma e a uma RNApol
                pm.getSpl().disconnect();
                this.Disconnect();
            } else if (pm.isRnaP() && !pm.isSpliceosome()) { //Se estiver ligado a uma RNApol
                this.Disconnect();
            } else if (!pm.isRnaP() && pm.isSpliceosome()) { //Se estiver ligado a um spliceosoma
                pm.getSpl().disconnect();
            }
            pmRNA.remove(pm);
        }
    }


    // Algorithm - sum
    @Override
    public double sum () {

        double sm = getKc() + getKp() + getKd() + getKdg();
        return sm;
    }

    // Setters & Getters
    public double getKc () {
        if( !gene.getAvailablePromotor()){
            return kc *0;
        } else {
            return kc;
        }
    }

    public double getKp() {
        if (connected) {
            int lim = currentGene.getConnectedRNApol().indexOf(this); //Posi��o do elemento na lista de RNApol ligados
            if (currentGene.getConnectedRNApol().indexOf(this)== currentGene.getConnectedRNApol().size()-1 ||
                    RNApolPosition <(currentGene.getConnectedRNApol().get(lim+1).RNApolPosition)- RNApolDimension) {
                return kp;
            }
        }
        return kp * 0;
    }

    public double getKd () {
        if (connected)
            return kd;
        return 0.0;
    }

    public double getKdg () {
        if (pmRNA.isEmpty())
            return kdg * 0;
        return kdg * (pmRNA.size());
    }

    public int getRNApolDimension() {
        return RNApolDimension;
    }

    public ArrayList<premRNA> getPmRNA() {
        return pmRNA;
    }

    public void setPmRNA(ArrayList<premRNA> pmRNA) {
        this.pmRNA = pmRNA;
    }

    public int getTranscriptionVelocity() {
        return transcriptionVelocity;
    }

    public void setTranscriptionVelocity(int transcriptionVelocity) {
        this.transcriptionVelocity = transcriptionVelocity;
    }

    public int getRNApolPosition() {
        return RNApolPosition;
    }

    public void setRNApolPosition(int RNApolPosition) {
        this.RNApolPosition = RNApolPosition;
    }

    public boolean isConnected() {
        return connected;
    }

    public void setConnected(boolean connected) {
        this.connected = connected;
    }

    public Gene getCurrentGene(){
        return currentGene;
    }

    public void setCurrent(Gene gen) {
        this.currentGene = gen;
    }

    public premRNA getCurrentP() {
        return currentP;
    }

    public void setCurrentP(premRNA currentP) {
        this.currentP = currentP;
    }

    // Override
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((currentGene == null) ? 0 : currentGene.hashCode());
        result = prime * result
                + ((currentP == null) ? 0 : currentP.hashCode());
        result = prime * result + RNApolDimension;
        long temp;
        temp = Double.doubleToLongBits(kp);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kd);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kc);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        result = prime * result + (connected ? 1231 : 1237);
        result = prime * result + ((pmRNA == null) ? 0 : pmRNA.hashCode());
        result = prime * result + RNApolPosition;
        result = prime * result + transcriptionVelocity;
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
        if (currentGene == null) {
            if (other.currentGene != null)
                return false;
        } else if (!currentGene.equals(other.currentGene))
            return false;
        if (currentP == null) {
            if (other.currentP != null)
                return false;
        } else if (!currentP.equals(other.currentP))
            return false;
        if (RNApolDimension != other.RNApolDimension)
            return false;
        if (Double.doubleToLongBits(kp) != Double.doubleToLongBits(other.kp))
            return false;
        if (Double.doubleToLongBits(kd) != Double.doubleToLongBits(other.kd))
            return false;
        if (Double.doubleToLongBits(kc) != Double.doubleToLongBits(other.kc))
            return false;
        if (connected != other.connected)
            return false;
        if (pmRNA == null) {
            if (other.pmRNA != null)
                return false;
        } else if (!pmRNA.equals(other.pmRNA))
            return false;
        if (RNApolPosition != other.RNApolPosition)
            return false;
        if (transcriptionVelocity != other.transcriptionVelocity)
            return false;
        return true;
    }


}
