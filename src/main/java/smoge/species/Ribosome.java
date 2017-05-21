package smoge.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Ribosome extends Element {

    private int translationVelocity;
    private int ribosomePosition;
    private int ribosomeDimension;
    private boolean connected = false;
    private MessengerRNA currentM;
    private Protein currentProt;
    private ArrayList<Protein> protArray;
    private Spliceosome[] splRef;

    // Kinetic constants
    private double kc; // connect
    private double kd; // disconnect
    private double kp; // progress
    private double kdg; // Protein Degradation

    /**
     * Constructor
     */
    public Ribosome (Spliceosome[] spl, double kc, double kp, double kd, double kdg) {
        translationVelocity = 15;
        ribosomePosition = -1;
        ribosomeDimension = 37;
        connected = false;
        currentM = null;
        currentProt = null;
        protArray = new ArrayList<Protein>();
        splRef=spl;
        // Kinetic constants
        this.kc = kc;
        this.kp = kp;
        this.kd = kd;
        this.kdg = kdg;
    }

    // connect Ribosome-MessengerRNA
    public void connect(MessengerRNA mRNA) {
        if (!connected){ // If is no longer connected
            mRNA.connectRib(this);

            ribosomePosition = 0;
            connected = true;
            currentM = mRNA;
            currentProt = new Protein(); // New protein
            protArray.add(0, currentProt); // Add to protein list
			/*kc = 0.0;
			kd = 0.000114;
			kp = 1000;*/

        }
    }

    // disconnect Ribosome-/-MessengerRNA
    // Check transcription abortion
    public void disconnect() {
        if (connected){
            if (ribosomePosition < currentM.getLength()) { // Transcription abortion
                protArray.remove(currentProt);
                if (currentM.getRibosomes().get(0).equals(this)) {
                    currentM.setAvailable(true);
                }
            }
            currentM.disconnectRib(this);
            ribosomePosition = -1;
            connected = false;
            currentM = null;
            currentProt=null;
			/*kc = 0.4;
			kd = 0.0;
			kp = 0.0;*/
        }
    }

    // progress
    public void progress(){
        if(connected){
            ribosomePosition += translationVelocity; // Time step equal to transcription velocity (?)

            if (ribosomePosition > ribosomeDimension) {
                currentM.setAvailable(true);
            }

            if (ribosomePosition >currentM.getLength()) {
                currentProt.finished();
                disconnect();
            }
        }
    }

    // Protein Degradation
    public void degradeProtein() {
        if (!protArray.isEmpty()) {
            Protein prot;
            prot = protArray.get(protArray.size()-1);
            if (prot.isComplete()) {
                protArray.remove(prot);
            } else {
                disconnect();
            }
        }
    }

    // Algorithm - sum
    @Override
    public double sum() {
        double sm = getKc() + getKd() + getKp() + getKdg();
        return sm;
    }

    // Getters & Setters
    public double getKc() {
        double a = 0;
        for (int i = 0; i < splRef.length; i++) {
            if (!splRef[i].getArrayM().isEmpty()) {
                for (int j = 0; j < splRef[i].getArrayM().size(); j++) {
                    if (splRef[i].getArrayM().get(j).isAvailable()) {
                        a += 1;
                    }
                }
            }
        }
        return kc * a;
    }

    public double getKp() {
        if (connected) {
            int i = currentM.getRibosomes().indexOf(this);
            if (currentM.getRibosomes().indexOf(this) == currentM.getRibosomes().size()-1 ||
                    ribosomePosition + translationVelocity < currentM.getRibosomes().get(i+1).getRibosomePosition() - ribosomeDimension) {
                return kp;
            }
        }
        return kp * 0;
    }

    public double getKd() {
        if (connected)
            return kd;
        return 0.0;
    }

    public double getKdg() {
        if (protArray.isEmpty())
            return kdg*0;
        return kdg * (protArray.size());
    }

    public int getTranslationVelocity() {
        return translationVelocity;
    }

    public void setTranslationVelocity(int translationVelocity) {
        this.translationVelocity = translationVelocity;
    }

    public int getRibosomePosition() {
        return ribosomePosition;
    }

    public void setRibosomePosition(int ribosomePosition) {
        this.ribosomePosition = ribosomePosition;
    }

    public int getRibosomeDimension() {
        return ribosomeDimension;
    }

    public void setRibosomeDimension(int ribosomeDimension) {
        this.ribosomeDimension = ribosomeDimension;
    }

    public boolean isConnected() {
        return connected;
    }

    public void setConnected(boolean connected) {
        this.connected = connected;
    }

    public MessengerRNA getCurrent() {
        return currentM;
    }

    public void setCurrent(MessengerRNA current) {
        this.currentM = current;
    }

    public MessengerRNA getCurrentM() {
        return currentM;
    }

    public void setCurrentM(MessengerRNA currentM) {
        this.currentM = currentM;
    }

    public Protein getCurrentProt() {
        return currentProt;
    }

    public void setCurrentProt(Protein currentProt) {
        this.currentProt = currentProt;
    }

    public void setKc(double kc) {

        this.kc = kc;
    }

    public void setKd(double kd) {
        this.kd = kd;
    }

    public void setKp(double kp) {
        this.kp = kp;
    }

    public void setKdg(double kdg) {
        this.kdg = kdg;
    }

    public ArrayList<Protein> getProtArray() {
        return protArray;
    }

    public void setProtArray(ArrayList<Protein> protArray) {
        this.protArray = protArray;
    }

    // Override
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((protArray == null) ? 0 : protArray.hashCode());
        result = prime * result
                + ((currentM == null) ? 0 : currentM.hashCode());
        result = prime * result
                + ((currentProt == null) ? 0 : currentProt.hashCode());
        result = prime * result + ribosomeDimension;
        long temp;
        temp = Double.doubleToLongBits(kp);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kd);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kc);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        result = prime * result + (connected ? 1231 : 1237);
        result = prime * result + ribosomePosition;
        result = prime * result + translationVelocity;
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
        Ribosome other = (Ribosome) obj;
        if (protArray == null) {
            if (other.protArray != null)
                return false;
        } else if (!protArray.equals(other.protArray))
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
        if (ribosomeDimension != other.ribosomeDimension)
            return false;
        if (Double.doubleToLongBits(kp) != Double.doubleToLongBits(other.kp))
            return false;
        if (Double.doubleToLongBits(kd) != Double.doubleToLongBits(other.kd))
            return false;
        if (Double.doubleToLongBits(kc) != Double.doubleToLongBits(other.kc))
            return false;
        if (connected != other.connected)
            return false;
        if (ribosomePosition != other.ribosomePosition)
            return false;
        if (translationVelocity != other.translationVelocity)
            return false;
        return true;
    }
}
