package smoge.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Spliceosome extends Element {

    private boolean connected;
    private PrecursormRNA currentP;
    private ArrayList<MessengerRNA> arrayM;
    private RNApolymerase[] polRef;

    // Kinetic constants
    private double kc; // connect
    private double ks; // splice
    private double kt; // Transport -> nuc-cit
    private double kdg; // MessengerRNA Degradation

    /**
     * Constructor
     */
    public Spliceosome(RNApolymerase[] rna, double kc, double ks, double kt, double kdg){
        this.kc = kc;
        this.ks = ks;
        this.kt = kt;
        this.kdg = kdg;
        connected = false;
        currentP = null;
        arrayM = new ArrayList<MessengerRNA>();
        polRef = rna;
    }

    //connect
    public boolean connect(PrecursormRNA pmR) {
        if (!connected && !pmR.isSpliceosome() && !pmR.isCitosol()) { //If not connected
            pmR.connectSpliceosome(this);
            connected = true;
            currentP = pmR;
            //kc = 0.0; // N�o se pode connect (ja esta)
            //ks = 125.5; //Pode ocorrer splicing
            return true;
        }
        return false;
    }

    //splice
    public void splice() {
        boolean a = currentP.splicing();
		/*if (a) { //All splicing sites done
			//ks = 0.0; // No more splicing
			//kt = 533; // Transport is possible
		}*/
    }

    // Transport n-c
    public void transport(int dg) {
        if (connected && currentP.getSplSitesTT() == currentP.getSplSitesD() && !currentP.isRnaP()) {
            currentP.transporte();
            disconnect();

            MessengerRNA m = new MessengerRNA(dg);
            arrayM.add(0, m);
        }
    }

    //disconnect
    public void disconnect() {
        currentP = null; // Not connected to any pre-MessengerRNA
        connected = false;
        //kc = 56.2; // May connect
        //ks = 0.0; // Can't do any splicing
        //kt = 0.0; // No transport
    }

    // MessengerRNA degradation
    public void degrademRNA() {
        if (!arrayM.isEmpty()) {
            MessengerRNA m;
            m = arrayM.get(arrayM.size()-1); // eldest MessengerRNA
            if (m.getRibosomes().isEmpty()) { // If no ribosome connected
                arrayM.remove(m);
            } else {
                for (int i = 0; i < m.getRibosomes().size(); i++) { // Goes trough all connected ribosomes
                    m.getRibosomes().get(i).disconnect(); // Disconnects ribosomes
                }
                arrayM.remove(m);
            }
        }
    }

    // Algorithm - sum
    @Override
    public double sum() {
        double sm = getKc()+ getKs() + getKt() + getKdg();
        return sm;
    }

    // Getters & Setters
    public double getKc() {
        int a = 0;
        for (int i = 0; i < polRef.length; i++) {
            if (!polRef[i].getPmRNA().isEmpty()) {
                for (int j=0; j < polRef[i].getPmRNA().size(); j++) {
                    if (!polRef[i].getPmRNA().get(j).isSpliceosome() && !polRef[i].getPmRNA().get(j).isCitosol()) {
                        a += 1;
                    }
                }
            }
        }
        double sm = kc * ((double)a);
        return sm;
    }

    public double getKs() {
        if (connected && currentP.getSplSitesA() > 0) {
            return ks * ((double)currentP.getSplSitesA());
        }
        return ks * 0;
    }

    public double getKt() {
        if (connected && currentP.getSplSitesTT()==currentP.getSplSitesD() && !currentP.isRnaP())
            return kt;
        return 0.0;
    }

    public double getKdg() {
        if(arrayM.isEmpty())
            return kdg * (double)0;
        return kdg *(arrayM.size());
    }

    public void setKc(double kc) {
        this.kc = kc;
    }

    public void setKs(double ks) {
        this.ks = ks;
    }

    public void setKt(double kt) {
        this.kt = kt;
    }

    public void setKdg(double kdg) {
        this.kdg = kdg;
    }

    public boolean isConnected() {
        return connected;
    }

    public void setConnected(boolean connected) {
        this.connected = connected;
    }

    public PrecursormRNA getCurrentP() {
        return currentP;
    }

    public void setCurrentP(PrecursormRNA currentP) {
        this.currentP = currentP;
    }

    public ArrayList<MessengerRNA> getArrayM() {
        return arrayM;
    }

    public void setArrayM(ArrayList<MessengerRNA> arrayM) {
        this.arrayM = arrayM;
    }

    // Override
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = super.hashCode();
        result = prime * result + ((arrayM == null) ? 0 : arrayM.hashCode());
        result = prime * result
                + ((currentP == null) ? 0 : currentP.hashCode());
        long temp;
        temp = Double.doubleToLongBits(kdg);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kc);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(ks);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(kt);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        result = prime * result + (connected ? 1231 : 1237);
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
        Spliceosome other = (Spliceosome) obj;
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
        if (Double.doubleToLongBits(kdg) != Double.doubleToLongBits(other.kdg))
            return false;
        if (Double.doubleToLongBits(kc) != Double.doubleToLongBits(other.kc))
            return false;
        if (Double.doubleToLongBits(ks) != Double.doubleToLongBits(other.ks))
            return false;
        if (Double.doubleToLongBits(kt) != Double.doubleToLongBits(other.kt))
            return false;
        if (connected != other.connected)
            return false;
        return true;
    }
}
