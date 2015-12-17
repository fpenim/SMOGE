package smoge.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class mRNA {
    private int length;
    private boolean available;
    private ArrayList<Ribosome> ribosomes;

    /**
     * Constructor
     */
    public mRNA (int dg) {
        length = dg / 3;
        ribosomes = new ArrayList<Ribosome>();
        available =true;
    }

    // connect mRNA-Ribosome
    public void connectRib (Ribosome rib) {
        ribosomes.add(0,rib);
        available = false;
    }

    // disconnect mRNA-/-Ribosome
    public void disconnectRib (Ribosome rib){ //Recebe o ribossoma a desligar
        ribosomes.remove(rib);
    }

    // Getters & Setters
    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public ArrayList<Ribosome> getRibosomes() {
        return ribosomes;
    }

    public void setRibosomes(ArrayList<Ribosome> ribosomes) {
        this.ribosomes = ribosomes;
    }

    public boolean isAvailable() {
        return available;
    }

    public void setAvailable(boolean available) {
        this.available = available;
    }

    // Override
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((ribosomes == null) ? 0 : ribosomes.hashCode());
        result = prime * result + length;
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
        mRNA other = (mRNA) obj;
        if (ribosomes == null) {
            if (other.ribosomes != null)
                return false;
        } else if (!ribosomes.equals(other.ribosomes))
            return false;
        if (length != other.length)
            return false;
        return true;
    }
}
