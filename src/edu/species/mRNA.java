package edu.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class mRNA {
    private int dim; // Dimens�o do mRNA
    private boolean disp;
    private ArrayList<Ribossoma> RibossomasLig; //Lista de ribossomas ligados ao mRNA

    /**
     * Contrutor da Classe
     */
    public mRNA(int dg){
        dim = dg/3;
        RibossomasLig = new ArrayList<Ribossoma>();
        disp=true;
    }

    //----------------------------------------- A��es do mRNA -----------------------------------------//

    //Ligar mRNA-Ribossoma
    public void LigarRib(Ribossoma rib){
        RibossomasLig.add(0,rib);
        disp=false;
        //return true;

    }

    // Desigar mRNA-/-Ribossoma
    public void DesligarRib(Ribossoma rib){ //Recebe o ribossoma a desligar
        RibossomasLig.remove(rib);
    }

    //-----------------------------------------Getters and Setters-----------------------------------------//
    public int getDim() {
        return dim;
    }

    public void setDim(int dim) {
        this.dim = dim;
    }

    public ArrayList<Ribossoma> getRibossomasLig() {
        return RibossomasLig;
    }

    public void setRibossomasLig(ArrayList<Ribossoma> ribossomasLig) {
        RibossomasLig = ribossomasLig;
    }

    public boolean isDisp() {
        return disp;
    }

    public void setDisp(boolean disp) {
        this.disp = disp;
    }

    //----------------------------------------- Override -----------------------------------------//
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((RibossomasLig == null) ? 0 : RibossomasLig.hashCode());
        result = prime * result + dim;
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
        if (RibossomasLig == null) {
            if (other.RibossomasLig != null)
                return false;
        } else if (!RibossomasLig.equals(other.RibossomasLig))
            return false;
        if (dim != other.dim)
            return false;
        return true;
    }
}
