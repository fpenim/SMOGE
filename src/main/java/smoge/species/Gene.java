package main.java.smoge.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Gene {

    private boolean availablePromotor;
    private int geneLength;
    private ArrayList<RNApol> connectedRNApol;

    public Gene (int geneLength) {
        availablePromotor = true; //Qd criado um gene tem sempre o promotor disponivel
        this.geneLength = geneLength; //Nr de pb do gene
        connectedRNApol = new ArrayList<RNApol>(); //Lista onde serao adicionadas todas as RNApol ligadas ao gene
    }

    //Verifica disponibilidade do promotor
	/*public boolean isAvaiable(){
		if (!availablePromotor)//Se o promotor estiver ocupado, devolve false
			return false;
		if (connectedRNApol.isEmpty())//Se a lista estiver vazia e o promotor disponivel, devolve true
			return true;
		//Se a lista tiver pelo menos um elemento
		//Verifica se a ultima RNApol adicionada ja avancou o suf
		if(connectedRNApol.get(0).getRNApolPosition()-((connectedRNApol.get(0).getRNApolDimension()/2)+1)>
			(connectedRNApol.get(0).getRNApolDimension()/2)){
				return true;
		}
		return false;
	}*/

    /**
     *
     * @param rna
     * @requires this.isAvaiable()==true
     */
    public void LigarRNA (RNApol rna) {
        connectedRNApol.add(0, rna); //Adiciona a rna na 1a posicao
        availablePromotor =false; //Passa promotor a ocupado
    }

    //disconnect
    public void Disconnect (RNApol rna) {
        if (connectedRNApol.get(0).equals(rna)) {
            this.availablePromotor = true;
        }
        connectedRNApol.remove(rna);
    }

    // Setters & Getters
    public void setAvailablePromotor(boolean availablePromotor){
        this.availablePromotor = availablePromotor;
    }

    public boolean getAvailablePromotor(){
        return availablePromotor;
    }

    public int getGeneLength(){
        return geneLength;
    }

    public ArrayList<RNApol> getConnectedRNApol() {
        return connectedRNApol;
    }

    public void setConnectedRNApol(ArrayList<RNApol> rNApolLigados) {
        connectedRNApol = rNApolLigados;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((connectedRNApol == null) ? 0 : connectedRNApol.hashCode());
        result = prime * result + geneLength;
        result = prime * result + (availablePromotor ? 1231 : 1237);
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
        Gene other = (Gene) obj;
        if (connectedRNApol == null) {
            if (other.connectedRNApol != null)
                return false;
        } else if (!connectedRNApol.equals(other.connectedRNApol))
            return false;
        if (geneLength != other.geneLength)
            return false;
        if (availablePromotor != other.availablePromotor)
            return false;
        return true;
    }
}
