package smoge.species;

import java.util.ArrayList;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Gene {

    private boolean promotorDisp;
    private int dimGene;
    private ArrayList<RNApol> RNApolLigados;

    public Gene(int dimGene){
        promotorDisp = true;//Qd criado um gene tem sempre o promotor disponivel
        this.dimGene = dimGene;//Nr de pb do gene
        RNApolLigados = new ArrayList<RNApol>();//Lista onde serao adicionadas todas as RNApol ligadas ao gene
    }

    //Verifica disponibilidade do promotor
	/*public boolean isAvaiable(){
		if (!promotorDisp)//Se o promotor estiver ocupado, devolve false
			return false;
		if (RNApolLigados.isEmpty())//Se a lista estiver vazia e o promotor disponivel, devolve true
			return true;
		//Se a lista tiver pelo menos um elemento
		//Verifica se a ultima RNApol adicionada ja avancou o suf
		if(RNApolLigados.get(0).getPosRNApol()-((RNApolLigados.get(0).getDimRNApol()/2)+1)>
			(RNApolLigados.get(0).getDimRNApol()/2)){
				return true;
		}
		return false;
	}*/

    /**
     *
     * @param rna
     * @requires this.isAvaiable()==true
     */
    public void LigarRNA(RNApol rna){//Ligar
        RNApolLigados.add(0, rna);//Adiciona a rna na 1a posicao
        promotorDisp=false;//Passa promotor a ocupado
    }

    //Desligar
    public void DesligarRNA(RNApol rna){
        if(RNApolLigados.get(0).equals(rna)){
            this.promotorDisp = true;
        }
        RNApolLigados.remove(rna);
    }

    //Setters and Getters
    public void setPromotorDisp(boolean promotorDisp){
        this.promotorDisp=promotorDisp;
    }

    public boolean getPromotorDisp(){
        return promotorDisp;
    }

    public int getDimGene(){
        return dimGene;
    }

    public ArrayList<RNApol> getRNApolLigados() {
        return RNApolLigados;
    }

    public void setRNApolLigados(ArrayList<RNApol> rNApolLigados) {
        RNApolLigados = rNApolLigados;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((RNApolLigados == null) ? 0 : RNApolLigados.hashCode());
        result = prime * result + dimGene;
        result = prime * result + (promotorDisp ? 1231 : 1237);
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
        if (RNApolLigados == null) {
            if (other.RNApolLigados != null)
                return false;
        } else if (!RNApolLigados.equals(other.RNApolLigados))
            return false;
        if (dimGene != other.dimGene)
            return false;
        if (promotorDisp != other.promotorDisp)
            return false;
        return true;
    }
}
