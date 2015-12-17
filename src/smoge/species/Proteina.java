package smoge.species;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Proteina {

    private boolean completa;

    /**
     * Construtor da Classe
     */
    public Proteina(){
        completa = false;
    }

    //Final da Tradução
    public boolean Fim(){
        completa = true;
        return completa;
    }

    //---------------------------------------Getters and Setters---------------------------------------//
    public boolean isCompleta() {
        return completa;
    }

    public void setCompleta(boolean completa) {
        this.completa = completa;
    }

    //------------------------------------------- Override -------------------------------------------//
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + (completa ? 1231 : 1237);
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
        Proteina other = (Proteina) obj;
        if (completa != other.completa)
            return false;
        return true;
    }

}
