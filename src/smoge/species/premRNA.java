package smoge.species;

/**
 * Created by fpenim on 12/12/2015.
 */
public class premRNA {

    private int splSitesTT;
    private int splSitesA;
    private int splSitesD;
    private boolean citosol;
    private boolean spliceosoma;
    private Spliceosoma spl;
    private boolean rnaP;

    /**
     * Construtor da Classe
     * @param spliceSitesTT nr de locais de spl total
     * @param splSitesA nr de locais de splice disponiveis
     * @param splSitesD nr de locais de splice efectuados
     */
    public premRNA(int spliceSitesTT,int splSitesA, int splSitesD){
        this.splSitesTT = spliceSitesTT;
        this.splSitesA = splSitesA;
        this.splSitesD = splSitesD;
        citosol = false;
        spliceosoma = false;
        spl = null;
        rnaP = true;
    }

    //Ligar
    public void ligarSpl(Spliceosoma s){ //Junto com o Spliceosoma
        spliceosoma = true;
        spl = s;
    }

    //Splicing
    public boolean splicing(){ //
        if(splSitesD<splSitesTT){
            splSitesD += 1;
            splSitesA -= 1;
        }
        if(splSitesD==splSitesTT){
            return true;
        }
        return false;
    }

    //Transporte
    public void transporte(){ // Junto com Slpiceosoma
        citosol = true;
        spliceosoma = false;
        spl = null;
    }

    //Getters and Setters
    public int getSplSitesTT() {
        return splSitesTT;
    }

    public int getSplSitesA() {
        return splSitesA;
    }

    public void setSplSitesA(int splSitesA) {
        this.splSitesA = splSitesA;
    }

    public int getSplSitesD() {
        return splSitesD;
    }

    public void setSplSitesD(int splSitesD) {
        this.splSitesD = splSitesD;
    }

    public boolean isCitosol() {
        return citosol;
    }

    public void setCitosol(boolean citosol) {
        this.citosol = citosol;
    }

    public boolean isSpliceosoma() {
        return spliceosoma;
    }

    public void setSpliceosoma(boolean spliceosoma) {
        this.spliceosoma = spliceosoma;
    }

    public Spliceosoma getSpl() {
        return spl;
    }

    public void setSpl(Spliceosoma spl) {
        this.spl = spl;
    }

    public boolean isRnaP() {
        return rnaP;
    }

    public void setRnaP(boolean rnaP) {
        this.rnaP = rnaP;
    }

    //----------------------------------------- Override -----------------------------------------//
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((spl == null) ? 0 : spl.hashCode());
        result = prime * result + splSitesA;
        result = prime * result + splSitesD;
        result = prime * result + splSitesTT;
        result = prime * result + (spliceosoma ? 1231 : 1237);
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
        premRNA other = (premRNA) obj;
        if (spl == null) {
            if (other.spl != null)
                return false;
        } else if (!spl.equals(other.spl))
            return false;
        if (splSitesA != other.splSitesA)
            return false;
        if (splSitesD != other.splSitesD)
            return false;
        if (splSitesTT != other.splSitesTT)
            return false;
        if (spliceosoma != other.spliceosoma)
            return false;
        return true;
    }
}
