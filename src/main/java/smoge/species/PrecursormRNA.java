package smoge.species;

/**
 * Created by fpenim on 12/12/2015.
 */
public class PrecursormRNA {

    private int splSitesTT;
    private int splSitesA;
    private int splSitesD;
    private boolean citosol;
    private boolean spliceosome;
    private Spliceosome spl;
    private boolean rnaP;

    /**
     * Construtor da Classe
     * @param spliceSitesTT total number of splicing sites
     * @param splSitesA available splicing sites
     * @param splSitesD done splicing sites
     */
    public PrecursormRNA(int spliceSitesTT, int splSitesA, int splSitesD) {
        this.splSitesTT = spliceSitesTT;
        this.splSitesA = splSitesA;
        this.splSitesD = splSitesD;
        citosol = false;
        spliceosome = false;
        spl = null;
        rnaP = true;
    }

    // Connect
    public void connectSpliceosome (Spliceosome s) { // Together with the Slpiceosome
        citosol = true;
        spliceosome = true;
        spl = s;
    }

    //splice
    public boolean splicing () {
        if (splSitesD < splSitesTT) {
            splSitesD += 1;
            splSitesA -= 1;
        }
        if (splSitesD == splSitesTT) {
            return true;
        }
        return false;
    }

    // Transport
    public void transporte () { // Together with the Slpiceosome
        citosol = true;
        spliceosome = false;
        spl = null;
    }

    // Getters & Setters
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

    public boolean isSpliceosome() {
        return spliceosome;
    }

    public void setSpliceosome(boolean spliceosome) {
        this.spliceosome = spliceosome;
    }

    public Spliceosome getSpl() {
        return spl;
    }

    public void setSpl(Spliceosome spl) {
        this.spl = spl;
    }

    public boolean isRnaP() {
        return rnaP;
    }

    public void setRnaP(boolean rnaP) {
        this.rnaP = rnaP;
    }

    // Override
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((spl == null) ? 0 : spl.hashCode());
        result = prime * result + splSitesA;
        result = prime * result + splSitesD;
        result = prime * result + splSitesTT;
        result = prime * result + (spliceosome ? 1231 : 1237);
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
        PrecursormRNA other = (PrecursormRNA) obj;
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
        if (spliceosome != other.spliceosome)
            return false;
        return true;
    }
}
