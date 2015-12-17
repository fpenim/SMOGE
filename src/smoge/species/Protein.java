package smoge.species;

/**
 * Created by fpenim on 12/12/2015.
 */
public class Protein {

    private boolean complete;

    /**
     * Constructor
     */
    public Protein(){
        complete = false;
    }

    // Translation finished
    public boolean finished(){
        complete = true;
        return complete;
    }

    // Getters & Setters
    public boolean isComplete() {
        return complete;
    }

    public void setComplete(boolean complete) {
        this.complete = complete;
    }

    // Override
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + (complete ? 1231 : 1237);
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
        Protein other = (Protein) obj;
        if (complete != other.complete)
            return false;
        return true;
    }

}
