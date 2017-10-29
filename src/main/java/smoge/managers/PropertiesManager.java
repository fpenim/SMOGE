package smoge.managers;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class PropertiesManager {

    private static PropertiesManager propertiesManager = null;
    private Properties properties = null;

    private PropertiesManager () {
        properties = new Properties();
        InputStream inputStream = this.getClass().getResourceAsStream("/application.properties");

        try {
            properties.load(inputStream);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public synchronized static PropertiesManager getInstance(){
        if (propertiesManager == null) {
            propertiesManager = new PropertiesManager();
        }
        return propertiesManager;
    }

    /**
     * General properties
     */
    public long getSimulationTime() {
        return Long.parseLong(this.properties.getProperty("simulationTime", "1000"));
    }

    public long getOutputStep() {
        return Long.parseLong(this.properties.getProperty("outputStep", "10"));
    }

    public int getGeneLength() {
        return Integer.parseInt(this.properties.getProperty("geneLength","1500"));
    }

    public int getRnaPolNumber() {
        return Integer.parseInt(this.properties.getProperty("rnaPolNumber","100"));
    }

    public int getSpliceosomeNumber() {
        return Integer.parseInt(this.properties.getProperty("spliceosomeNumber","100"));
    }

    public int getRibosomeNumber() {
        return Integer.parseInt(this.properties.getProperty("ribosomeNumber","100"));
    }

    public int getSpliceSitesNumber() {
        return Integer.parseInt(this.properties.getProperty("splicesSiteNumber", "0"));
    }

    /**
     * RNA Polymerase related properties
     */
    public double getPolymeraseKc() {
        return Double.parseDouble(this.properties.getProperty("polKc", "2"));
    }

    public double getPolymeraseKp() {
        return Double.parseDouble(this.properties.getProperty("polKp", "2"));
    }

    public double getPolymeraseKd() {
        return Double.parseDouble(this.properties.getProperty("polKd", "2"));
    }

    public double getPolymeraseKdg() {
        return Double.parseDouble(this.properties.getProperty("polKdg", "2"));
    }

    /**
     * Spliceosome related properties
     */
    public double getSpliceosomeKc() {
        return Double.parseDouble(this.properties.getProperty("splKc", "0.5"));
    }

    public double getSpliceosomeKs() {
        return Double.parseDouble(this.properties.getProperty("splKs", "0.5"));
    }

    public double getSpliceosomeKt() {
        return Double.parseDouble(this.properties.getProperty("splKt", "0.5"));
    }

    public double getSpliceosomeKdg() {
        return Double.parseDouble(this.properties.getProperty("splKdg", "0.5"));
    }

    /**
     * Ribosome related properties
     */
    public double getRibosomeKc() {
        return Double.parseDouble(this.properties.getProperty("ribKc", "1.5"));
    }

    public double getRibosomeKp() {
        return Double.parseDouble(this.properties.getProperty("ribKp", "1.5"));
    }

    public double getRibosomeKd() {
        return Double.parseDouble(this.properties.getProperty("ribKd", "1.5"));
    }

    public double getRibosomeKdg() {
        return Double.parseDouble(this.properties.getProperty("ribKdg", "1.5"));
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("PropertiesManager{");
        sb.append("properties=").append(properties);
        sb.append('}');
        return sb.toString();
    }
}