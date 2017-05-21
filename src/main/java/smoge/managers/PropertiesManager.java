package smoge.managers;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

public class PropertiesManager {
    private static final Logger log = Logger.getLogger(PropertiesManager.class.getName());

    private static PropertiesManager propertiesManager = null;
    private Properties properties = null;

    private PropertiesManager () {
        properties = new Properties();
        FileInputStream file;

        String filePath = "/Users/flaviapenim/IdeaProjects/SMOGE/src/main/resources/application.properties";
        try {
            file = new FileInputStream(filePath);
            properties.load(file);
            file.close();

        } catch (FileNotFoundException e) {
            log.log(Level.SEVERE, "Properties file not found: ", e);
            System.exit(0);
        } catch (IOException e) {
            log.log(Level.SEVERE, "Error reading propertiesManager file: ", e);
            System.exit(0);
        }
    }

    public synchronized static PropertiesManager getInstance(){
        if (propertiesManager == null) {
            propertiesManager = new PropertiesManager();
        }
        return propertiesManager;
    }

    public int getGeneLength() {
        String s = this.properties.getProperty("geneLength","1500");
        return Integer.parseInt(s);
    }

    public int getRnaPolNumber() {
        String s = this.properties.getProperty("rnaPolNumber","100");
        return Integer.parseInt(s);
    }

    public int getSpliceosomeNumber() {
        String s = this.properties.getProperty("spliceosomeNumber","100");
        return Integer.parseInt(s);
    }

    public int getRibosomeNumber() {
        String s = this.properties.getProperty("ribosomeNumber","100");
        return Integer.parseInt(s);
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("PropertiesManager{");
        sb.append("properties=").append(properties);
        sb.append('}');
        return sb.toString();
    }
}