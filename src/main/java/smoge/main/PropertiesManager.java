package main.java.smoge.main;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

public class PropertiesManager {
    private static final Logger log = Logger.getLogger(PropertiesManager.class.getName());

    private static PropertiesManager propertiesManager = null;
    private static Properties properties = null;

    private PropertiesManager () {
        properties = new Properties();
        FileInputStream file;

        String filePath = "./smoge.properties";
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

    /**
     * Returns a pointer to the sigleton properties.
     * @return properties
     */
    public synchronized static Properties getProperties() {
        if (propertiesManager == null) {
            propertiesManager = new PropertiesManager();
        }
        return properties;
    }

    //TODO get methods for properties
}
