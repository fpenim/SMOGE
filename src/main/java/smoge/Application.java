package smoge;

import org.apache.log4j.Logger;
import smoge.managers.PropertiesManager;

import static smoge.Smoge.initializeEnvironment;
import static smoge.Smoge.startSimulation;

/**
 * Created by flaviapenim on 15/Jul/2017.
 */
public class Application {
    public static final Logger log = Logger.getLogger(Smoge.class);

    public static void main ( String[] args ) {
        log.info("Starting application...");

        PropertiesManager propertiesManager = PropertiesManager.getInstance();

        initializeEnvironment(propertiesManager);

        startSimulation(propertiesManager);

        log.info("Application terminated.");
    }
}
