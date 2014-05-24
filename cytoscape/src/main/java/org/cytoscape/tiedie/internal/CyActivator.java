package org.cytoscape.tiedie.internal;

import java.util.Properties;
import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.swing.CySwingApplication;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.service.util.AbstractCyActivator;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.osgi.framework.BundleContext;

/**
 *
 * @author SrikanthB
 */

public class CyActivator extends AbstractCyActivator {

    public CyApplicationManager cyApplicationManager;
    public CySwingApplication cyDesktopService;
    public CyServiceRegistrar cyServiceRegistrar;
    public TieDieMenuAction menuaction;
    public static CyNetworkFactory networkFactory;
    public static CyNetworkManager networkManager;
    public static CyNetworkViewFactory networkViewFactory;
    public static CyNetworkViewManager networkViewManager;

    @Override
    public void start(BundleContext context) throws Exception {
        String version = new String(" 1.0");
        System.out.println("TieDIE app. is loading");
        this.networkViewManager = getService(context, CyNetworkViewManager.class);
        this.networkViewFactory = getService(context, CyNetworkViewFactory.class);
        this.networkFactory = getService(context, CyNetworkFactory.class);
        this.networkManager = getService(context, CyNetworkManager.class);
        this.cyApplicationManager = getService(context, CyApplicationManager.class);
        this.cyDesktopService = getService(context, CySwingApplication.class);
        this.cyServiceRegistrar = getService(context, CyServiceRegistrar.class);
        menuaction = new TieDieMenuAction(cyApplicationManager, "TieDIE" + version, this);
        Properties properties = new Properties();
        registerAllServices(context, menuaction, properties);
    }

    public CyServiceRegistrar getcyServiceRegistrar() {
        return cyServiceRegistrar;
    }

    public CyApplicationManager getcyApplicationManager() {
        return cyApplicationManager;
    }

    public CySwingApplication getcytoscapeDesktopService() {
        return cyDesktopService;
    }
    
    public TieDieMenuAction getmenuaction() {
        return menuaction;
    }
  
}
