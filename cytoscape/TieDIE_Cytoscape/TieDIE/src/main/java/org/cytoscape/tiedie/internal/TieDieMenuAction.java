package org.cytoscape.tiedie.internal;

import java.awt.event.ActionEvent;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.application.swing.CySwingApplication;

/**
 *  @author SrikanthB
 * 
 *  A new menu item named "TieDIE" is created under Apps section of Cytoscape 
 */

public class TieDieMenuAction extends AbstractCyAction{
    
    public CyApplicationManager cyApplicationManager;
    public CySwingApplication cyDesktopService;
    public CyActivator cyactivator;

    public TieDieMenuAction(CyApplicationManager cyApplicationManager, final String menuTitle, CyActivator cyactivator) {
        super(menuTitle, cyApplicationManager, null, null);
        setPreferredMenu("Apps");
        this.cyactivator = cyactivator;
        this.cyApplicationManager = cyApplicationManager;
        this.cyDesktopService = cyactivator.getcytoscapeDesktopService();
    }

    public void actionPerformed(ActionEvent e) {
        TieDieCore tiediecore = new TieDieCore(cyactivator);
    }
    
    
}
