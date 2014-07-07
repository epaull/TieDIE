package org.cytoscape.tiedie.internal.logic;

import java.util.List;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyTable;
import org.cytoscape.view.model.CyNetworkView;


/**
 * @author SrikanthB
 * 
 * 
 * 
 * 
 */


public class TieDieLogicThread extends Thread {
    
    public CyNetwork currentnetwork;
    public CyNetworkView currentnetworkview;
     
    public TieDieLogicThread(CyNetwork currentnetwork, CyNetworkView currentnetworkview) {
        this.currentnetwork = currentnetwork;
        this.currentnetworkview = currentnetworkview;
    }

    @Override
    public void run(){
        
        List<CyNode> nodeList = currentnetwork.getNodeList();
        int totalnodecount = nodeList.size();
        CyTable edgeTable = currentnetwork.getDefaultEdgeTable();
        CyTable nodeTable = currentnetwork.getDefaultNodeTable();
        
        double[][] adjacencyMatrixOfNetwork = Kernel.createAdjMatrix(currentnetwork, nodeList, edgeTable, totalnodecount);
        double[][] degreeMatrixOfNetwork = Kernel.createDegMatrix(currentnetwork, nodeList, totalnodecount);
        double[][] laplacianMatrixOfNetwork = Kernel.createLapMatrix(adjacencyMatrixOfNetwork, degreeMatrixOfNetwork, totalnodecount);
        double[][] diffusionKernel = Kernel.createRequiredExponential(laplacianMatrixOfNetwork);
        
        HeatVector upstreamheatVector = new HeatVector(totalnodecount);
        upstreamheatVector = upstreamheatVector.extractHeatVector("upstreamheat",nodeList,nodeTable);
        HeatVector downstreamheatVector = new HeatVector(totalnodecount);
        downstreamheatVector = downstreamheatVector.extractHeatVector("downstreamheat",nodeList,nodeTable);
             
        HeatVector upstreamheatVectorDiffused = new HeatVector(totalnodecount);
        upstreamheatVectorDiffused =  Kernel.diffuse(upstreamheatVector, diffusionKernel);
        HeatVector downstreamheatVectorDiffused = new HeatVector(totalnodecount);
        downstreamheatVectorDiffused = Kernel.diffuse(downstreamheatVector, diffusionKernel); 
        
        double sizeFactor = 1;
        double linker_cutoff = TieDieUtil.findLinkerCutoff(nodeList, upstreamheatVector.nodeHeatList, downstreamheatVector.nodeHeatList, upstreamheatVectorDiffused,downstreamheatVectorDiffused, sizeFactor);// nodeList is the extra parameter to existing tiedie
        List<CyNode> filteredNodeList = TieDieUtil.filterLinkers(upstreamheatVectorDiffused, downstreamheatVectorDiffused, linker_cutoff);
        
        
        
        
        
        
    } 
    
     
    
}
