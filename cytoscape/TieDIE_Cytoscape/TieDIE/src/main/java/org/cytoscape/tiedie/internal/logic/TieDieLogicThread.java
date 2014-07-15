package org.cytoscape.tiedie.internal.logic;

import java.util.List;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyTable;
import org.cytoscape.view.model.CyNetworkView;


/**
 * @author SrikanthB
  
  TieDieLogicThread.java  keeps the ball rolling according to our execution flow and 
  this will be called from start button of GUI. 
 */


public class TieDieLogicThread extends Thread {
    
    /*
        About instance variables :
    
        current network is current selected network by user.
        currentnetworkview is current network view of corresponding network
    */
    public CyNetwork currentnetwork;
    public CyNetworkView currentnetworkview;
     
    public TieDieLogicThread(CyNetwork currentnetwork, CyNetworkView currentnetworkview) {
        this.currentnetwork = currentnetwork;
        this.currentnetworkview = currentnetworkview;
    }

    @Override
    public void run(){
        
        double sizeFactor, linker_cutoff;
        List<CyNode> nodeList = currentnetwork.getNodeList(); // get list of all nodes using currentnetwork
        int totalnodecount = nodeList.size();
        CyTable edgeTable = currentnetwork.getDefaultEdgeTable(); // get corresponding nodetable
        CyTable nodeTable = currentnetwork.getDefaultNodeTable(); // get corresponding edgetable
        
        /*
        Get adjacency matrix A , Degree matrix D
        Laplacian matrix L = D-A
        Required exponentiation  e^(-t*L)   where t is time of diffusion
        */
        double[][] adjacencyMatrixOfNetwork = Kernel.createAdjMatrix(currentnetwork, nodeList, edgeTable, totalnodecount);
        double[][] degreeMatrixOfNetwork = Kernel.createDegMatrix(currentnetwork, nodeList, totalnodecount);
        double[][] laplacianMatrixOfNetwork = Kernel.createLapMatrix(adjacencyMatrixOfNetwork, degreeMatrixOfNetwork, totalnodecount);
        double[][] diffusionKernel = Kernel.createRequiredExponential(laplacianMatrixOfNetwork);
        
        /*
        Create upstreamheatVector, downstreamheatVector for 2-way diffusion
        "Extract them using extractHeatVector"
        */
        HeatVector upstreamheatVector = new HeatVector(totalnodecount);
        upstreamheatVector = upstreamheatVector.extractHeatVector("upstreamheat",nodeList,nodeTable);
        HeatVector downstreamheatVector = new HeatVector(totalnodecount);
        downstreamheatVector = downstreamheatVector.extractHeatVector("downstreamheat",nodeList,nodeTable);
        
        // Get the diffused heat vectors which spread all over the network
        HeatVector upstreamheatVectorDiffused = new HeatVector(totalnodecount);
        upstreamheatVectorDiffused =  Kernel.diffuse(upstreamheatVector, diffusionKernel);
        HeatVector downstreamheatVectorDiffused = new HeatVector(totalnodecount);
        downstreamheatVectorDiffused = Kernel.diffuse(downstreamheatVector, diffusionKernel); 
        
        sizeFactor = 1;
        linker_cutoff = TieDieUtil.findLinkerCutoff(nodeList, upstreamheatVector.getnodeHeatSet(), downstreamheatVector.getnodeHeatSet(), upstreamheatVectorDiffused, downstreamheatVectorDiffused, sizeFactor);
        // nodeList is the extra parameter to existing tiedie  https://github.com/epaull/TieDIE/blob/master/lib/tiedie_util.py#L336
        
        
        
        
        
        
        
    } 
    
     
    
}
