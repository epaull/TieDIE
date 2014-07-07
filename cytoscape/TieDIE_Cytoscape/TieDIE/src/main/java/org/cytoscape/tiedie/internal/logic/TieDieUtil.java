package org.cytoscape.tiedie.internal.logic;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.cytoscape.model.CyNode;

/**
 *
 * @author SrikanthB
 */


public class TieDieUtil {

   
  
    
    public TieDieUtil(){
    
    }
    
    public static double findLinkerCutoff(List<CyNode> nodeList, List<CyNode> upstreamnodeheatList, List<CyNode> downstreamnodeheatList,HeatVector upstreamheatVectorDiffused,HeatVector downstreamheatVectorDiffused, double sizeFactor){
        double linker_cutoff;
        if(downstreamnodeheatList == null){
            linker_cutoff = findLinkerCutoffSingle(nodeList, upstreamnodeheatList, upstreamheatVectorDiffused, sizeFactor);
        }
           
        else {
            linker_cutoff = findLinkerCutoffMulti(nodeList, upstreamnodeheatList, downstreamnodeheatList, upstreamheatVectorDiffused, downstreamheatVectorDiffused, sizeFactor);
        }
         
        return linker_cutoff;
    }
    
    
    
    
    
    public static double findLinkerCutoffSingle(List<CyNode> nodeList, List<CyNode> upstreamnodeheatList, HeatVector upstreamheatVectorDiffused, double sizeFactor) {
        double linker_cutoff=0;
        double target_size;
        double EPSILON = 0.0001;
        target_size = (sizeFactor)*(upstreamnodeheatList.size());
        
        Map nodeDiffusedScoreMap;
        nodeDiffusedScoreMap = Kernel.getnodeDiffusedScoreMap(upstreamheatVectorDiffused, nodeList);
        
        
        
        
        
        
        return linker_cutoff;
    }
    
    
    
    public static double findLinkerCutoffMulti(List<CyNode> nodeList, List<CyNode> upstreamnodeheatList, List<CyNode> downstreamnodeheatList, HeatVector upstreamheatVectorDiffused, HeatVector downstreamheatVectorDiffused, double sizeFactor) {
        double linker_cutoff=0;
        
        
        
        
        
        
        return linker_cutoff;
    }

    
    
   
    public static List<CyNode> filterLinkers(HeatVector upstreamheatVectorDiffused, HeatVector downstreamheatVectorDiffused, double linker_cutoff){
        List<CyNode> filteredNodeList = null;
    
    
        return filteredNodeList;
    }
    
    
    
    
   
    
    
      
}
