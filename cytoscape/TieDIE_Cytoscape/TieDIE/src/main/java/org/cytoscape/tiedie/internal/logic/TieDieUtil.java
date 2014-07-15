package org.cytoscape.tiedie.internal.logic;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.cytoscape.model.CyNode;

import static org.cytoscape.tiedie.internal.logic.Kernel.getnodeDiffusedScoreMap;

/**
 *
 * @author SrikanthB
 * "TieDieUtil.java"  will have all utility functions used by "TieDieLogicThread.java"
 */

// Check that variables are not declared inside loop & check all maps,sets declarations
public class TieDieUtil {
    
    public static double findLinkerCutoff(List<CyNode> nodeList, Set<CyNode> upstreamnodeheatSet, Set<CyNode> downstreamnodeSet,HeatVector upstreamheatVectorDiffused,HeatVector downstreamheatVectorDiffused, double sizeFactor){
        double linker_cutoff;
        if(downstreamnodeSet == null){
            linker_cutoff = findLinkerCutoffSingle(nodeList, upstreamnodeheatSet, upstreamheatVectorDiffused, sizeFactor);
        }
           
        else {
            linker_cutoff = findLinkerCutoffMulti(nodeList, upstreamnodeheatSet, downstreamnodeSet, upstreamheatVectorDiffused, downstreamheatVectorDiffused, sizeFactor);
        }
         
        return linker_cutoff;
    }
    
    
    public static double findLinkerCutoffSingle(List<CyNode> nodeList, Set<CyNode> upstreamnodeSet, HeatVector upstreamheatVectorDiffused, double sizeFactor) {
        double linker_cutoff=0;
        double target_size;
        double EPSILON = 0.0001;
        
        target_size = (sizeFactor)*(upstreamnodeSet.size());
        Map nodeDiffusedScoreMap, nodeDiffusedScoreMapSorted;
        Set<CyNode> diffused_node_set;
        
        nodeDiffusedScoreMap = Kernel.getnodeDiffusedScoreMap(upstreamheatVectorDiffused, nodeList);
        nodeDiffusedScoreMapSorted = MapUtil.sortByValue(nodeDiffusedScoreMap);
        
        diffused_node_set = new HashSet<CyNode>();
        Iterator<Map.Entry<CyNode, Double>> iterator = nodeDiffusedScoreMapSorted.entrySet().iterator() ;
        
        while(iterator.hasNext()){
            Entry<CyNode, Double> entry = iterator.next();
            linker_cutoff = entry.getValue()+EPSILON ;
            diffused_node_set.add(entry.getKey());
            if((upstreamnodeSet.size())-(diffused_node_set.size()) > target_size)
            break;
        }
        
        return linker_cutoff;
        
        /**
        Set<Map.Entry<CyNode, Double>> entrySet = nodeDiffusedScoreMapSorted.entrySet();
        Set<CyNode> diffused_node_set = new HashSet<CyNode>();
             
        for (Entry entry : entrySet) {
            linker_cutoff = (Double)entry.getValue()+EPSILON ;
            diffused_node_set.add((CyNode)entry.getKey());
            if((upstreamnodeheatSet.size())-(diffused_node_set.size()) > target_size)
            break;
        } 
        return linker_cutoff;
        */
        
    }
    
    
    
    public static double findLinkerCutoffMulti(List<CyNode> nodeList, Set<CyNode> upstreamnodeSet, Set<CyNode> downstreamnodeSet, HeatVector upstreamheatVectorDiffused, HeatVector downstreamheatVectorDiffused, double sizeFactor) {
        double linker_cutoff=0;
        double EPSILON = 0.0001;
        double h, size_frac;
        
        Map linkers_NodeScoreMapSorted,linkers_NodeScoreMap, filtered_linkersNodeScoreMap;
        Map upstreamnodeDiffusedScoreMap,downstreamnodeDiffusedScoreMap;
        Map upstreamnodeDiffusedScoreMapSorted, downstreamnodeDiffusedScoreMapSorted;
        
        upstreamnodeDiffusedScoreMap = Kernel.getnodeDiffusedScoreMap(upstreamheatVectorDiffused, nodeList);
        downstreamnodeDiffusedScoreMap = Kernel.getnodeDiffusedScoreMap(downstreamheatVectorDiffused, nodeList);
        upstreamnodeDiffusedScoreMapSorted = MapUtil.sortByValue(upstreamnodeDiffusedScoreMap);
        downstreamnodeDiffusedScoreMapSorted = MapUtil.sortByValue(downstreamnodeDiffusedScoreMap);
        
        linkers_NodeScoreMap = findLinkersMap(nodeList, upstreamheatVectorDiffused, downstreamheatVectorDiffused);
        filtered_linkersNodeScoreMap = findFilteredLinkersMap(linkers_NodeScoreMap, 1);
        
        linkers_NodeScoreMapSorted = MapUtil.sortByValue(linkers_NodeScoreMap);
        Iterator<Map.Entry<CyNode, Double>> iterator = linkers_NodeScoreMapSorted.entrySet().iterator() ;
        
  
        while(iterator.hasNext()){
            Entry<CyNode, Double> entry = iterator.next();
            h = entry.getValue();
            linker_cutoff = h-EPSILON;
            //size_frac = scoreLinkers(upstreamnodeDiffusedScoreMap, upstreamnodeDiffusedScoreMapSorted,downstreamnodeDiffusedScoreMap, downstreamnodeDiffusedScoreMapSorted ,upstreamnodeSet,downstreamnodeSet , linker_cutoff , sizeFactor);
        }
        
        
        return linker_cutoff;
    }
       
    /**
         "Z function" is used to combine score vectors for two input sets according to literature
         "filterLinkers" is the "Z function" of TieDIE python implementation & is done in 2 functions here
            1. findLinkerMap returns all linker nodes,scores
            2. findFilteredMap returns all linker nodes whose scores > cutoff
    */
   
    public static Map findLinkersMap(List<CyNode> nodeList, HeatVector upstreamheatVectorDiffused, HeatVector downstreamheatVectorDiffused){
       
        Map LinkersNodeScoreMap;
        double min_heat,x,y;
        
        if(downstreamheatVectorDiffused == null){
            LinkersNodeScoreMap = getnodeDiffusedScoreMap(upstreamheatVectorDiffused, nodeList);
            return LinkersNodeScoreMap;
        }
        LinkersNodeScoreMap = new HashMap<CyNode,Double>();
        for(CyNode root : upstreamheatVectorDiffused.nodeHeatSet)
        {
            if(downstreamheatVectorDiffused.nodeHeatSet.contains(root)==false)
            continue;
            x = (Double)getnodeDiffusedScoreMap(upstreamheatVectorDiffused, nodeList).get(root);
            y = (Double)getnodeDiffusedScoreMap(downstreamheatVectorDiffused, nodeList).get(root);
            min_heat = Math.min(x,y);
            LinkersNodeScoreMap.put(root,min_heat);
        }
        return LinkersNodeScoreMap;
    }
    
   
    public static Map findFilteredLinkersMap(Map LinkersNodeScoreMap, double cutoff){// see why not linker_
        Map filtered_linkersNodeScoreMap = new HashMap<CyNode, Double>();
        Iterator<Map.Entry<CyNode, Double>> iterator = LinkersNodeScoreMap.entrySet().iterator() ;
        
        while(iterator.hasNext()){
            Entry<CyNode, Double> entry = iterator.next();
            if(entry.getValue()> cutoff)
            {
                filtered_linkersNodeScoreMap.put(entry.getKey(), entry.getValue());
            }
        }
        return filtered_linkersNodeScoreMap;
    }
   
    
    
    public static double scoreLinkers(){
        return 0;
    
    
    
    }
    
    
    
    
    
    
    
    
    
}
