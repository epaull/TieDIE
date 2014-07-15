package org.cytoscape.tiedie.internal.logic;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.cytoscape.model.CyNode;


/**
 *
 * @author SrikanthB
 * "TieDieUtil.java"  will have all utility functions used by "TieDieLogicThread.java"
 */

// Check that variables are not declared inside loop & check all maps,sets declarations
public class TieDieUtil {
    
      public static double findLinkerCutoff(List<CyNode> nodeList, Set<CyNode> upstreamnodeheatSet, Set<CyNode> downstreamnodeSet, Map upnodeScoreMapDiffused, Map downnodeScoreMapDiffused, double sizeFactor){
          double linker_cutoff;
          
          if(downnodeScoreMapDiffused == null) // hotnet algorithm
          {
              linker_cutoff = findLinkerCutoffSingle(nodeList, upstreamnodeheatSet, upnodeScoreMapDiffused, sizeFactor);
          }
          else
          {
              linker_cutoff = findLinkerCutoffMulti(nodeList, upstreamnodeheatSet, downstreamnodeSet, upnodeScoreMapDiffused, downnodeScoreMapDiffused, sizeFactor);
          }
          
          return linker_cutoff;
      } 
    
    
      public static double findLinkerCutoffSingle(List<CyNode> nodeList, Set<CyNode> upstreamnodeSet, Map upnodeScoreMapDiffused, double sizeFactor) {
          double linker_cutoff=0;
          double targetSize;
          double EPSILON = 0.0001;
          Map nodeDiffusedScoreMapSorted;
          Set<CyNode> diffusedSet;
      
          targetSize = (sizeFactor)*(upstreamnodeSet.size());
          nodeDiffusedScoreMapSorted = MapUtil.sortByValue(upnodeScoreMapDiffused);
          diffusedSet = new HashSet<CyNode>();
          
          Iterator<Map.Entry<CyNode, Double>> iterator = nodeDiffusedScoreMapSorted.entrySet().iterator() ;
          while(iterator.hasNext()){
              Entry<CyNode, Double> entry = iterator.next();
              linker_cutoff = entry.getValue()+EPSILON ;
              diffusedSet.add(entry.getKey());
              if((upstreamnodeSet.size())-(diffusedSet.size()) > targetSize)
              break;
          }
          return linker_cutoff;
   
      }
    
    
       public static double findLinkerCutoffMulti(List<CyNode> nodeList, Set<CyNode> upstreamnodeSet, Set<CyNode> downstreamnodeSet, Map upnodeScoreMapDiffused, Map downnodeScoreMapDiffused, double sizeFactor){
           double linker_cutoff=0;
           double EPSILON = 0.0001;
           double h, size_frac;
           Map upScoreMapDiffusedSorted,downScoreMapDiffusedSorted;
           //First linkers are found and are then filtered linkers
           Map linkers_nodeScoreMapSorted,linkers_nodeScoreMap, filtered_linkersNodeScoreMap;
           upScoreMapDiffusedSorted = MapUtil.sortByValue(upnodeScoreMapDiffused);
           downScoreMapDiffusedSorted = MapUtil.sortByValue(downnodeScoreMapDiffused);
           
           linkers_nodeScoreMap = findLinkersMap(nodeList, upstreamnodeSet, downstreamnodeSet, upnodeScoreMapDiffused, downnodeScoreMapDiffused);
           filtered_linkersNodeScoreMap = findFilteredLinkersMap(linkers_nodeScoreMap, 1);
           
           linkers_nodeScoreMapSorted = MapUtil.sortByValue(linkers_nodeScoreMap);
           Iterator<Map.Entry<CyNode, Double>> iterator = linkers_nodeScoreMapSorted.entrySet().iterator() ;
           
           while(iterator.hasNext()){
               Entry<CyNode, Double> entry = iterator.next();
               h = entry.getValue();
               linker_cutoff = h-EPSILON;
               //size_frac = scoreLinkers(upstreamnodeDiffusedScoreMap, upstreamnodeDiffusedScoreMapSorted,downstreamnodeDiffusedScoreMap, downstreamnodeDiffusedScoreMapSorted ,upstreamnodeSet,downstreamnodeSet , linker_cutoff , sizeFactor);
           }
           
           return linker_cutoff;
       }
       /*
          "Z function" is used to combine score vectors for two input sets according to literature
          "filterLinkers" is the "Z function" of TieDIE python implementation & is done in 2 functions here
            1. findLinkerMap returns all linker nodes,scores
            2. findFilteredMap returns all linker nodes whose scores > cutoff
       */
       public static Map findLinkersMap(List<CyNode> nodeList,Set<CyNode> upstreamnodeSet, Set<CyNode> downstreamnodeSet,Map upnodeScoreMapDiffused, Map downnodeScoreMapDiffused){
           Map linkers_nodeScoreMap;
           double min_heat,x,y;
           
           linkers_nodeScoreMap = new HashMap<CyNode,Double>();
           if(downnodeScoreMapDiffused == null){
               return upnodeScoreMapDiffused;
           }
           for(CyNode root :upstreamnodeSet ){
               if(downstreamnodeSet.contains(root)==false)
               continue;
               x = (Double)upnodeScoreMapDiffused.get(root);
               x = (double)x;
               y = (Double)downnodeScoreMapDiffused.get(root);
               y = (double)y;
               min_heat = Math.min(x,y);
               linkers_nodeScoreMap.put(root, min_heat);
               
            // if(downstreamheatVectorDiffused.nodeHeatSet.contains(root)==false)
            // continue;
            // x = (Double)getnodeDiffusedScoreMap(upstreamheatVectorDiffused, nodeList).get(root);
            // y = (Double)getnodeDiffusedScoreMap(downstreamheatVectorDiffused, nodeList).get(root);
            // min_heat = Math.min(x,y);
            // linkers_nodeScoreMap.put(root,min_heat);
           }
           return linkers_nodeScoreMap;
       }
    
       public static Map findFilteredLinkersMap(Map linkers_nodeScoreMap, double cutoff){// see why not linker_
           Map filtered_linkersNodeScoreMap = new HashMap<CyNode, Double>();
           Iterator<Map.Entry<CyNode, Double>> iterator = linkers_nodeScoreMap.entrySet().iterator() ;
        
           while(iterator.hasNext()){
               Entry<CyNode, Double> entry = iterator.next();
               if(entry.getValue()> cutoff){
                   filtered_linkersNodeScoreMap.put(entry.getKey(), entry.getValue());
               }
           }
        return filtered_linkersNodeScoreMap;
       }
  
       
       
       
}
