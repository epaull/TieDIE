
package org.cytoscape.tiedie.internal.logic;

import Jama.Matrix;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;


/**
 * @author SrikanthB
 * 
 */

public class Kernel {
  
    private static final double t = 0.1;
    
    
    
    public Kernel() {
       
    }

    public double getTime(){
        return t;
    }
    
    
    public static double[][] createAdjMatrix(CyNetwork currentnetwork, List<CyNode> nodeList, CyTable edgeTable, int totalnodecount) {
        //make an adjacencymatrix for the current network
        double[][] adjacencyMatrixOfNetwork = new double[totalnodecount][totalnodecount];
        int k = 0;
        for (CyNode root : nodeList) {
            List<CyNode> neighbors = currentnetwork.getNeighborList(root, CyEdge.Type.OUTGOING);
            for (CyNode neighbor : neighbors) {
                List<CyEdge> edges = currentnetwork.getConnectingEdgeList(root, neighbor, CyEdge.Type.DIRECTED);
                if (edges.size() > 0) {
                    CyRow row = edgeTable.getRow(edges.get(0).getSUID());
                    try {
                        adjacencyMatrixOfNetwork[k][nodeList.indexOf(neighbor)] = Double.parseDouble(row.get(CyEdge.INTERACTION, String.class));
                    } catch (Exception ex) {
                    }
                }
            }
            k++;
        }
        return adjacencyMatrixOfNetwork;
    }

    
    
    
    public static double[][] createDegMatrix(CyNetwork currentnetwork, List<CyNode> nodeList, int totalnodecount) {

        double[][] degreeMatrixOfNetwork = new double[totalnodecount][totalnodecount];

        int k = 0;
        for (CyNode root : nodeList) {
            List<CyNode> myNeighbourList = currentnetwork.getNeighborList(root, CyEdge.Type.ANY);
            int myNodeDegree = myNeighbourList.size();
            degreeMatrixOfNetwork[k][k] = myNodeDegree;
            k++;
        }

        return degreeMatrixOfNetwork;
    }

    
    public static double[][] createLapMatrix(double[][] adjacencyMatrixOfNetwork, double[][] degreeMatrixOfNetwork, int totalnodecount) {

        double[][] laplacianMatrixOfNetwork;
        Matrix D = new Matrix(degreeMatrixOfNetwork);
        Matrix A = new Matrix(adjacencyMatrixOfNetwork);
        Matrix L = D.minus(A);
        laplacianMatrixOfNetwork = L.getArrayCopy();
        
        return laplacianMatrixOfNetwork;
    }

    
    public static double[][] createRequiredExponential(double[][] laplacianMatrixOfNetwork){
        
        double[][] minusOftL;
        double[][] diffusionKernel;
        DoubleMatrix diffusionKernelMatrix;
        
        Matrix C = new Matrix(laplacianMatrixOfNetwork);
        C = C.timesEquals(-t);  //  (-t)*L
        minusOftL = C.getArrayCopy();
        
        diffusionKernelMatrix = new DoubleMatrix(minusOftL);
        diffusionKernelMatrix = MatrixFunctions.expm(diffusionKernelMatrix); // exponentiation 
        diffusionKernel = diffusionKernelMatrix.toArray2();
        
        return diffusionKernel;
    
    }
    
    
    
    public static HeatVector diffuse(HeatVector inputVector, double[][] diffusionKernel){
        Matrix diffusedVectorMatrix; 
        HeatVector diffusedOutputRowVector;
                 
        Matrix dKernelMatrix = new Matrix(diffusionKernel);
        diffusedVectorMatrix = inputVector.heatVectorOfScores.times(dKernelMatrix);
        diffusedOutputRowVector= new HeatVector(diffusedVectorMatrix);
        return diffusedOutputRowVector;
    }
    
    public static Map getnodeDiffusedScoreMap(HeatVector diffusedOutputRowVector, List<CyNode> nodeList){
        
        Map nodeDiffusedScoreMap; 
        nodeDiffusedScoreMap = new HashMap<Double, CyNode>();
        int count=0;
        for(CyNode root : nodeList){
            nodeDiffusedScoreMap.put(diffusedOutputRowVector.heatVectorOfScores.get(0,count), root);
            count++;
        }
    
        return nodeDiffusedScoreMap;
    }
    
   
    
 
}

    
    
    
    
    
    

