
package org.cytoscape.tiedie.internal.logic;

import Jama.Matrix;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;


/**
 * @author SrikanthB
      "HeatKernel as a substitution for pageRank"
    -> Heat diffusion is actually the average of all random walks.
    -> PageRank is also a random walk model, but it continues until an equilibrium is reached whereas 
       the heat diffusion stops at that timepoint t, which we set to 0.1 (based on empirical evidence)
 
 */

public class Kernel {
 
    /*
        About instance variable :
        
        Time parameter should be 0.1. This is basically time for diffusion.
        t=0.1 is the optimal parameter from
        Reference : http://genome.cshlp.org/content/18/12/1991.full
    */
    
    private static final double t = 0.1;
    

    public double getTime(){
        return t;
    }
    
     /*
     About methods :
    
        1. "createAdjMatrix" creates adjacency matrix of input gene network using nodetable. The possible
        entries of adjacency matrix are +1, -1 ,0.
        2. "createDegMatrix" creates degree matrix of network which is a diagonal matrix
        3. "createLapMatrix" creates Laplacian matrix of network. 
            LaplacianMatrixOfNetwork = DegreeMatrix - AdjacencyMatrix
            L = D-A  
            Reference : http://en.wikipedia.org/wiki/Laplacian_matrix
        4. "createRequiredExponential" creates the required exponential 
            which is e ^(-t*L). 
            where t is time of diffusion = 0.1,
            L is laplacian matrix of network
            
            e^matrix is a matrix again!  jblas library has been used for exponentiation of matrix.
            Reference : http://mikiobraun.github.io/jblas/javadoc/org/jblas/MatrixFunctions.html#expm(org.jblas.DoubleMatrix)
            A scaled Pade approximation algorithm is used by the library.
            Running time of the Pade approximation is O(n^2) worst case, 
            and this this is going to be the slowest part of the TieDIE algorithm. 
    
        5. "diffuse" method takes inputheatvector and returns outputrowvector.This is the core to TieDIE algorithm.
            The input sets only have some nodes, but the output heat vector is over all nodes in the network.
            It's because the inputs from other methods only give confidence values for some of the genes.
            The advantage of tiedie is to take that vector, and make it into a diffused vector that has continuous
            values over all genes. it's 'fixing' that deficiency in the input data.
    
        6. Now that we have diffused scores for each node , "getnodeDiffusedScoreMap" returns 
           map <CyNode, diffusedHeatScore> 
    
    
    */
    
    public static double[][] createAdjMatrix(CyNetwork currentnetwork, List<CyNode> nodeList, CyTable edgeTable, int totalnodecount) {
        //make an adjacencymatrix for the current network
        double[][] adjacencyMatrixOfNetwork = new double[totalnodecount][totalnodecount];
        String natureOfInteraction;
        CyRow row;
        int k = 0;
        for (CyNode root : nodeList) {
            List<CyNode> neighbors = currentnetwork.getNeighborList(root, CyEdge.Type.OUTGOING);
            for (CyNode neighbor : neighbors) {
                List<CyEdge> edges = currentnetwork.getConnectingEdgeList(root, neighbor, CyEdge.Type.DIRECTED);
                if (edges.size() > 0) {
                    row = edgeTable.getRow(edges.get(0).getSUID());
                    try {
                        /*  CyEdge.INTERACTION ----  A String column created by default for every CyEdge that holds the interaction description of the edge.
                            static final String INTERACTION
                        
                            Conventions :
                            "a = post-transcriptional, t = transcriptional; '>' = activating, '|' = inactivating"
                        */
                        natureOfInteraction = row.get(CyEdge.INTERACTION, String.class);
                       
                        if("-a>".equals(natureOfInteraction) || "-t>".equals(natureOfInteraction)){
                            adjacencyMatrixOfNetwork[k][nodeList.indexOf(neighbor)] = 1;
                        }
                        else if("-a|".equals(natureOfInteraction) || "-t|".equals(natureOfInteraction)){
                            adjacencyMatrixOfNetwork[k][nodeList.indexOf(neighbor)] = -1;
                        }
                        else
                        {
                            adjacencyMatrixOfNetwork[k][nodeList.indexOf(neighbor)] = 0;
                        }
                         
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
        nodeDiffusedScoreMap = new LinkedHashMap<CyNode,Double>();
        int count=0;
        for(CyNode root : nodeList){
            nodeDiffusedScoreMap.put(root, diffusedOutputRowVector.heatVectorOfScores.get(0,count));
            count++;
        }
   
        return nodeDiffusedScoreMap;
    }
    
   
}

    
    
    
    
    
    

