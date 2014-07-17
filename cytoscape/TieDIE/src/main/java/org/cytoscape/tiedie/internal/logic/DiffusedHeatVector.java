package org.cytoscape.tiedie.internal.logic;

import Jama.Matrix;

/**
 *
 * @author SrikanthB
 */
public class DiffusedHeatVector {
    private Matrix heatVectorOfScores;
    private int numOfColumns;
    
    public DiffusedHeatVector(int numOfColumns) {
        this.numOfColumns = numOfColumns;
        this.heatVectorOfScores = new Matrix(1, numOfColumns);
    }

    public DiffusedHeatVector(Matrix rowVector) {
        this.heatVectorOfScores = rowVector;
        this.numOfColumns = rowVector.getColumnDimension();
    }
    
    public Matrix getVectorOfScores(){
        return heatVectorOfScores;
    }
    
    public DiffusedHeatVector extractDiffusedHeatVector(HeatVector inputVector, Kernel heatDiffusionKernel){
        return heatDiffusionKernel.diffuse(inputVector);
    }    
    
}
