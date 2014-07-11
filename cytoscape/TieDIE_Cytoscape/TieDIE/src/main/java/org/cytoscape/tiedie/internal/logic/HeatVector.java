package org.cytoscape.tiedie.internal.logic;

import Jama.Matrix;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;

/**
 *
 * @author SrikanthB
 */
public class HeatVector {

    Matrix heatVectorOfScores;
    int numOfColumns;
    List<CyNode> nodeHeatList;
    Map<CyNode, Double> nodeScoreMap;

    public HeatVector(int numOfColumns) {
        this.numOfColumns = numOfColumns;
        this.heatVectorOfScores = new Matrix(1, numOfColumns);
    }

    public HeatVector(Matrix rowVector) {
        this.heatVectorOfScores = rowVector;
        this.numOfColumns = rowVector.getColumnDimension();
    }

    public HeatVector extractHeatVector(String columnName, List<CyNode> nodeList, CyTable nodeTable) {
        int counter = 0;
        double heatscore;
        nodeScoreMap = new HashMap<CyNode, Double>();

        for (CyNode root : nodeList) {
            CyRow row = nodeTable.getRow(root.getSUID());
            if (row.get(columnName, Double.class) != null) {
                heatscore = row.get(columnName, Double.class);
                heatVectorOfScores.set(0, counter, heatscore);
                nodeHeatList.add(root); // nodeHeatList has all nodes corresponding to that heat column
                nodeScoreMap.put(root, heatscore);
            }

            counter++;
        }

        return this;
    }

}
