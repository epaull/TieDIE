package org.cytoscape.tiedie.internal.logic;

import Jama.Matrix;
import org.cytoscape.model.CyEdge;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.NetworkTestSupport;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author SrikanthB
 */
public class KernelTest {
    
    static CyNetwork testNetwork =null;
    static double[][] A = {{0,1},{1,0}};
    static double[][] D = {{1,0},{0,1}};
    static double[][] L = {{1,-1},{-1,1}};
    static double[][] K = {{0.909365,0.090635},{0.090635,0.909365}};
    
    public KernelTest() {
    }
    
    @BeforeClass
    public static void oneTimeSetUpClass() {
    }
    
    @AfterClass
    public static void oneTimeTearDownClass() {
    }
    
    @Before
    public void setUp() {
        final NetworkTestSupport nts = new NetworkTestSupport();
        testNetwork = nts.getNetwork();
        // Now you have a CyNetwork!
        // Set name for network
        testNetwork.getRow(testNetwork).set(CyNetwork.NAME, "Test Network");
        // Add two nodes to the network
        CyNode node1 = testNetwork.addNode();
        CyNode node2 = testNetwork.addNode();
        // Set name for new nodes
        testNetwork.getRow(node1).set(CyNetwork.NAME, "Node1");
        testNetwork.getRow(node1).set(CyNetwork.NAME, "Node2");
        // Add an edge
        CyEdge e = testNetwork.addEdge(node1, node2, true);
        testNetwork.getDefaultEdgeTable().getRow(e.getSUID()).set(CyEdge.INTERACTION, "-a>");;
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of getTime method, of class Kernel.
     */
    @Test
    public void testGetTime() {
        System.out.println("getTime");
        Kernel instance = new Kernel(testNetwork);
        double expResult = 0.1;
        double result = instance.getTime();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of getadjacencyMatrixOfNetwork method, of class Kernel.
     */
    @Test
    public void testGetadjacencyMatrixOfNetwork() {
        System.out.println("getadjacencyMatrixOfNetwork");
        Kernel instance = new Kernel(testNetwork);
        double[][] expResult = A;
        double[][] result = instance.getadjacencyMatrixOfNetwork();
        assertArrayEquals(expResult, result);
    }

    /**
     * Test of getdiffusionKernelOfNetwork method, of class Kernel.
     */
    @Test
    public void testGetdiffusionKernelOfNetwork() {
        System.out.println("getdiffusionKernelOfNetwork");
        Kernel instance = new Kernel(testNetwork);
        double[][] expResult = K;
        double[][] r = instance.getdiffusionKernelOfNetwork();
        double rounded[][] = {{(double) Math.round(r[0][0] * 1000000) / 1000000,(double) Math.round(r[0][1] * 1000000) / 1000000},{(double) Math.round(r[1][0] * 1000000) / 1000000,(double) Math.round(r[1][1] * 1000000) / 1000000}};
        assertArrayEquals(expResult, rounded);
    }

    /**
     * Test of createAdjMatrix method, of class Kernel.
     */
    @Test
    public void testCreateAdjMatrix() {
        System.out.println("createAdjMatrix");
        Kernel instance = new Kernel(testNetwork);
        double[][] expResult = A;
        double[][] result = instance.createAdjMatrix();
        assertArrayEquals(expResult, result);
    }

    /**
     * Test of createDegMatrix method, of class Kernel.
     */
    @Test
    public void testCreateDegMatrix() {
        System.out.println("createDegMatrix");
        Kernel instance = new Kernel(testNetwork);
        double[][] expResult = D;
        double[][] result = instance.createDegMatrix();
        assertArrayEquals(expResult, result);
    }

    /**
     * Test of createLapMatrix method, of class Kernel.
     */
    @Test
    public void testCreateLapMatrix() {
        System.out.println("createLapMatrix");
        Kernel instance = new Kernel(testNetwork);
        instance.createAdjMatrix();
        instance.createDegMatrix();
        double[][] resMatrix = L;
        Matrix result = instance.createLapMatrix();
        assertArrayEquals(resMatrix, result.getArray());
    }

    /**
     * Test of createRequiredExponential method, of class Kernel.
     */
    @Test
    public void testCreateRequiredExponential() {
        System.out.println("createRequiredExponential");
        Kernel instance = new Kernel(testNetwork);
        double[][] expResult = K;
        double[][] r = instance.createRequiredExponential();
        double rounded[][] = {{(double) Math.round(r[0][0] * 1000000) / 1000000,(double) Math.round(r[0][1] * 1000000) / 1000000},{(double) Math.round(r[1][0] * 1000000) / 1000000,(double) Math.round(r[1][1] * 1000000) / 1000000}};
        assertArrayEquals(expResult, rounded);
    }

    /**
     * Test of diffuse method, of class Kernel.
     */
    @Test
    public void testDiffuse() {//33.53776945690825,56.46223054309174
        System.out.println("diffuse");
        double [][] vector = {{31,59}};
        HeatVector inputVector = new HeatVector(new Matrix(vector));
        Kernel instance = new Kernel(testNetwork);
        double[][] expResult = {{33.537769,56.462231}};
        DiffusedHeatVector result = instance.diffuse(inputVector);
        double r[][] = result.getVectorOfScores().getArray();
        double rounded[][] = {{(double) Math.round(r[0][0] * 1000000) / 1000000,(double) Math.round(r[0][1] * 1000000) / 1000000}};
        
        assertArrayEquals(expResult, rounded);
    }
    
}
