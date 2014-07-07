package logictest;

import java.util.List;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;

import org.cytoscape.tiedie.internal.CyActivator;
import org.cytoscape.tiedie.internal.logic.Kernel;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 *
 * @author SrikanthB
 */
 
public class KernelTest {
    
    static CyNetwork currentnetwork =null;
    
    public KernelTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
//        final NetworkTestSupport nts = new NetworkTestSupport();
//        final CyNetwork network = nts.getNetwork();
        // Now you have a CyNetwork!
        // Create a new network
        currentnetwork = CyActivator.networkFactory.createNetwork();
        // Set name for network
        currentnetwork.getRow(currentnetwork).set(CyNetwork.NAME, "Test");
        // Add two nodes to the network
        CyNode node1 = currentnetwork.addNode();
        CyNode node2 = currentnetwork.addNode();
        // Set name for new nodes
        currentnetwork.getRow(node1).set(CyNetwork.NAME, "Node1");
        currentnetwork.getRow(node1).set(CyNetwork.NAME, "Node2");
        // Add an edge
        currentnetwork.addEdge(node1, node2, true);
        // Add the network to Cytoscape
        CyActivator.networkManager.addNetwork(currentnetwork);
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of getTime method, of class Kernel.
     */
    @Test
    public void DISABLED_testGetTime() {
        System.out.println("getTime");
        Kernel instance = new Kernel();
        double expResult = 0.0;
        double result = instance.getTime();
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of createAdjMatrix method, of class Kernel.
     */
    @Test
    public void testCreateAdjMatrix() {
        System.out.println("createAdjMatrix");
//        List<CyNode> nodeList = null;
//        CyTable edgeTable = null;
//        int totalnodecount = 0;
//        double[][] expResult = null;
//        double[][] result = Kernel.createAdjMatrix(currentnetwork, nodeList, edgeTable, totalnodecount);
//        assertArrayEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
    }

    /**
     * Test of createDegMatrix method, of class Kernel.
     */
    @Test
    public void DISABLED_testCreateDegMatrix() {
        System.out.println("createDegMatrix");
        CyNetwork currentnetwork = null;
        List<CyNode> nodeList = null;
        int totalnodecount = 0;
        double[][] expResult = null;
        double[][] result = Kernel.createDegMatrix(currentnetwork, nodeList, totalnodecount);
        assertArrayEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of createLapMatrix method, of class Kernel.
     */
    @Test
    public void DISABLED_testCreateLapMatrix() {
        System.out.println("createLapMatrix");
        double[][] adjacencyMatrixOfNetwork = null;
        double[][] degreeMatrixOfNetwork = null;
        int totalnodecount = 0;
        double[][] expResult = null;
        double[][] result = Kernel.createLapMatrix(adjacencyMatrixOfNetwork, degreeMatrixOfNetwork, totalnodecount);
        assertArrayEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of createRequiredExponential method, of class Kernel.
     */
    @Test
    public void DISABLED_testCreateRequiredExponential() {
        System.out.println("createRequiredExponential");
        double[][] laplacianMatrixOfNetwork = null;
        double[][] expResult = null;
        double[][] result = Kernel.createRequiredExponential(laplacianMatrixOfNetwork);
        assertArrayEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}
