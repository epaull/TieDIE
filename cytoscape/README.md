TieDIE_Cytoscape
================



## For users interested in TieDIE app. for Cytoscape
*Stay tuned*



## For developers interested in TieDIE app. for Cytoscape (How to set up Netbeans IDE)
1. Install Netbeans [from here] (https://netbeans.org/downloads/)
2. Once installed ,Go to  File-> Open Project and this dialog appears [dialog like here] (http://s29.postimg.org/dqseabrhj/Capture.png)
3. You just have to select our TieDIE maven project which should be already downloaded(cloned) on your computer in advance from here [clone this repo](https://github.com/srikanthBezawada/TieDIE_Cytoscape.git)
4. Cytoscape javadocs are available here [3.0.0](http://chianti.ucsd.edu/cytoscape-3.0.0/API/)      [3.0.1](http://chianti.ucsd.edu/cytoscape-3.0.1/API/)
           [3.2.0](http://code.cytoscape.org/jenkins/job/cytoscape-3-javadoc/javadoc/)
5. More references for Cytoscape app. development [Cytoscape app. developer](http://wiki.cytoscape.org/Cytoscape_3/AppDeveloper/)
6. [Cytoscape App. cookbook] (http://wiki.cytoscape.org/Cytoscape_3/AppDeveloper/Cytoscape_3_App_Cookbook)
7. Once you are done with developing the project, you can clean and build the project and you'll find a jar file in the target folder of sourcecode. This jar file is used to install TieDIE app. through Cytoscape. You should download Cytoscape 3.x from here [Download Cytoscape](http://www.cytoscape.org/download.html)


## More about TieDIE 
1. ["TieDIE" journal](http://bioinformatics.oxfordjournals.org/content/29/21/2757.short)

2. [Existing TieDIE software](https://sysbiowiki.soe.ucsc.edu/tiedie)

3. [Existing TieDIE repo](https://github.com/epaull/TieDIE) (There is a cytoscape branch which has the code for tiedie app. hosted there as well)


## About this project

=======
[TieDIE  algorithm] (http://bioinformatics.oxfordjournals.org/content/29/21/2757.short) finds "subnetworks" connecting genomic perturbations to transcriptional changes in large gene interaction networks. This algorithm is useful in optimizing cancer treatment by viewing cancer from a gene network perspective thus enhancing our understanding of disease initiation, progression and therapy. This project deploys TieDIE algorithm into cytoscape as an app. 

This is a [GSOC 2014 project](https://www.google-melange.com/gsoc/project/details/google/gsoc2014/srikanthh/5733935958982656) developed by [Srikanth.B]
(https://github.com/srikanthBezawada) under the supervision of [Evan Paull](https://github.com/epaull)


