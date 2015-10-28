/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.chiapet;

import LGL.data.AnchorCluster;
import LGL.data.CLUSTER;
import LGL.data.REGION;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Vector;

/*
 * Difference from PetCluster2 to PetCluster
 * 1) the input file needs to be sorted by
 *    a) the chromosome in the read1,
 *    b) the chromosome in the read2,
 *    c) the coordinate in the read1, ###### not required in current implementation
 *    d) the coordinate in the read2  ###### not required in current implementation
 *
 */
/**
 *
 * @author ligl
 */
public class PetCluster2 {

    int debugLevel = 4;
    int sortingLabel = 1; // ascending order [default]
    Calendar rightNow = Calendar.getInstance();
    // The PETs are sorted by chromosomes and coordindates.
    // The head must be smaller than the tail; if not, exchange them with reverse
    // function
    // With sorted PETs, it is easier to merge the PETs
    Vector<CLUSTER> clusters = new Vector<CLUSTER>();
    ArrayList anchorClusters = new ArrayList<AnchorCluster>();

    public PetCluster2(String preclusterFile, String clusterFileName, int sortingLabel) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start PetCluster2 ... ");

        this.sortingLabel = sortingLabel;

        PrintWriter clusterFileOut = new PrintWriter(new BufferedWriter(new FileWriter(clusterFileName)));

        // load inter-ligation PETs
        String chrom1 = "";
        String chrom2 = "";
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(preclusterFile))));
        String line;
        int nPETs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 0) // skip the short lines
            {
                continue;
            }
            String[] fields = line.split("\t");
            REGION head = new REGION(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
            REGION tail = new REGION(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5]));
            double weight = 1.0;
            if (fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            CLUSTER cluster = new CLUSTER(head, tail, weight, this.sortingLabel);
            if ((!chrom1.equalsIgnoreCase(cluster.getHead().getChrom())) || (!chrom2.equalsIgnoreCase(cluster.getTail().getChrom()))) {
                generateClusters(clusters);
                outputClusters(clusterFileOut);
                clusters.clear();
                chrom1 = cluster.getHead().getChrom();
                chrom2 = cluster.getTail().getChrom();
            }

            clusters.add(cluster);
            nPETs++;
            if (nPETs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000) + "M PETs read from " + preclusterFile);
            }
        }
        fileIn.close();
        generateClusters(clusters);
        outputClusters(clusterFileOut);

        clusterFileOut.close();
        System.out.println("number of clusters = " + anchorClusters.size());
    }

    public void outputClusters(PrintWriter clusterFileOut) throws IOException {
        for (int i = 0; i < anchorClusters.size(); i++) {
            AnchorCluster anchorCluster = ((ArrayList<AnchorCluster>) anchorClusters).get(i);
            clusterFileOut.println(anchorCluster.toString());
            if (i % 100000 == 0) {
                clusterFileOut.flush();
            }
        }
        clusterFileOut.flush();
    }

    public void generateClusters(Vector<CLUSTER> clusters) {
        anchorClusters.clear();
        for (int iCluster = 0; iCluster < clusters.size(); iCluster++) {
            anchorClusters.add(new AnchorCluster(clusters.elementAt(iCluster), (int) clusters.elementAt(iCluster).getWeight(), sortingLabel));
            //System.out.println(clusters.elementAt(iCluster) + " //  " +(int)clusters.elementAt(iCluster).getWeight() +  "  @@  " + sortingLabel);
        }

        // generate the anchorCluster list
        int oldAnchorClusterSize = -1;
        int newAnchorClusterSize = anchorClusters.size();
        while (oldAnchorClusterSize != newAnchorClusterSize) {
            oldAnchorClusterSize = newAnchorClusterSize;
            anchorClusters = mergeAnchorClusters(anchorClusters); // multiple merges are required, since some clusters after merge can be overlapped
            newAnchorClusterSize = anchorClusters.size();
        }
    }

    ArrayList<AnchorCluster> mergeAnchorClusters(ArrayList<AnchorCluster> oldAnchorClusters) {
        if (oldAnchorClusters.size() <= 1) {
            return oldAnchorClusters;
        }
        Collections.sort(oldAnchorClusters);

        // initialize the searched lables
        for (int i = 0; i < oldAnchorClusters.size(); i++) {
            oldAnchorClusters.get(i).setbSearched(0);
        }

        ArrayList<AnchorCluster> newAnchorClusters = new ArrayList<AnchorCluster>();

        int iAnchorCluster = 0;
        while (iAnchorCluster < oldAnchorClusters.size()) {
            if (debugLevel > 2) {
                if (iAnchorCluster % 100000 == 0) {
                    rightNow = Calendar.getInstance();
                    System.out.println("[" + rightNow.getTime().toString() + "] " + (iAnchorCluster / 1000) + "K anchor clusters (of " + oldAnchorClusters.size() + ") processed to combine anchor clusters");
                }
            }
            AnchorCluster anchorCluster = oldAnchorClusters.get(iAnchorCluster);
            anchorCluster.setbSearched(1);
            for (int i = iAnchorCluster + 1; i < oldAnchorClusters.size(); i++) {
                // no overlap between the two head anchors ==> new clusters
                if (anchorCluster.getHead().overlap(oldAnchorClusters.get(i).getHead()) == false) {
                    break;
                }

                // skip the ones already marked
                if (oldAnchorClusters.get(i).getbSearched() == 1) {
                    continue;
                }

                if (anchorCluster.getTail().overlap(oldAnchorClusters.get(i).getTail()) == true) {
                    // head overlap, and tail overlap
                    // combine head
                    try {
                        anchorCluster.getHead().combine(oldAnchorClusters.get(i).getHead());
                    } catch (Exception e) {
                        System.out.println("Exception " + e);
                    }
                    // combine tail
                    try {
                        anchorCluster.getTail().combine(oldAnchorClusters.get(i).getTail());
                    } catch (Exception e) {
                        System.out.println("Exception " + e);
                    }
                    // combine the iPets_between
                    anchorCluster.setiPets_between(anchorCluster.getiPets_between() + oldAnchorClusters.get(i).getiPets_between());
                    // set bSearched
                    oldAnchorClusters.get(i).setbSearched(1);
                }
                // no else, because the tail anchor in the next cluster can be far far away
            }
            newAnchorClusters.add(anchorCluster);

            // skip the clusters already marked
            while (iAnchorCluster < oldAnchorClusters.size()) {
                if (oldAnchorClusters.get(iAnchorCluster).getbSearched() == 1) {
                    iAnchorCluster++;
                } else {
                    break;
                }
            }
        }

        return newAnchorClusters;
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 2) {
            int sortingLabel = 1; // ascending order [default]
            new PetCluster2(args[0], args[1], sortingLabel);
        } else if (args.length == 3) {
            int sortingLabel = Integer.parseInt(args[2]);
            new PetCluster2(args[0], args[1], sortingLabel);
        } else {
            System.out.println("Usage: java PetCluster2 <pre_cluster_file> <clusterFileName> [<sortingLabel>]");
            System.out.println("[<sortingLabel>] is optional. The default is sorting the head and tail in ascending order");
            System.out.println("sortingLabel:  1: ascending order [default]");
            System.out.println("              -1: descending order");
            System.out.println("              otherwise: no sorting");
            System.exit(1);
        }
    }
}
