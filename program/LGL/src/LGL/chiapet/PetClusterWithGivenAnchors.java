/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.chiapet;

import LGL.data.ANCHOR;
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
import java.util.Calendar;
import java.util.Collections;
import java.util.Vector;

/*
 * Output for debug
 * 1) densities in the bins, normalized by the number of iPets and the number of dummy iPets
 * 2) density matrix from 1)
 * 3) spanDensities
 * 4) disance density matrix
 * 5) prob matrix
 *
 */
/**
 *
 * @author ligl
 */
public class PetClusterWithGivenAnchors {

    int debugLevel = 4;
    Calendar rightNow = Calendar.getInstance();
    // The PETs are sorted by chromosomes and coordindates.
    // The head must be smaller than the tail; if not, exchange them with reverse function
    // With sorted PETs, it is easier to merge the PETs
    Vector<ANCHOR> anchors = new Vector<ANCHOR>();
    Vector anchorClusters = new Vector<AnchorCluster>();
    String prefix = null;
    int sortingLabel = 0;
    Vector<CLUSTER> clusters = new Vector<CLUSTER>();

    public PetClusterWithGivenAnchors(String preclusterFile, String anchorRegionFile, String prefix, int sortingLabel) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start PetClusterWithGivenAnchors ... ");

        this.sortingLabel = sortingLabel;
        // load anchors
        anchors = loadAnchors(anchorRegionFile);
        // output file
        this.prefix = prefix;
        String clusterFileName = new String(prefix + ".cluster.txt");
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
    }

    Vector<ANCHOR> loadAnchors(String regionFile) throws IOException {
        Vector<ANCHOR> anchorsTmp = new Vector<ANCHOR>();
        Vector<REGION> regions = REGION.load(regionFile);

        for (int i = 0; i < regions.size(); i++) {
            ANCHOR anchor = new ANCHOR(regions.get(i));
            anchorsTmp.add(anchor);
        }

        Collections.sort(anchorsTmp);

        return anchorsTmp;
    }

    public void outputClusters(PrintWriter clusterFileOut) throws IOException {
        for (int i = 0; i < anchorClusters.size(); i++) {
            AnchorCluster cluster = ((Vector<AnchorCluster>) anchorClusters).elementAt(i);
            clusterFileOut.println(cluster.toString());
        }
        clusterFileOut.flush();
    }

    public void generateClusters(Vector<CLUSTER> clusters) {
        anchorClusters.clear();

        // generate the anchorCluster list
        AnchorCluster anchorCluster;
        for (int i = 0; i < clusters.size(); i++) {
            if (debugLevel > 2) {
                if (i % 100000 == 0) {
                    rightNow = Calendar.getInstance();
                    System.out.println("[" + rightNow.getTime().toString() + "] " + (i / 1000) + "K iPets processed to generate anchor clusters");
                }
            }

            ANCHOR headAnchor = searchAnchor(anchors, new ANCHOR(clusters.elementAt(i).getHead()));
            ANCHOR tailAnchor = searchAnchor(anchors, new ANCHOR(clusters.elementAt(i).getTail()));

            if ((headAnchor != null) && (tailAnchor != null)) {
                anchorCluster = new AnchorCluster(headAnchor, tailAnchor, 1);
                anchorClusters.add(anchorCluster);
            }
        }
        anchorClusters = mergeAnchorClusters((Vector<AnchorCluster>) anchorClusters);
    }

    Vector<AnchorCluster> mergeAnchorClusters(Vector<AnchorCluster> anchorClusters) {
        if (anchorClusters.size() <= 1) {
            return anchorClusters;
        }

        Collections.sort(anchorClusters);
        Vector tempAnchorClusters = new Vector<AnchorCluster>();
        AnchorCluster anchorCluster = anchorClusters.elementAt(0);
        for (int i = 1; i < anchorClusters.size(); i++) {
            if (debugLevel > 2) {
                if (i % 100000 == 0) {
                    rightNow = Calendar.getInstance();
                    System.out.println("[" + rightNow.getTime().toString() + "] " + (i / 1000) + "K anchor clusters (of " + anchorClusters.size() + ") processed to combine anchor clusters");
                }
            }
            if (anchorCluster.compareTo(anchorClusters.elementAt(i)) == 0) {
                anchorCluster.setiPets_between(anchorCluster.getiPets_between() + anchorClusters.elementAt(i).getiPets_between());
            } else {
                tempAnchorClusters.add(anchorCluster);
                anchorCluster = anchorClusters.elementAt(i);
            }
        }
        tempAnchorClusters.add(anchorCluster);

        return tempAnchorClusters;
    }

    // find the anchor overlapped with the current region
    ANCHOR searchAnchor(Vector anchors, ANCHOR anchor) {
        int index = Collections.binarySearch(anchors, anchor);
        // the anchor is in the anchor list
        if (index >= 0) {
            return ((Vector<ANCHOR>) anchors).elementAt(index);
        }

        // the anchor is not in the anchor list
        index = -index - 1;
        if (index >= anchors.size()) {
            index = anchors.size() - 1;
        }

        //index--;
        while (index >= 0) {
            ANCHOR tmpAnchor = ((Vector<ANCHOR>) anchors).elementAt(index);
            if (tmpAnchor.overlap(anchor)) {
                return tmpAnchor;
            }
            if ((tmpAnchor.getChrom().compareTo(anchor.getChrom()) == 0) && (anchor.getStart() > tmpAnchor.getEnd())) {
                break;
            }
            index--;
        }

        // System.out.println("Error in anchor list\n" + "index = 0\t" + hit2.toString());
        return null;
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 4) {
            new PetClusterWithGivenAnchors(args[0], args[1], args[2], Integer.valueOf(args[3]));
        } else {
            System.out.println("Usage: java PetClusterWithGivenAnchors <preclusterFile> <given_anchor_file> <prefix> <sortingLabel>");
            System.out.println("  sortingLabel:  1 - sort the head and tail anchors in ascending order");
            System.out.println("                -1 - sort the head and tail anchors in ascending order");
            System.out.println("         otherwise - no sort (default)");
            System.exit(1);
        }
    }
}
