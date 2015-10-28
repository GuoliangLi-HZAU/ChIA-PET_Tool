/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.List;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class AnchorCluster implements Comparable {

    private ANCHOR head;
    private ANCHOR tail;
    private int iPets_between;
    private double score;
    private int bSearched = 0;
    private String petIndexes = "";

    public AnchorCluster(ANCHOR head, ANCHOR tail, int iPets_between) {
        this.head = head;
        this.tail = tail;
        this.iPets_between = iPets_between;
        this.score = 0.0;
        this.bSearched = 0;
        this.petIndexes = "";
        sortHeadTail();
    }

    public AnchorCluster(CLUSTER cluster, int iPets_between) {
        this.head = new ANCHOR(cluster.getHead());
        this.tail = new ANCHOR(cluster.getTail());
        this.iPets_between = iPets_between;
        this.score = 0.0;
        this.bSearched = 0;
        this.setPetIndexes(cluster.getPetIndexes());
        sortHeadTail();
    }

    public AnchorCluster(CLUSTER cluster, int iPets_between, int sortingLabel) {
        this.head = new ANCHOR(cluster.getHead());
        this.tail = new ANCHOR(cluster.getTail());
        this.iPets_between = iPets_between;
        this.score = 0.0;
        this.bSearched = 0;
        this.setPetIndexes(cluster.getPetIndexes());

        if (sortingLabel == 1) {
            sortHeadTail();
        } else if (sortingLabel == -1) {
            sortHeadTail_descending();
        }
    }

    private void sortHeadTail() {
        if (getHead().compareTo(getTail()) > 0) {
            ANCHOR tempAnchor = getHead();
            setHead(getTail());
            setTail(tempAnchor);
        }
    }

    private void sortHeadTail_descending() {
        if (getHead().compareTo(getTail()) < 0) {
            ANCHOR tempAnchor = getHead();
            setHead(getTail());
            setTail(tempAnchor);
        }
    }

    @Override
    public String toString() {
        int sameChrom = 0;
        if (getHead().getChrom().compareTo(getTail().getChrom()) == 0) {
            sameChrom = 1;
        }
        return (getHead().toString() + "\t" + getTail().toString() + "\t" + this.getiPets_between() + "\t" + sameChrom + "\t" + ANCHOR.calculateDistance(getHead(), getTail()) + "\t" + this.getScore() + "\t" + this.getPetIndexes());
    }

    public int compareTo(Object anotherAnchorCluster) throws ClassCastException {
        if (!(anotherAnchorCluster instanceof AnchorCluster)) {
            throw new ClassCastException("An AnchorCluster object expected.");
        }
        int result = this.getHead().compareTo(((AnchorCluster) anotherAnchorCluster).getHead());
        if (result == 0) {
            result = this.getTail().compareTo(((AnchorCluster) anotherAnchorCluster).getTail());
        }
        return result;
    }

    public static Vector<AnchorCluster> load(String AnchorClusterFile) throws IOException {
        Vector<AnchorCluster> anchorClusters = new Vector<AnchorCluster>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(AnchorClusterFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            ANCHOR head = new ANCHOR(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
            ANCHOR tail = new ANCHOR(fields[5], Integer.parseInt(fields[6]), Integer.parseInt(fields[7]));
            AnchorCluster anchorCluster = new AnchorCluster(head, tail, Integer.parseInt(fields[10]));
            anchorClusters.add(anchorCluster);
        }
        fileIn.close();
        return anchorClusters;
    }

    public static void save(List<AnchorCluster> anchorClusters, String clusterFileName) throws IOException {
        PrintWriter clusterFileOut = new PrintWriter(new BufferedWriter(new FileWriter(clusterFileName)));
        for (AnchorCluster cluster : anchorClusters) {
            clusterFileOut.println(cluster.toString());
        }
        clusterFileOut.close();
    }

    /**
     * @return the head
     */
    public ANCHOR getHead() {
        return head;
    }

    /**
     * @param head the head to set
     */
    public void setHead(ANCHOR head) {
        this.head = head;
    }

    /**
     * @return the tail
     */
    public ANCHOR getTail() {
        return tail;
    }

    /**
     * @param tail the tail to set
     */
    public void setTail(ANCHOR tail) {
        this.tail = tail;
    }

    /**
     * @return the iPets_between
     */
    public int getiPets_between() {
        return iPets_between;
    }

    /**
     * @param iPets_between the iPets_between to set
     */
    public void setiPets_between(int iPets_between) {
        this.iPets_between = iPets_between;
    }

    /**
     * @return the score
     */
    public double getScore() {
        return score;
    }

    /**
     * @param score the score to set
     */
    public void setScore(double score) {
        this.score = score;
    }

    /**
     * @return the bSearched
     */
    public int getbSearched() {
        return bSearched;
    }

    /**
     * @param bSearched the bSearched to set
     */
    public void setbSearched(int bSearched) {
        this.bSearched = bSearched;
    }

    /**
     * @return the petIndexes
     */
    public String getPetIndexes() {
        return petIndexes;
    }

    /**
     * @param petIndexes the petIndexes to set
     */
    public void setPetIndexes(String petIndexes) {
        this.petIndexes = petIndexes;
    }
}
