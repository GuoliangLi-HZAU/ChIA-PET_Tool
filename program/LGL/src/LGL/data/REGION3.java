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
import java.util.Collections;
import java.util.HashSet;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class REGION3 extends REGION {

    private int nTags = 0;
    private int bSearched = 0;
    private HashSet<REGION3> linkedRegions = new HashSet<REGION3>();
    private HashSet<CLUSTER> linkedClusters = new HashSet<CLUSTER>();
    private int alternativeSplicing = 0;

    public REGION3(String chrom, int start, int end, int nTags) {
        this.setChrom(chrom);
        this.setStart(start);
        this.setEnd(end);
        this.setnTags(nTags);
        this.setAlternativeSplicing(0);
    }

    public REGION3(REGION region, int nTags) {
        this.setChrom(region.getChrom());
        this.setStart(region.getStart());
        this.setEnd(region.getEnd());
        this.setnTags(nTags);
        this.setAlternativeSplicing(0);
    }

    public REGION3(REGION3 region3) {
        this.setChrom(region3.getChrom());
        this.setStart(region3.getStart());
        this.setEnd(region3.getEnd());
        this.setnTags(region3.getnTags());
        this.setAlternativeSplicing(0);
    }

    // the separator is TAB
    @Override
    public String toString() {
        return new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getnTags());
    }

    @Override
    public String toString(int mode) {
        switch (mode) {
            case 1:
                return new String(this.getChrom() + ":" + this.getStart() + "-" + this.getEnd() + "(" + this.getnTags() + "," + this.getAlternativeSplicing() + ")");

            case 2:
                return new String(this.getChrom() + ":" + this.getStart() + "-" + this.getEnd() + "(" + this.getnTags() + ")");
            case 3:
                return new String(this.getChrom() + ":" + this.getStart() + "-" + this.getEnd() + "(" + this.getSpan() + ")");

            default:
                return new String(this.getChrom() + ":" + this.getStart() + "-" + this.getEnd() + "(" + this.getnTags() + "," + this.getAlternativeSplicing() + ")");
        }
    }

    /**
     * @return the nTags
     */
    public int getnTags() {
        return nTags;
    }

    /**
     * @param nTags the nTags to set
     */
    public void setnTags(int nTags) {
        this.nTags = nTags;
    }

    // cannot override the method load in REGION
    public static Vector<REGION3> loadRegion3(String Region3File) throws IOException {
        Vector<REGION3> region3s = new Vector<REGION3>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(Region3File))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.trim().length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.trim().split("\t");
            if (fields.length < 3) {
                continue;
            }
            int nTags = 1;
            if (fields.length >= 4) {
                nTags = Integer.parseInt(fields[3]);
            }
            REGION3 region3 = new REGION3(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), nTags);
            if (region3.getStart() > region3.getEnd()) {
                System.out.println("start is larger than end: " + region3.toString());
                region3.setStart(Integer.parseInt(fields[2]));
                region3.setEnd(Integer.parseInt(fields[1]));
                System.out.println("exchanged start and end");
            }
            region3s.add(region3);
        }
        fileIn.close();

        return region3s;
    }

    // cannot override the method load in REGION
    public static void saveRegion3(String RegionFile, Vector<REGION3> regions) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(RegionFile)));

        for (int i = 0; i < regions.size(); i++) {
            fileOut.println(regions.elementAt(i).toString());
        }

        fileOut.close();
    }

    public void combine(REGION3 anotherRegion3) throws Exception {
        if (this.getChrom().compareTo(anotherRegion3.getChrom()) != 0) {
            throw new Exception("Two region3s are not in the same chromosome");
        }

        int newStart = (this.getStart() < anotherRegion3.getStart()) ? this.getStart() : anotherRegion3.getStart();
        int newEnd = (this.getEnd() > anotherRegion3.getEnd()) ? this.getEnd() : anotherRegion3.getEnd();
        this.setStart(newStart);
        this.setEnd(newEnd);
        this.setnTags(this.getnTags() + anotherRegion3.getnTags());
    }

    public static Vector<REGION3> combineOverlappedRegion3s(Vector<REGION3> region3s) {
        Vector<REGION3> newRegion3s = new Vector<REGION3>();
        if (region3s.size() <= 0) {
            return newRegion3s;
        }
        Collections.sort(region3s);
        REGION3 currentRegion3 = region3s.elementAt(0);
        for (int iRegion3 = 1; iRegion3 < region3s.size(); iRegion3++) {
            REGION3 region3 = region3s.elementAt(iRegion3);
            if (currentRegion3.overlap(region3) == true) {
                try {
                    currentRegion3.combine(region3);
                } catch (Exception e) {
                    System.out.println("Exception: " + e.toString());
                }
            } else {
                newRegion3s.add(currentRegion3);
                currentRegion3 = region3;
            }
        }
        newRegion3s.add(currentRegion3);

        return newRegion3s;
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
     * @return the alternativeSplicing
     */
    public int getAlternativeSplicing() {
        return alternativeSplicing;
    }

    /**
     * @param alternativeSplicing the alternativeSplicing to set
     */
    public void setAlternativeSplicing(int alternativeSplicing) {
        this.alternativeSplicing = alternativeSplicing;
    }

    /**
     * @return the linkedRegions
     */
    public HashSet<REGION3> getLinkedRegions() {
        return linkedRegions;
    }

    /**
     * @param linkedRegions the linkedRegions to set
     */
    public void setLinkedRegions(HashSet<REGION3> linkedRegions) {
        this.linkedRegions = linkedRegions;
    }

    /**
     * @return the linkedClusters
     */
    public HashSet<CLUSTER> getLinkedClusters() {
        return linkedClusters;
    }

    /**
     * @param linkedClusters the linkedClusters to set
     */
    public void setLinkedClusters(HashSet<CLUSTER> linkedClusters) {
        this.linkedClusters = linkedClusters;
    }
}
