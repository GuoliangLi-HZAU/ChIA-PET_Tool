/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.shortReads;

import LGL.data.REGION;
import LGL.dataConcise.RegionConcise;
import LGL.util.SeqUtil;
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
import java.util.Hashtable;
import java.util.Vector;

/**
 *
 * To keep the order of the regions in the input file, data structure REGION is
 * used To increase the search efficiency, the REGIONs are converted to
 * REGIONconcise
 *
 * @author ligl
 */
public class TagCountInGivenRegions {

    int debugLevel = 4;
    Calendar rightNow = Calendar.getInstance();
    Vector<REGION> regions = null;
    Hashtable<RegionConcise, REGION> hashRegionConcise2Region = new Hashtable<RegionConcise, REGION>();
    Hashtable<String, Integer> hashChrom2maxRegionSpans = new Hashtable<String, Integer>();
    Hashtable<String, ArrayList<RegionConcise>> hashChrom2regionConcise = null;
    Hashtable<REGION, Integer> hashRegion2tagCount = new Hashtable<REGION, Integer>();
    int extensionLength = 200;
    int extensionMode = 1;
    int minOverlapSize = 1;

    public TagCountInGivenRegions(String tagFile, String regionFile, String outputFile, int extensionLength, int extensionMode) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start TagCountInGivenRegions ... ");
        this.extensionLength = extensionLength;
        this.extensionMode = extensionMode;
        // load regions
        regions = REGION.load(regionFile);
        // convert regions to regionConcise
        // and the hashmap to the regions
        generateRegionConciseHash(regions);
        // generate the tag counts for each region: the number of tags in the regions
        generateTagCount(tagFile);
        // output the region clusters
        outputRegionTagCount(outputFile);
    }

    void outputRegionTagCount(String outputFile) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

        for (int i = 0; i < regions.size(); i++) {
            fileOut.println(regions.elementAt(i).toString() + "\t" + regions.elementAt(i).getAnnotation() + "\t" + hashRegion2tagCount.get(regions.elementAt(i)).intValue());
        }

        fileOut.close();
    }

    void generateTagCount(String tagFile) throws IOException {

        hashRegion2tagCount.clear();
        for (REGION region : regions) {
            hashRegion2tagCount.put(region, new Integer(0));
        }

        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(tagFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 0) {  // skip the short lines
                continue;
            }
            String fields[] = line.split("\t");
            updateTagCount(fields[0], Integer.parseInt(fields[1]), fields[2].charAt(0));
        }
        fileIn.close();
    }

    void updateTagCount(String chrom, int loci, char strand) {
        ArrayList regionConciseList = hashChrom2regionConcise.get(chrom);
        if (regionConciseList == null) {
            return;
        }
        int start = 0, end = 0;
        if (this.extensionMode == 2) { // extension in both directions
            start = loci - this.extensionLength;
            end = loci + this.extensionLength;
        } else { // extension from 5' to 3'
            if (SeqUtil.isForwardStrand(strand) == true) {
                start = loci;
                end = loci + this.extensionLength;
            } else {
                start = loci - this.extensionLength;
                end = loci;
            }
        }
        RegionConcise regionConcise = new RegionConcise(start, end);

        int index = Collections.binarySearch(regionConciseList, regionConcise);
        if (index < 0) {
            index = -index - 1;
        }
        if (index >= regionConciseList.size()) {
            index = regionConciseList.size() - 1;
        }

        if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(index)) >= minOverlapSize) {
            updateTagCount(((ArrayList<RegionConcise>) regionConciseList).get(index));
        }

        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (start - ((ArrayList<RegionConcise>) regionConciseList).get(iRegion).getStart() > hashChrom2maxRegionSpans.get(chrom).intValue()) {
                break;
            }
            if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(iRegion)) >= minOverlapSize) {
                updateTagCount(((ArrayList<RegionConcise>) regionConciseList).get(iRegion));
            }
        }
        for (int iRegion = index + 1; iRegion < regionConciseList.size(); iRegion++) {
            if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(iRegion)) >= minOverlapSize) {
                updateTagCount(((ArrayList<RegionConcise>) regionConciseList).get(iRegion));
            } else {
                break;
            }
        }
    }

    void updateTagCount(RegionConcise regionConcise) {
        REGION region = hashRegionConcise2Region.get(regionConcise);
        hashRegion2tagCount.put(region, new Integer(hashRegion2tagCount.get(region).intValue() + 1));
    }

    void generateRegionConciseHash(Vector<REGION> regions) {
        hashChrom2regionConcise = new Hashtable<String, ArrayList<RegionConcise>>();

        for (REGION region : regions) {
            updateChrom2regionConciseHash(region);
        }

        for (String chrom : hashChrom2regionConcise.keySet()) {
            // sort concise regions
            ArrayList<RegionConcise> regionConciseList = hashChrom2regionConcise.get(chrom);
            Collections.sort(regionConciseList);

            // get the max region span for each chromosome
            int maxRegionSpan = Integer.MIN_VALUE;
            for (RegionConcise region : regionConciseList) {
                int regionSpan = region.getSpan();
                if (maxRegionSpan < regionSpan) {
                    maxRegionSpan = regionSpan;
                }
            }
            hashChrom2maxRegionSpans.put(chrom, new Integer(maxRegionSpan));
        }
    }

    void updateChrom2regionConciseHash(REGION region) {
        String chrom = region.getChrom();
        ArrayList<RegionConcise> regionList = hashChrom2regionConcise.get(chrom);
        if (regionList == null) {
            regionList = new ArrayList<RegionConcise>();
            hashChrom2regionConcise.put(chrom, regionList);
        }
        RegionConcise regionConcise = new RegionConcise(region.getStart(), region.getEnd());
        regionList.add(regionConcise);

        hashRegionConcise2Region.put(regionConcise, region);
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 5) {
            new TagCountInGivenRegions(args[0], args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]));
        } else {
            System.out.println("Usage: java TagCountInGivenRegions <input_aln_file> <region_file> <output_file> <extensionLength> <extensionMode>");
            System.out.println("       <extensionMode>: 1 - from 5' to 3';  2 - in both directions");
            System.exit(1);
        }
    }
}
