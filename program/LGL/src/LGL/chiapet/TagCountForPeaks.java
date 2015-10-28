/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.chiapet;

import LGL.data.PEAK;
import LGL.dataConcise.RegionConcise;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author ligl
 */
public class TagCountForPeaks {

    int debugLevel = 4;
    Calendar rightNow = Calendar.getInstance();
    List<PEAK> peaks;
    HashMap<RegionConcise, PEAK> hashRegionConcise2Peak = new HashMap<RegionConcise, PEAK>();
    HashMap<String, Integer> maxRegionSpans = new HashMap<String, Integer>();
    int minOverlapSize = 1;

    public TagCountForPeaks(String peakFile, String petFile, String outputFile) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start TagCountForPeaks ... ");
        // read clusters
        peaks = PEAK.loadPeaks(peakFile);
        // generate concise regions, and the hashmap to the peaks
        HashMap<String, ArrayList<RegionConcise>> peakHash = generatePeakHash(peaks);
        // update the peaks with the pet file: the number of tags in the peak regions
        updatePeaks(petFile, peakHash);
        // output the peak clusters
        PEAK.save(peaks, outputFile, 2);
    }

    void updatePeaks(String petFile, HashMap<String, ArrayList<RegionConcise>> peakHash) throws IOException {
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(petFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 0) {  // skip the short lines
                continue;
            }
            String fields[] = line.split("\t");
            if (fields[0].compareToIgnoreCase(fields[3]) == 0) { // intra-chromosomal
                updatePeaks(fields[0], Integer.parseInt(fields[1]), peakHash, 0);
                updatePeaks(fields[3], Integer.parseInt(fields[4]), peakHash, 0);
            } else {
                updatePeaks(fields[0], Integer.parseInt(fields[1]), peakHash, 1);
                updatePeaks(fields[3], Integer.parseInt(fields[4]), peakHash, 1);
            }
        }
        fileIn.close();
    }

    void updatePeaks(String chrom, int loci, HashMap<String, ArrayList<RegionConcise>> peakHash, int mode) {
        ArrayList regionConciseList = peakHash.get(chrom);
        if (regionConciseList == null) {
            return;
        }
        RegionConcise regionConcise = new RegionConcise(loci, loci + 1);

        int index = Collections.binarySearch(regionConciseList, regionConcise);
        if (index < 0) {
            index = -index - 1;
        }
        if (index >= regionConciseList.size()) {
            index = regionConciseList.size() - 1;
        }

        if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(index)) >= minOverlapSize) {
            updatePeaks(((ArrayList<RegionConcise>) regionConciseList).get(index), mode);
        }

        int start = loci;
        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (start - ((ArrayList<RegionConcise>) regionConciseList).get(iRegion).getStart() > maxRegionSpans.get(chrom).intValue()) {
                break;
            }
            if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(iRegion)) >= minOverlapSize) {
                updatePeaks(((ArrayList<RegionConcise>) regionConciseList).get(iRegion), mode);
            }
        }
        for (int iRegion = index + 1; iRegion < regionConciseList.size(); iRegion++) {
            if (regionConcise.overlappedSize(((ArrayList<RegionConcise>) regionConciseList).get(iRegion)) >= minOverlapSize) {
                updatePeaks(((ArrayList<RegionConcise>) regionConciseList).get(iRegion), mode);
            } else {
                break;
            }
        }
    }

    void updatePeaks(RegionConcise regionConciseTemp, int mode) {
        PEAK peak = hashRegionConcise2Peak.get(regionConciseTemp);
        if (mode == 0) { // intra-chromosomal
            peak.setiPets_intra(peak.getiPets_intra() + 1);
        } else {
            peak.setiPets_inter(peak.getiPets_inter() + 1);
        }
    }

    HashMap<String, ArrayList<RegionConcise>> generatePeakHash(List anchorClusters) {
        HashMap<String, ArrayList<RegionConcise>> peakHash = new HashMap<String, ArrayList<RegionConcise>>();

        for (PEAK peak : ((ArrayList<PEAK>) anchorClusters)) {
            updatePeakHash(peakHash, peak);
        }

        for (String chrom : peakHash.keySet()) {
            ArrayList<RegionConcise> regionConciseList = peakHash.get(chrom);
            Collections.sort(regionConciseList);
            int maxRegionSpan = Integer.MIN_VALUE;
            for (RegionConcise region : regionConciseList) {
                int regionSpan = region.getSpan();
                if (maxRegionSpan < regionSpan) {
                    maxRegionSpan = regionSpan;
                }
            }
            maxRegionSpans.put(chrom, new Integer(maxRegionSpan));
        }

        return peakHash;
    }

    void updatePeakHash(HashMap<String, ArrayList<RegionConcise>> peakHash, PEAK peak) {
        String chrom = peak.getChrom();
        ArrayList<RegionConcise> regionList = peakHash.get(chrom);
        if (regionList == null) {
            regionList = new ArrayList<RegionConcise>();
            peakHash.put(chrom, regionList);
        }
        RegionConcise regionConcise = new RegionConcise(peak.getStart(), peak.getEnd());
        regionList.add(regionConcise);
        hashRegionConcise2Peak.put(regionConcise, peak);
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 3) {
            new TagCountForPeaks(args[0], args[1], args[2]);
        } else {
            System.out.println("Usage: java TagCountForPeaks <peak_file> <pet_file> <output_file>");
            System.exit(1);
        }
    }
}
