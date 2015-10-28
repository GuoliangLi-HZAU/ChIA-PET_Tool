/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import LGL.dataConcise.Hit3Concise;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class PEAK extends REGION {

    private int summitStart;
    private int summitEnd;
    private int nTags;
    private int nTags_control = 0;
    private double fold_enrichment = 0;
    private double score = 0;  // = -10 * log10(pvalue)
    private int filtered;
    private int id = -1;
    private int iPets_intra;
    private int iPets_inter;

    public PEAK(String chrom, int start, int end, int summitStart, int summitEnd, int nTags) {
        super(chrom, start, end);
        this.setSummitStart(summitStart);
        this.setSummitEnd(summitEnd);
        this.setnTags(nTags);
        this.setiPets_inter(0);
        this.setiPets_intra(0);
    }

    public PEAK(String chrom, int start, int end, int summitStart, int summitEnd, int nTags, int id) {
        super(chrom, start, end);
        this.setSummitStart(summitStart);
        this.setSummitEnd(summitEnd);
        this.setnTags(nTags);
        this.setId(id);
        this.setiPets_inter(0);
        this.setiPets_intra(0);
    }

    // the separator is TAB
    @Override
    public String toString() {
        return new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getSummitStart() + "\t" + this.getSummitEnd() + "\t" + this.getnTags());
    }

    @Override
    public String toString(int mode) {
        return new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getSummitStart() + "\t" + this.getSummitEnd() + "\t" + this.getnTags() + "\t" + this.getnTags_control() + "\t" + this.getFold_enrichment());
    }

    /**
     *
     * @param mode 1: output format for peak calling 2: output format for tag
     * counts
     * @return
     */
    public String toString2(int mode) {
        String result = "";
        switch (mode) {
            case 1:
                result = new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getSummitStart() + "\t" + this.getSummitEnd() + "\t" + this.getnTags() + "\t" + this.getnTags_control() + "\t" + this.getFold_enrichment());
                break;
            case 2:
                result = new String(this.getId() + "\t" + this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getnTags() + "\t" + this.getiPets_intra() + "\t" + this.getiPets_inter());
                break;
            default:
                result = new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getSummitStart() + "\t" + this.getSummitEnd() + "\t" + this.getnTags() + "\t" + this.getnTags_control() + "\t" + this.getFold_enrichment());
                break;
        }
        return result;
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

    /**
     * @return the fold_enrichment
     */
    public double getFold_enrichment() {
        return fold_enrichment;
    }

    /**
     * @param fold_enrichment the fold_enrichment to set
     */
    public void setFold_enrichment(double fold_enrichment) {
        this.fold_enrichment = fold_enrichment;
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
     * @return the start
     */
    public int getSummitStart() {
        return summitStart;
    }

    /**
     * @param start the start to set
     */
    public void setSummitStart(int summitStart) {
        this.summitStart = summitStart;
    }

    /**
     * @return the end
     */
    public int getSummitEnd() {
        return summitEnd;
    }

    /**
     * @param end the end to set
     */
    public void setSummitEnd(int summitEnd) {
        this.summitEnd = summitEnd;
    }

    /**
     * @return the filtered
     */
    public int getFiltered() {
        return filtered;
    }

    /**
     * @param filtered the filtered to set
     */
    public void setFiltered(int filtered) {
        this.filtered = filtered;
    }

    public static Vector<PEAK> filterNearbyPeaks(Vector<PEAK> peaks, int minimumDistance) {
        Vector<PEAK> peaksTemp = peaks;
        int oldPeakNumber = -1;
        int newPeakNumber = peaks.size();
        while (oldPeakNumber != newPeakNumber) {
            oldPeakNumber = newPeakNumber;
            peaksTemp = PEAK.filterNearbyPeaksOnce(peaksTemp, minimumDistance);
            newPeakNumber = peaksTemp.size();
        }
        peaksTemp = PEAK.combineNearbyPeaks(peaksTemp, minimumDistance);
        return peaksTemp;
    }

    public static Vector<PEAK> filterNearbyPeaksOnce(Vector<PEAK> peaks, int minimumDistance) {
        if (peaks.size() <= 1) {
            return peaks;
        }

        // initialize the filtered label
        Collections.sort(peaks);
        int nPeaks = peaks.size();
        for (int i = 0; i < nPeaks; i++) {
            peaks.elementAt(i).setFiltered(0);
        }

        //forward filtering
        for (int i = 0; i < nPeaks - 1; i++) {
            for (int j = i + 1; j < nPeaks; j++) {
                if ((peaks.elementAt(i).getChrom().compareTo(peaks.elementAt(j).getChrom()) == 0) && (peaks.elementAt(j).getStart() - peaks.elementAt(i).getEnd() < minimumDistance)) {
                    if (peaks.elementAt(i).getnTags() < peaks.elementAt(j).getnTags()) {
                        peaks.elementAt(i).setFiltered(1);
                    } else if (peaks.elementAt(i).getnTags() > peaks.elementAt(j).getnTags()) {
                        peaks.elementAt(j).setFiltered(1);
                    }
                } else {
                    break; // more than the minimum distance, skip the comparison of the remaining peaks
                }
            }
        }

        // create new peak list
        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        for (int i = 0; i < nPeaks; i++) {
            if (peaks.elementAt(i).getFiltered() == 0) {
                peaksTemp.add(peaks.elementAt(i));
            }
        }

        return peaksTemp;
    }

    public static Vector<PEAK> combineNearbyPeaks(Vector<PEAK> peaks, int minimumDistance) {
        if (peaks.size() <= 1) {
            return peaks;
        }

        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        Collections.sort(peaks);

        int nPeaks = peaks.size();
        PEAK peak = peaks.elementAt(0);
        for (int i = 1; i < nPeaks; i++) {
            if ((peak.getChrom().compareTo(peaks.elementAt(i).getChrom()) == 0) && (peaks.elementAt(i).getStart() - peak.getEnd() < minimumDistance) && (peaks.elementAt(i).getnTags() == peak.getnTags())) {
                try {
                    peak.combine(peaks.elementAt(i));
                } catch (Exception e) {
                    System.out.println("Exception: " + e);
                }
            } else {
                peaksTemp.add(peak);
                peak = peaks.elementAt(i);
            }
        }

        return peaksTemp;
    }

    // all the regions with the coverage more than the minimum coverage are peak regions
    // bindingSiteCalling:  for peak region calling
    // bindingSiteCalling2: for peak summit calling
    // bindingSiteCalling3: for peak summit calling, with concise data structure
    // bindingSiteCalling4: for peak region calling, with concise data structure
    public static Vector<PEAK> bindingSiteCalling(Vector<HIT3> hit3s, int minimumCoverage, int minimumPeakDistance) {
        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        Collections.sort(hit3s);

        int nHit3s = hit3s.size();
        int insideApeak = 0; //1: inside a peak; 0: not inside a peak
        int coverage = 0;
        int maxCoverage = 0; //max coverage in a peak region
        String currentChrom = "";
        int hit3Start = 0;
        int hit3End = 0;
        int inSummit = 0;
        int summitStart = 0;
        int summitEnd = 0;
        for (int iHit3s = 0; iHit3s < nHit3s; iHit3s++) {
            HIT3 hit3 = hit3s.elementAt(iHit3s);
            if (currentChrom.compareTo(hit3.getChrom()) != 0) {
                // a new chromosome
                if (insideApeak == 1) {
                    PEAK.addPeak(peaksTemp, currentChrom, hit3Start, hit3End, summitStart, summitEnd, maxCoverage);
                }
                insideApeak = 0;
                coverage = 1;
                maxCoverage = 1;
                currentChrom = hit3.getChrom();
                hit3Start = 0;
                hit3End = 0;
                inSummit = 0;
                summitStart = 0;
                summitEnd = 0;
            } else {
                // the existing chromosome
                if (hit3.getStartEndLabel() >= 1) // start of a tag, increase the coverage
                {
                    coverage++;
                    if (insideApeak == 1) {
                        hit3End = hit3.getLoci();
                        if (coverage > maxCoverage) {
                            inSummit = 1;
                            maxCoverage = coverage;
                            summitStart = hit3.getLoci();
                            summitEnd = summitStart;
                        }
                    } else {
                        // a new peak region
                        if (coverage >= minimumCoverage) {
                            // identify a peak region
                            insideApeak = 1;
                            maxCoverage = coverage;
                            hit3Start = hit3.getLoci();
                            hit3End = hit3Start;
                            inSummit = 1;
                            summitStart = hit3Start;
                            summitEnd = hit3Start;
                        }
                    }

                } else if (hit3.getStartEndLabel() <= 0) { // end of a tag, decrease the coverage
                    coverage--;
                    if (insideApeak == 1) {
                        if (inSummit == 1) {
                            summitEnd = hit3.getLoci();
                            inSummit = 0;
                        }
                        hit3End = hit3.getLoci();
                        if (coverage < minimumCoverage) {
                            insideApeak = 0;
                            PEAK.addPeak(peaksTemp, currentChrom, hit3Start, hit3End, summitStart, summitEnd, maxCoverage);
                        }
                    }
                }
            }
        }
        if (insideApeak == 1) {
            PEAK.addPeak(peaksTemp, currentChrom, hit3Start, hit3End, summitStart, summitEnd, maxCoverage);
        }

        peaksTemp = PEAK.filterNearbyPeaks(peaksTemp, minimumPeakDistance);

        return peaksTemp;
    }

    // all the regions with the local maximum coverage are peak regions
    // bindingSiteCalling:  for peak region calling
    // bindingSiteCalling2: for peak summit calling
    // bindingSiteCalling3: for peak summit calling, with concise data structure
    // bindingSiteCalling4: for peak region calling, with concise data structure
    public static Vector<PEAK> bindingSiteCalling2(Vector<HIT3> hit3s, int minimumCoverage, int minimumPeakDistance) {
        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        Collections.sort(hit3s);

        int nHit3s = hit3s.size();
        int coverage = 0;
        String currentChrom = "";
        int inSummit = 0;
        int summitStart = 0;
        int summitEnd = 0;
        for (int iHit3s = 0; iHit3s < nHit3s; iHit3s++) {
            HIT3 hit3 = hit3s.elementAt(iHit3s);
            if (currentChrom.compareTo(hit3.getChrom()) != 0) {
                // a new chromosome
                if (inSummit == 1) {
                    PEAK.addPeak(peaksTemp, currentChrom, hit3.getLoci(), summitEnd, summitStart, summitEnd, coverage);
                }
                coverage = 1;
                currentChrom = hit3.getChrom();
                inSummit = 0;
                summitStart = 0;
                summitEnd = 0;
            } else {
                // the existing chromosome
                if (hit3.getStartEndLabel() == 1) // start of a tag, increase the coverage
                {
                    coverage++;
                    inSummit = 1;
                    summitStart = hit3.getLoci();
                    summitEnd = summitStart;
                } else { // end of a tag, decrease the coverage
                    if (inSummit == 1) {
                        summitEnd = hit3.getLoci();
                        if (coverage >= minimumCoverage) {
                            PEAK.addPeak(peaksTemp, currentChrom, summitStart, summitEnd, summitStart, summitEnd, coverage);
                        }
                        inSummit = 0;
                    }
                    coverage--;
                }
            }
        }
        if (inSummit == 1) {
            PEAK.addPeak(peaksTemp, currentChrom, summitStart, summitEnd, summitStart, summitEnd, coverage);
        }

        peaksTemp = PEAK.filterNearbyPeaks(peaksTemp, minimumPeakDistance);

        return peaksTemp;
    }

    // bindingSiteCalling:  for peak region calling
    // bindingSiteCalling2: for peak summit calling
    // bindingSiteCalling3: for peak summit calling, with concise data structure
    // bindingSiteCalling4: for peak region calling, with concise data structure
    public static Vector<PEAK> bindingSiteCalling3(Hashtable<String, Vector<Hit3Concise>> hit2Hash, int minimumCoverage, int minimumPeakDistance) {
        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        Set<String> chromNames = hit2Hash.keySet();
        if (chromNames != null) {
            /*
             Iterator<String> itr = chromNames.iterator();
             while (itr.hasNext()) {
             String chromName = itr.next();
             Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
             peaksTemp.addAll(PEAK.bindingSiteCalling3(chromName, hitsTemp, minimumCoverage, minimumPeakDistance));
             }*/
            Vector<String> chromNameV = new Vector<String>(chromNames);
            Collections.sort(chromNameV);
            for (int i = 0; i < chromNameV.size(); i++) {
                String chromName = chromNameV.elementAt(i);
                Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
                peaksTemp.addAll(PEAK.bindingSiteCalling3(chromName, hitsTemp, minimumCoverage, minimumPeakDistance));
            }
        }

        return peaksTemp;
    }

    // all the regions with the local maximum coverage are peak regions
    public static Vector<PEAK> bindingSiteCalling3(String chrom, Vector<Hit3Concise> hit3s, int minimumCoverage, int minimumPeakDistance) {
        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        Collections.sort(hit3s);

        int nHit3s = hit3s.size();
        int coverage = 0;
        int inSummit = 0;
        int summitStart = 0;
        int summitEnd = 0;
        for (int iHit3s = 0; iHit3s < nHit3s; iHit3s++) {
            Hit3Concise hit3 = hit3s.elementAt(iHit3s);
            if (hit3.getStartEndLabel() == 1) // start of a tag, increase the coverage
            {
                coverage++;
                inSummit = 1;
                summitStart = hit3.getLoci();
                summitEnd = summitStart;
            } else { // end of a tag, decrease the coverage
                if (inSummit == 1) {
                    summitEnd = hit3.getLoci();
                    if (coverage >= minimumCoverage) {
                        PEAK.addPeak(peaksTemp, chrom, summitStart, summitEnd, summitStart, summitEnd, coverage);
                    }
                    inSummit = 0;
                }
                coverage--;
            }
        }

        if (inSummit == 1) {
            PEAK.addPeak(peaksTemp, chrom, summitStart, summitEnd, summitStart, summitEnd, coverage);
        }

        peaksTemp = PEAK.filterNearbyPeaks(peaksTemp, minimumPeakDistance);

        return peaksTemp;
    }

    // bindingSiteCalling:  for peak region calling
    // bindingSiteCalling2: for peak summit calling
    // bindingSiteCalling3: for peak summit calling, with concise data structure
    // bindingSiteCalling4: for peak region calling, with concise data structure
    public static Vector<PEAK> bindingSiteCalling4(Hashtable<String, Vector<Hit3Concise>> hit2Hash, int minimumCoverage, int minimumPeakDistance) {
        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        Set<String> chromNames = hit2Hash.keySet();
        if (chromNames != null) {
            /*
             Iterator<String> itr = chromNames.iterator();
             while (itr.hasNext()) {
             String chromName = itr.next();
             Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
             peaksTemp.addAll(PEAK.bindingSiteCalling3(chromName, hitsTemp, minimumCoverage, minimumPeakDistance));
             }*/
            Vector<String> chromNameV = new Vector<String>(chromNames);
            Collections.sort(chromNameV);
            for (int i = 0; i < chromNameV.size(); i++) {
                String chromName = chromNameV.elementAt(i);
                Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
                System.out.println(chromName + ":\t" + hitsTemp.size() + " regions");
                peaksTemp.addAll(PEAK.bindingSiteCalling4(chromName, hitsTemp, minimumCoverage, minimumPeakDistance));
            }
        }

        return peaksTemp;
    }

    // all the regions with the coverage more than the minimum coverage are peak regions
    public static Vector<PEAK> bindingSiteCalling4(String chrom, Vector<Hit3Concise> hit3s, int minimumCoverage, int minimumPeakDistance) {
        Vector<PEAK> peaksTemp = new Vector<PEAK>();
        Collections.sort(hit3s);

        int nHit3s = hit3s.size();
        int coverage = 0;
        int maxCoverage = 0;
        int inEnrichedRegion = 0;
        int inSummit = 0;
        int start = 0;
        int end = 0;
        int summitStart = 0;
        int summitEnd = 0;
        for (int iHit3s = 0; iHit3s < nHit3s; iHit3s++) {
            Hit3Concise hit3 = hit3s.elementAt(iHit3s);
            if (hit3.getStartEndLabel() >= 1) // start of a tag, increase the coverage
            {
                coverage += hit3.getStartEndLabel();
                if ((inEnrichedRegion == 0) && (coverage >= minimumCoverage)) {
                    inEnrichedRegion = 1;
                    inSummit = 1;
                    start = hit3.getLoci();
                    end = start;
                    summitStart = start;
                    summitEnd = summitStart;
                    maxCoverage = coverage;
                } else if (inEnrichedRegion == 1) {
                    if (maxCoverage < coverage) {
                        maxCoverage = coverage;
                        summitStart = hit3.getLoci();
                        summitEnd = summitStart;
                        inSummit = 1;
                    }
                }
            } else if (hit3.getStartEndLabel() < 0) { // end of a tag, decrease the coverage
                if (inSummit == 1) {
                    summitEnd = hit3.getLoci();
                    inSummit = 0;
                }
                if (inEnrichedRegion == 1) {
                    end = hit3.getLoci();
                    if ((coverage >= minimumCoverage) && (coverage + hit3.getStartEndLabel() < minimumCoverage)) {
                        PEAK.addPeak(peaksTemp, chrom, start, end, summitStart, summitEnd, maxCoverage);
                        inEnrichedRegion = 0;
                    }
                }
                coverage += hit3.getStartEndLabel();
                if (coverage < 0) {
                    coverage = 0;
                }
            }
        }

        if (inEnrichedRegion == 1) {
            PEAK.addPeak(peaksTemp, chrom, start, end, summitStart, summitEnd, maxCoverage);
        }

        //peaksTemp = PEAK.filterNearbyPeaks(peaksTemp, minimumPeakDistance);
        return peaksTemp;
    }

    public static void addPeak(Vector<PEAK> peaksTemp, String chrom, int start, int end, int summitStart, int summitEnd, int nTags) {
        PEAK peak = new PEAK(chrom, start, end, summitStart, summitEnd, nTags);
        peaksTemp.add(peak);
    }

    public static ArrayList<PEAK> loadPeaks(String PeakFile) throws IOException {
        ArrayList<PEAK> peaks = new ArrayList<PEAK>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(PeakFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            // version before 2015-01-17
            //PEAK peak = new PEAK(fields[1], Integer.parseInt(fields[2]), Integer.parseInt(fields[3]), Integer.parseInt(fields[2]), Integer.parseInt(fields[3]), Integer.parseInt(fields[4]), Integer.parseInt(fields[0]));
            // version for 2015-01-17: without peak index
            int nTags = 0;
            if (fields.length >= 6) {
                nTags = Integer.parseInt(fields[5]);
            }
            PEAK peak = new PEAK(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), Integer.parseInt(fields[3]), Integer.parseInt(fields[4]), nTags);
            peaks.add(peak);
        }
        fileIn.close();
        return peaks;
    }

    public static ArrayList summedCoverage(Hashtable<String, Vector<Hit3Concise>> hit2Hash) {
        double summedCoverage = 0.0;
        Hashtable<Integer, Double> coverageDistribution = new Hashtable<Integer, Double>();

        Set<String> chromNames = hit2Hash.keySet();
        if (chromNames != null) {
            Vector<String> chromNameV = new Vector<String>(chromNames);
            Collections.sort(chromNameV);
            for (int i = 0; i < chromNameV.size(); i++) {
                String chromName = chromNameV.elementAt(i);
                Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
                summedCoverage += summedCoverage(hitsTemp, coverageDistribution);
            }
        }

        ArrayList coverageSummary = new ArrayList();
        coverageSummary.add(new Double(summedCoverage));
        coverageSummary.add(coverageDistribution);

        return coverageSummary;
    }

    // generate the summed coverage = sum(coverage * covered_length)
    public static double summedCoverage(Vector<Hit3Concise> hit3s, Hashtable<Integer, Double> coverageDistribution) {
        double summedCoverage = 0.0;
        int nHit3s = hit3s.size();

        if (nHit3s <= 1) {
            return summedCoverage;
        }

        Collections.sort(hit3s);

        Hit3Concise hit3 = hit3s.elementAt(0);
        int coverage = hit3.getStartEndLabel();
        int currentLoci = hit3.getLoci();
        for (int iHit3s = 1; iHit3s < nHit3s; iHit3s++) {
            hit3 = hit3s.elementAt(iHit3s);
            int coveredLength = hit3.getLoci() - currentLoci;
            summedCoverage += (coverage * coveredLength);
            Integer coverageInt = new Integer(coverage);
            Double totalCovered = coverageDistribution.get(coverageInt);
            if (totalCovered == null) {
                totalCovered = new Double(0);
            }
            totalCovered = new Double(totalCovered.doubleValue() + coveredLength);
            coverageDistribution.put(coverageInt, totalCovered);
            currentLoci = hit3.getLoci();
            coverage += hit3.getStartEndLabel();
            if (coverage < 0) {
                coverage = 0;
            }
        }

        return summedCoverage;
    }

    public static void save(Vector<PEAK> peaks, String peakFile, int mode) throws IOException {
        save(peaks.subList(0, peaks.size()), peakFile, mode);
    }

    public static void save(List<PEAK> peaks, String peakFile, int mode) throws IOException {
        if (peaks == null) {
            return;
        }
        PrintWriter peakFileOut = new PrintWriter(new BufferedWriter(new FileWriter(peakFile)));
        for (int iPeak = 0; iPeak < peaks.size(); iPeak++) {
            peakFileOut.println(peaks.get(iPeak).toString2(mode));
        }
        peakFileOut.close();
    }

    /**
     * @return the nTags_control
     */
    public int getnTags_control() {
        return nTags_control;
    }

    /**
     * @param nTags_control the nTags_control to set
     */
    public void setnTags_control(int nTags_control) {
        this.nTags_control = nTags_control;
    }

    /**
     * @return the id
     */
    public int getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(int id) {
        this.id = id;
    }

    /**
     * @return the iPets_intra
     */
    public int getiPets_intra() {
        return iPets_intra;
    }

    /**
     * @param iPets_intra the iPets_intra to set
     */
    public void setiPets_intra(int iPets_intra) {
        this.iPets_intra = iPets_intra;
    }

    /**
     * @return the iPets_inter
     */
    public int getiPets_inter() {
        return iPets_inter;
    }

    /**
     * @param iPets_inter the iPets_inter to set
     */
    public void setiPets_inter(int iPets_inter) {
        this.iPets_inter = iPets_inter;
    }
}
