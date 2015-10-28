/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.dataConcise;

/**
 *
 * @author ligl
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class RegionConcise implements Comparable {

    private int start;
    private int end;
    private String annotation;

    public RegionConcise() {
        this.setStart(-1);
        this.setEnd(-1);
        this.setAnnotation("---");
    }

    public RegionConcise(int start, int end) {
        this.setStart(start);
        this.setEnd(end);
        this.setAnnotation("---");
    }

    public void combine(RegionConcise anotherRegionConcise) throws Exception {
        int newStart = (this.getStart() < anotherRegionConcise.getStart()) ? this.getStart() : anotherRegionConcise.getStart();
        int newEnd = (this.getEnd() > anotherRegionConcise.getEnd()) ? this.getEnd() : anotherRegionConcise.getEnd();
        this.setStart(newStart);
        this.setEnd(newEnd);
    }

    // combine the regions in the vector as a list of non-overlapped regions
    public static Vector<RegionConcise> combine(Vector<RegionConcise> regions) throws Exception {
        // If there is less than or equal to 1 region, no need combination and return the original vector
        if(regions.size() <= 1) {
            return regions;
        }

        Collections.sort(regions);
        Vector<RegionConcise> newRegions = new Vector<RegionConcise>();
        RegionConcise region = regions.elementAt(0);
        for (int i = 1; i < regions.size(); ++i) {
            if(region.overlap(regions.elementAt(i))) {
                region.combine(regions.elementAt(i));
            }
            else {
                newRegions.add(region);
                region = regions.elementAt(i);
            }
        }
        newRegions.add(region);
        return newRegions;
    }

    // calculate the coverage by the non-overlapped region list
    public static int coverage(Vector<RegionConcise> regions) {
        int coveredSize = 0;
        for(int i=0; i<regions.size(); ++i) {
            coveredSize += regions.elementAt(i).getSpan();
        }
        return coveredSize;
    }

    public int compareTo(Object anotherRegionConcise) throws ClassCastException {
        if (!(anotherRegionConcise instanceof RegionConcise)) {
            throw new ClassCastException("A RegionConcise object expected.");
        }
        int result = this.getStart() - ((RegionConcise) anotherRegionConcise).getStart();
        if (result == 0) // same start
        {
            result = this.getEnd() - ((RegionConcise) anotherRegionConcise).getEnd();
        }

        return result;
    }

    public int getSpan() {
        return Math.abs(this.getEnd() - this.getStart());
    }

    public boolean overlap(RegionConcise anotherRegionConcise) {
        boolean overlapped = false;

        if (((this.getStart() >= anotherRegionConcise.getStart()) && (this.getStart() < anotherRegionConcise.getEnd())) || ((anotherRegionConcise.getStart() >= this.getStart()) && (anotherRegionConcise.getStart() < this.getEnd()))) {
            overlapped = true;
        }

        return overlapped;
    }

    // Assume that there are overlapped regions inside the "regions"
    // Assume that "regions" is already sorted in ascending order
    public boolean overlap(Vector regions) {
        boolean overlapped = false;
        int index = Collections.binarySearch(regions, this);
        if (index >= 0) {
            overlapped = true;
            return overlapped;
        }

        if (index < 0) {
            index = -index - 1;
        }
        if (index >= regions.size()) {
            index = regions.size() - 1;
        }

        if (this.overlap(((Vector<RegionConcise>) regions).elementAt(index)) == true) {
            overlapped = true;
            return overlapped;
        }

        // exhaustive search, in case that there is an extreme region to cover the whole chromosome
        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (this.overlap((RegionConcise) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            }
        }
        for (int iRegion = index + 1; iRegion < regions.size(); iRegion++) {
            if (this.overlap((RegionConcise) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            } else {
                break;
            }
        }

        return overlapped;
    }

    // Assume that there is no region overlapped in "regions"
    // Assume that "regions" is already sorted in ascending order
    public boolean overlap2(Vector regions) {
        boolean overlapped = false;
        int index = Collections.binarySearch(regions, this);
        if (index >= 0) {
            overlapped = true;
            return overlapped;
        }

        if (index < 0) {
            index = -index - 1;
        }
        if (index >= regions.size()) {
            index = regions.size() - 1;
        }

        if (this.overlap(((Vector<RegionConcise>) regions).elementAt(index)) == true) {
            overlapped = true;
            return overlapped;
        }

        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (this.overlap((RegionConcise) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            } else {
                break;
            }
        }
        for (int iRegion = index + 1; iRegion < regions.size(); iRegion++) {
            if (this.overlap((RegionConcise) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            } else {
                break;
            }
        }

        return overlapped;
    }

    public int overlappedSize(RegionConcise anotherRegionConcise) {
        int commonRegionSize = -1;
        int overlappedEnd = (this.getEnd() < anotherRegionConcise.getEnd()) ? this.getEnd() : anotherRegionConcise.getEnd();
        if ((this.getStart() >= anotherRegionConcise.getStart()) && (this.getStart() < anotherRegionConcise.getEnd())) {
            commonRegionSize = overlappedEnd - this.getStart();
        } else if ((anotherRegionConcise.getStart() >= this.getStart()) && (anotherRegionConcise.getStart() < this.getEnd())) {
            commonRegionSize = overlappedEnd - anotherRegionConcise.getStart();
        }

        return commonRegionSize;
    }

    public RegionConcise overlappedRegion(RegionConcise anotherRegionConcise) {
        int overlappedStart = (this.getStart() > anotherRegionConcise.getStart()) ? this.getStart() : anotherRegionConcise.getStart();
        int overlappedEnd = (this.getEnd() < anotherRegionConcise.getEnd()) ? this.getEnd() : anotherRegionConcise.getEnd();

        // if two regions have no overlap, return null
        if (overlappedStart > overlappedEnd) {
            return null;
        }

        RegionConcise newRegion = new RegionConcise(overlappedStart, overlappedEnd);
        return newRegion;
    }

    // input: two regions
    // output: a region with the smaller start and the larger end: the maximum coverage,
    //         even if the two input regions are not overlapped
    public static RegionConcise maximumCoverageRegion(RegionConcise regionconcise1, RegionConcise regionConcise2) {
        int overlappedStart = (regionconcise1.getStart() < regionConcise2.getStart()) ? regionconcise1.getStart() : regionConcise2.getStart();
        int overlappedEnd = (regionconcise1.getEnd() > regionConcise2.getEnd()) ? regionconcise1.getEnd() : regionConcise2.getEnd();

        // if two regions have no overlap, return null
        if (overlappedStart > overlappedEnd) {
            return null;
        }

        RegionConcise newRegion = new RegionConcise(overlappedStart, overlappedEnd);
        return newRegion;
    }

    public static int calculateCenterDistance(RegionConcise region1, RegionConcise region2) {
        int distance = Integer.MAX_VALUE;
        distance = Math.abs(region1.getStart() + region1.getEnd() - region2.getStart() - region2.getEnd()) / 2;
        return distance;
    }

    // if distance < 0
    // it means that two regions are overlapped; the absolute value of the distance is the overlapped size
    public static int calculateBoundaryDistance(RegionConcise region1, RegionConcise region2) {
        int distance = Integer.MAX_VALUE;
        RegionConcise R1, R2;
        // to make sure R1.start < R2.start
        if(region1.getStart() < region2.getStart()) {
            R1 = region1;
            R2 = region2;
        }
        else {
            R1 = region2;
            R2 = region1;
        }
        if(R1.getEnd() < R2.getStart()) {
            distance = R2.getStart() - R1.getEnd(); // R1 is upstream of R2
        }
        else if(R1.getEnd() < R2.getEnd()) {
            distance = R2.getStart() - R1.getEnd(); // R1 is partially intersected with R2
        }
        else {
            distance = R2.getStart() - R2.getEnd(); // R2 is inside R1
        }
        return distance;
    }

    // search the region list which may have overlapped regions itself
    // Assume the reference regions are sorted
    public static RegionConcise searchOverlappedRegion(Vector regions, RegionConcise regionToSearch) {
        int index = Collections.binarySearch(regions, regionToSearch);
        // the region is in the region list
        if (index >= 0) {
            return ((Vector<RegionConcise>) regions).elementAt(index);
        }

        // the regionToSearch is not in the region list
        RegionConcise region = null;
        for (int iRegion = 0; iRegion < regions.size(); iRegion++) {
            // exhaustive search, in case the regions are overlapped
            // if the regions are not overlapped, binary search can be faster
            region = ((Vector<RegionConcise>) regions).elementAt(iRegion);
            if (region.overlap(regionToSearch) == true) {
                return region;
            }
        }

        return null;
    }

    // search the region list which donot have overlapped regions itself
    // Assume the reference regions are sorted
    public static RegionConcise searchOverlappedRegion2(Vector regions, RegionConcise regionToSearch) {
        int index = Collections.binarySearch(regions, regionToSearch);
        // the region is in the region list
        if (index >= 0) {
            return ((Vector<RegionConcise>) regions).elementAt(index);
        }

        if (index < 0) {
            index = -index - 1;
        }
        if (index >= regions.size()) {
            index = regions.size() - 1;
        }

        if (regionToSearch.overlap(((Vector<RegionConcise>) regions).elementAt(index)) == true) {
            return ((Vector<RegionConcise>) regions).elementAt(index);
        }

        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (regionToSearch.overlap((RegionConcise) (regions.elementAt(iRegion))) == true) {
                return ((Vector<RegionConcise>) regions).elementAt(iRegion);
            } else {
                break;
            }
        }
        for (int iRegion = index + 1; iRegion < regions.size(); iRegion++) {
            if (regionToSearch.overlap((RegionConcise) (regions.elementAt(iRegion))) == true) {
                return ((Vector<RegionConcise>) regions).elementAt(iRegion);
            } else {
                break;
            }
        }

        return null;
    }

    public static Hashtable<String, Vector<RegionConcise>> load(String RegionFile) throws IOException {
        Hashtable<String, Vector<RegionConcise>> region2Hash = new Hashtable<String, Vector<RegionConcise>>();

        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(RegionFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            // create the region3
            String fields[] = line.split("\t");
            String chromName = fields[0];
            RegionConcise regionConcise = new RegionConcise(Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
            if (fields.length >= 4) {
                regionConcise.setAnnotation(fields[3]);
            }
            if (regionConcise.getStart() > regionConcise.getEnd()) {
                System.out.println("start is larger than end: " + chromName + "\t" + regionConcise.getStart() + "\t" + regionConcise.getEnd());
                regionConcise.setStart(Integer.parseInt(fields[2]));
                regionConcise.setEnd(Integer.parseInt(fields[1]));
                System.out.println("exchanged start and end");
            }

            // add the new region3 to the hashtable
            Vector<RegionConcise> regionsTemp = region2Hash.get(chromName);
            if (regionsTemp == null) {
                regionsTemp = new Vector<RegionConcise>();
                region2Hash.put(chromName, regionsTemp);
            }
            regionsTemp.add(regionConcise);
        }
        fileIn.close();

        // sort the regions in each chromosome
        Set<String> chromNames = region2Hash.keySet();
        if (chromNames != null) {
            Iterator<String> itr = chromNames.iterator();
            while (itr.hasNext()) {
                String chromName = itr.next();
                Vector<RegionConcise> regionsTemp = (Vector<RegionConcise>) region2Hash.get(chromName);
                Collections.sort(regionsTemp);
                region2Hash.put(chromName, regionsTemp);
            }
        }

        return region2Hash;
    }

    public static void save(String RegionFile, Hashtable<String, Vector<RegionConcise>> region2Hash) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(RegionFile)));

        Set<String> chromNameSet = region2Hash.keySet();
        if (chromNameSet != null) {
            Vector<String> chromNames = new Vector<String>(chromNameSet);
            Collections.sort(chromNames);
            Iterator<String> itr = chromNames.iterator();
            while (itr.hasNext()) {
                String chromName = itr.next();
                Vector<RegionConcise> regionsTemp = (Vector<RegionConcise>) region2Hash.get(chromName);
                for (int i = 0; i < regionsTemp.size(); i++) {
                    fileOut.println(chromName + "\t" + regionsTemp.elementAt(i).getStart() + "\t" + regionsTemp.elementAt(i).getEnd());
                }
            }
        }

        fileOut.close();
    }

    public static Hashtable<String, Integer> getLongestLength(Hashtable<String, Vector<RegionConcise>> region2Hash) throws IOException {
            // get the length of the longest region in each chromosome
        Hashtable<String, Integer> maxRegionSize = new Hashtable<String, Integer>();
        Set<String> chromNameSet = region2Hash.keySet();
        if (chromNameSet != null) {
            Vector<String> chromNames = new Vector<String>(chromNameSet);
            Collections.sort(chromNames);
            Iterator<String> itr = chromNames.iterator();
            while (itr.hasNext()) {
                String chromName = itr.next();
                Vector<RegionConcise> regionsTemp = (Vector<RegionConcise>) region2Hash.get(chromName);
                int maxLength = -1;
                for (int i = 0; i < regionsTemp.size(); i++) {
                    int length = Math.abs(regionsTemp.elementAt(i).getEnd() - regionsTemp.elementAt(i).getStart());
                    if (maxLength < length) {
                        maxLength = length;
                    }
                    //System.out.println(chromName + "\t" + regionsTemp.elementAt(i).getStart() + "\t" + regionsTemp.elementAt(i).getEnd());
                }
                maxRegionSize.put(chromName, new Integer(maxLength));
            }
        }
        return maxRegionSize;
}

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * @return the annotation
     */
    public String getAnnotation() {
        return annotation;
    }

    /**
     * @param annotation the annotation to set
     */
    public void setAnnotation(String annotation) {
        this.annotation = annotation;
    }
}
