/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

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
import java.util.Random;
import java.util.Vector;

public class REGION implements Comparable {

    private String chrom;
    private int start;
    private int end;
    private String annotation;

    public REGION() {
        this.setChrom("");
        this.setStart(-1);
        this.setEnd(-1);
        this.setAnnotation("---");
    }

    public REGION(String chrom, int start, int end) {
        this.setChrom(chrom);
        this.setStart(start);
        this.setEnd(end);
        this.setAnnotation("---");
    }

    public REGION(REGION anotherRegion) {
        this.setChrom(anotherRegion.getChrom());
        this.setStart(anotherRegion.getStart());
        this.setEnd(anotherRegion.getEnd());
        this.setAnnotation(anotherRegion.getAnnotation());
    }

    // the separator is TAB
    @Override
    public String toString() {
        return new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd());
    }

    public String toString(int mode) {
        return new String(this.getChrom() + ":" + this.getStart() + "-" + this.getEnd());
    }
    // the separator is specified by the user

    public String toString(String sep) {
        return new String(this.getChrom() + sep + this.getStart() + sep + this.getEnd());
    }

    // the separator is specified by the user
    // the output has the annotation information
    public String toString2(String sep) {
        return new String(this.getChrom() + sep + this.getStart() + sep + this.getEnd() + sep + this.getSpan() + sep + this.getAnnotation());
    }

    public void combine(REGION anotherRegion) throws Exception {
        if (this.getChrom().compareTo(anotherRegion.getChrom()) != 0) {
            throw new Exception("Two regions are not in the same chromosome");
        }

        int newStart = (this.getStart() < anotherRegion.getStart()) ? this.getStart() : anotherRegion.getStart();
        int newEnd = (this.getEnd() > anotherRegion.getEnd()) ? this.getEnd() : anotherRegion.getEnd();
        this.setStart(newStart);
        this.setEnd(newEnd);
    }

    public int compareTo(Object anotherRegion) throws ClassCastException {
        if (!(anotherRegion instanceof REGION)) {
            throw new ClassCastException("A REGION object expected.");
        }
        int result = this.getChrom().compareTo(((REGION) anotherRegion).getChrom());
        if (result == 0) // same chromosome
        {
            result = this.getStart() - ((REGION) anotherRegion).getStart();
            if (result == 0) // same chromosome and same start
            {
                result = this.getEnd() - ((REGION) anotherRegion).getEnd();
            }
        }
        return result;
    }

    /**
     * @return the chrom
     */
    public String getChrom() {
        return chrom;
    }

    /**
     * @param chrom the chrom to set
     */
    public void setChrom(String chrom) {
        this.chrom = chrom;
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

    public int getSpan() {
        return Math.abs(getStart() - getEnd());
    }

    public boolean overlap(REGION anotherRegion) {
        boolean overlapped = false;
        if (this.getChrom().compareTo(anotherRegion.getChrom()) == 0) {
            if (((this.getStart() >= anotherRegion.getStart()) && (this.getStart() <= anotherRegion.getEnd())) || ((anotherRegion.getStart() >= this.getStart()) && (anotherRegion.getStart() <= this.getEnd()))) {
                overlapped = true;
            }
        }

        return overlapped;
    }

    // Assume that there is no region overlapped in "regions"
    // Assume that the regions are sorted
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

        if (this.overlap(((Vector<REGION>) regions).elementAt(index)) == true) {
            overlapped = true;
            return overlapped;
        }

        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (this.overlap((REGION) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            }
        }
        for (int iRegion = index + 1; iRegion < regions.size(); iRegion++) {
            if (this.overlap((REGION) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            } else {
                break;
            }
        }

        return overlapped;
    }

    public boolean overlap(Vector regions, int maxSpan) {
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

        if (this.overlap(((Vector<REGION>) regions).elementAt(index)) == true) {
            overlapped = true;
            return overlapped;
        }

        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (this.overlap((REGION) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            } else if (this.getChrom().compareToIgnoreCase(((REGION) (regions.elementAt(iRegion))).getChrom()) != 0) {
                break; // two regions from different chromosomes
            } else if (this.getStart() - ((REGION) (regions.elementAt(iRegion))).getStart() > maxSpan) {
                break; // two regions from the same chromosome, and already far way from each other (and the regions with smaller indexes are further away)
            }
        }
        for (int iRegion = index + 1; iRegion < regions.size(); iRegion++) {
            if (this.overlap((REGION) (regions.elementAt(iRegion))) == true) {
                overlapped = true;
                return overlapped;
            } else {
                break;
            }
        }

        return overlapped;
    }

    public int overlappedSize(REGION anotherRegion) {
        int commonRegionSize = -1;
        if (this.getChrom().compareTo(anotherRegion.getChrom()) == 0) {
            int overlappedEnd = (this.getEnd() < anotherRegion.getEnd()) ? this.getEnd() : anotherRegion.getEnd();
            if ((this.getStart() >= anotherRegion.getStart()) && (this.getStart() < anotherRegion.getEnd())) {
                commonRegionSize = overlappedEnd - this.getStart();
            } else if ((anotherRegion.getStart() >= this.getStart()) && (anotherRegion.getStart() < this.getEnd())) {
                commonRegionSize = overlappedEnd - anotherRegion.getStart();
            }
        }

        return commonRegionSize;
    }

    public static int calculateCenterDistance(REGION region1, REGION region2) {
        int distance = Integer.MAX_VALUE;
        if (region1.getChrom().compareTo(region2.getChrom()) == 0) // the same chromosome
        {
            distance = Math.abs(region1.getStart() + region1.getEnd() - region2.getStart() - region2.getEnd()) / 2;
        }
        return distance;
    }

    public static int calculateBoundaryDistance(REGION region1, REGION region2) {
        int distance = Integer.MAX_VALUE;
        // if two regions are overlapped, the distance between them is 0
        if (region1.overlap(region2)) {
            distance = 0;
        } else if (region1.getChrom().compareTo(region2.getChrom()) == 0) // the same chromosome
        {
            // check the 4 combinations of the starts and ends from two regions
            // actually, only 2 combinations are valid after checking the overlap of the two regions
            distance = Math.abs(region1.getStart() - region2.getEnd());
            if (distance > Math.abs(region1.getEnd() - region2.getStart())) {
                distance = Math.abs(region1.getEnd() - region2.getStart());
            }
        }
        return distance;
    }

    // search the region list which may have overlapped regions itself
    // Assume the reference regions are sorted
    public static REGION searchOverlappedRegion(Vector regions, REGION regionToSearch) {
        int index = Collections.binarySearch(regions, regionToSearch);
        // the region is in the region list
        if (index >= 0) {
            return ((Vector<REGION>) regions).elementAt(index);
        }

        // the regionToSearch is not in the region list
        REGION region = null;
        for (int iRegion = 0; iRegion < regions.size(); iRegion++) {
            // exhaustive search, in case the regions are overlapped
            // if the regions are not overlapped, binary search can be faster
            region = ((Vector<REGION>) regions).elementAt(iRegion);
            if (region.overlap(regionToSearch) == true) {
                return region;
            }
        }

        return null;
    }

    // search the region list which donot have overlapped regions itself
    // Assume the reference regions are sorted
    public static REGION searchOverlappedRegion2(Vector regions, REGION regionToSearch) {
        int index = Collections.binarySearch(regions, regionToSearch);
        // the region is in the region list
        if (index >= 0) {
            return ((Vector<REGION>) regions).elementAt(index);
        }

        if (index < 0) {
            index = -index - 1;
        }
        if (index >= regions.size()) {
            index = regions.size() - 1;
        }

        if (regionToSearch.overlap(((Vector<REGION>) regions).elementAt(index)) == true) {
            return ((Vector<REGION>) regions).elementAt(index);
        }

        for (int iRegion = index - 1; iRegion >= 0; iRegion--) {
            if (regionToSearch.overlap((REGION) (regions.elementAt(iRegion))) == true) {
                return ((Vector<REGION>) regions).elementAt(iRegion);
            } else {
                break;
            }
        }
        for (int iRegion = index + 1; iRegion < regions.size(); iRegion++) {
            if (regionToSearch.overlap((REGION) (regions.elementAt(iRegion))) == true) {
                return ((Vector<REGION>) regions).elementAt(iRegion);
            } else {
                break;
            }
        }


        /*
         // the regionToSearch is not in the region list
         REGION region = null;
         for (int iRegion = -index - 1; iRegion >= 0; iRegion--) {
         // exhaustive search, in case the regions are overlapped
         // if the regions are not overlapped, binary search can be faster
         region = ((Vector<REGION>) regions).elementAt(iRegion);
         if (region.overlap(regionToSearch) == true) {
         return region;
         }
         if ((region.getEnd() < regionToSearch.getStart()) || (region.getStart() > regionToSearch.getEnd())) {
         return null;
         }
         }
         */
        return null;
    }

    public static Vector<REGION> load(String RegionFile) throws IOException {
        Vector<REGION> regions = new Vector<REGION>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(RegionFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            REGION region = new REGION(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
            if (fields.length >= 4) {
                region.setAnnotation(fields[3]);
            }
            if (region.getStart() > region.getEnd()) {
                System.out.println("start is larger than end: " + region.toString());
                region.setStart(Integer.parseInt(fields[2]));
                region.setEnd(Integer.parseInt(fields[1]));
                System.out.println("exchanged start and end");
            }
            regions.add(region);
        }
        fileIn.close();

        //Collections.sort(regions);
        return regions;
    }

    public static void save(String RegionFile, Vector<REGION> regions) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(RegionFile)));

        for (int i = 0; i < regions.size(); i++) {
            fileOut.println(regions.elementAt(i).toString());
        }

        fileOut.close();
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

    public static REGION randomRegion(Random rand, Vector<ChromInfo> chromInfos, long genomeLength, int extensionLength) {
        REGION region = new REGION();
        long genomeLoci = (long) (rand.nextDouble() * genomeLength);
        for (int i = 0; i < chromInfos.size(); i++) {
            if (genomeLoci < chromInfos.get(i).getLength()) {
                region.setChrom(chromInfos.get(i).getChromName());
                region.setStart((int) genomeLoci);
                region.setEnd((int) genomeLoci + extensionLength);
                break;
            } else {
                genomeLoci -= chromInfos.get(i).getLength();
            }
        }
        return region;
    }
}
