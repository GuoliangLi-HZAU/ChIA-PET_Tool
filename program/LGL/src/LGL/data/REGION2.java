/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import LGL.util.SeqUtil;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class REGION2 extends REGION {

    private char strand;

    public REGION2() {
        super();
        this.setStrand('+');
    }

    public REGION2(String chrom, int start, int end, char strand) {
        this.setChrom(chrom);
        this.setStart(start);
        this.setEnd(end);
        this.setStrand(strand);
    }

    public void combine(REGION2 anotherRegion) throws Exception {
        if (this.getChrom().compareTo(anotherRegion.getChrom()) != 0) {
            throw new Exception("Two region2s are not in the same chromosome: " + "\n" + this.toString() + "\n" + anotherRegion.toString());
        }
        if (this.getStrand() != anotherRegion.getStrand()) {
            throw new Exception("Two region2s are not in the same strand: " + "\n" + this.toString() + "\n" + anotherRegion.toString());
        }
        int newStart = (this.getStart() < anotherRegion.getStart()) ? this.getStart() : anotherRegion.getStart();
        int newEnd = (this.getEnd() > anotherRegion.getEnd()) ? this.getEnd() : anotherRegion.getEnd();
        this.setStart(newStart);
        this.setEnd(newEnd);
    }

    // cannot override the method load in REGION
    public static Vector<REGION2> loadRegion2(String Region2File) throws IOException {
        Vector<REGION2> region2s = new Vector<REGION2>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(Region2File))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            REGION2 region2 = new REGION2(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), fields[3].charAt(0));
            if (fields.length >= 5) {
                region2.setAnnotation(fields[4]);
            }
            if (region2.getStart() > region2.getEnd()) {
                System.out.println("start is larger than end: " + region2.toString());
                region2.setStart(Integer.parseInt(fields[2]));
                region2.setEnd(Integer.parseInt(fields[1]));
                System.out.println("exchanged start and end");
            }
            region2s.add(region2);
        }
        fileIn.close();

        //Collections.sort(region2s);
        return region2s;
    }

    public static void saveRegion2(String Region2File, Vector<REGION2> region2s) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(Region2File)));

        for (int i = 0; i < region2s.size(); i++) {
            fileOut.println(region2s.elementAt(i).toString());
        }

        fileOut.close();
    }

    @Override
    public String toString() {
        return new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getStrand());
    }

    @Override
    public int compareTo(Object anotherRegion) throws ClassCastException {
        if (!(anotherRegion instanceof REGION2)) {
            throw new ClassCastException("A REGION2 object expected.");
        }
        int result = this.getChrom().compareTo(((REGION2) anotherRegion).getChrom());
        if (result == 0) // same chromosome
        {
            if (this.getStrand() != ((REGION2) anotherRegion).getStrand()) {
                if (this.getStrand() == '+') {
                    result = -1;
                } else {
                    result = 1;
                }
            }
            if (result == 0) // same chromosome and same strand
            {
                result = this.getStart() - ((REGION2) anotherRegion).getStart();
                if (result == 0) // same chromosome, same strand, same start
                {
                    result = this.getEnd() - ((REGION2) anotherRegion).getEnd();
                }
            }
        }
        return result;
    }

    /**
     * @return the strand
     */
    public char getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(char strand) {
        this.strand = SeqUtil.getStrand(strand);
    }
}
