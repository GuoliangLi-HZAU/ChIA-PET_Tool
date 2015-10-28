/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.chiapet;

import LGL.data.HIT3;
import LGL.data.PEAK;
import LGL.util.SeqUtil;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Calendar;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class BindingSitesFromPETs {

    int debugLevel = 3;
    Calendar rightNow = Calendar.getInstance();
    int extensionLength = 500;
    int selfLigationCutoff = 8000;
    int minimumCoverage = 3;
    int peakMode = 1; // 1: the peak is a region;   2: the peak is a local summit
    int minimumPeakDistance = 200;
    Vector<HIT3> hit3s = null;
    Vector<PEAK> peaks = null;

    public BindingSitesFromPETs(String inFile, String bindingSiteFile, int extensionLength, int selfLigationCutoff, int minimumCoverage, int peakMode, int minimumPeakDistance) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start BindingSites ... ");

        this.extensionLength = extensionLength;
        this.minimumCoverage = minimumCoverage;
        this.selfLigationCutoff = selfLigationCutoff;
        this.minimumPeakDistance = minimumPeakDistance;
        this.peakMode = peakMode;
        System.out.println("this.peakMode = " + this.peakMode);
        hit3s = load(inFile);
        if (this.peakMode == 1) {
            System.out.println("IN peakMode 1");
            peaks = PEAK.bindingSiteCalling(hit3s, minimumCoverage, minimumPeakDistance);
        } else {
            System.out.println("IN peakMode 2");
            peaks = PEAK.bindingSiteCalling2(hit3s, minimumCoverage, minimumPeakDistance);
        }
        System.out.println("Number of peaks: " + peaks.size());

        PrintWriter BSfileOut = new PrintWriter(new BufferedWriter(new FileWriter(bindingSiteFile)));
        for (int iPeak = 0; iPeak < peaks.size(); iPeak++) {
            BSfileOut.println(peaks.elementAt(iPeak).toString());
        }
        BSfileOut.close();
    }

    Vector<HIT3> load(String inFile) throws IOException {
        Vector<HIT3> hit3sTemp = new Vector<HIT3>();
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] " + "start loading PETs from " + inFile);

        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inFile))));

        String line;
        int nPETs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String[] fields = line.split("\t");
            String chrom1 = fields[0];
            int loci1 = Integer.parseInt(fields[1]);
            char strand1 = fields[2].charAt(0);
            String chrom2 = fields[3];
            int loci2 = Integer.parseInt(fields[4]);
            char strand2 = fields[5].charAt(0);

            if (chrom1.compareTo(chrom2) != 0) {
                // inter-ligations - different chromosomes
                addHit3s(hit3sTemp, chrom1, loci1, strand1);
                addHit3s(hit3sTemp, chrom2, loci2, strand2);
            } else if (Math.abs(loci1 - loci2) > selfLigationCutoff) {
                // inter-ligations - same chromosome, long distance
                addHit3s(hit3sTemp, chrom1, loci1, strand1);
                addHit3s(hit3sTemp, chrom2, loci2, strand2);
            } else if (strand1 == strand2) {
                // other pets
                // same chromosome, short distance, same strands
                addHit3s(hit3sTemp, chrom1, loci1, strand1);
                addHit3s(hit3sTemp, chrom2, loci2, strand2);
            } else if (SeqUtil.isForwardStrand(strand1) == true) {
                // same chromosome, short distance, different strands
                // read1 from plus strand, read2 from minus strand,
                if (loci1 < loci2) {
                    // other pets
                    // same chromosome, short distance, different strands
                    // read1 from plus strand, read2 from minus strand,
                    // read1 loci smaller than read2 loci
                    addHit3s(hit3sTemp, chrom1, loci1, strand1);
                    addHit3s(hit3sTemp, chrom2, loci2, strand2);
                } else {
                    // self ligation
                    // same chromosome, short distance, different strands
                    // read1 from plus strand, read2 from minus strand,
                    // read1 loci larger than or equal to read2 loci
                    HIT3 hit = new HIT3(chrom2, loci2, 1); // start of a tag
                    hit3sTemp.add(hit);
                    hit = new HIT3(chrom1, loci1, 0); // end of a tag
                    hit3sTemp.add(hit);
                }

            } else {
                // same chromosome, short distance, different strands
                // read1 from minus strand, read2 from plus strand,
                if (loci1 > loci2) {
                    // other pets
                    // same chromosome, short distance, different strands
                    // read1 from minus strand, read2 from plus strand,
                    // read1 loci larger than read2 loci
                    addHit3s(hit3sTemp, chrom1, loci1, strand1);
                    addHit3s(hit3sTemp, chrom2, loci2, strand2);
                } else {
                    // self ligation
                    // same chromosome, short distance, different strands
                    // read1 from minus strand, read2 from plus strand,
                    // read1 loci smaller than or equal to read2 loci
                    HIT3 hit = new HIT3(chrom1, loci1, 1); // start of a tag
                    hit3sTemp.add(hit);
                    hit = new HIT3(chrom2, loci2, 0); // end of a tag
                    hit3sTemp.add(hit);
                }
            }

            nPETs++;
            if (nPETs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000) + "M PETs read from " + inFile);
            }
        }
        fileIn.close();
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] totally " + (nPETs / 1000000) + "M PETs read from " + inFile);

        return hit3sTemp;
    }

    void addHit3s(Vector<HIT3> hit3sTemp, String chrom, int loci, char strand) {
        HIT3 hit;
        if (SeqUtil.isForwardStrand(strand) == true) {
            hit = new HIT3(chrom, loci - extensionLength, 1); // start of a tag
            hit3sTemp.add(hit);
            hit = new HIT3(chrom, loci, 0); // end of a tag
            hit3sTemp.add(hit);
        } else {
            hit = new HIT3(chrom, loci, 1); // start of a tag
            hit3sTemp.add(hit);
            hit = new HIT3(chrom, loci + extensionLength, 0); // end of a tag
            hit3sTemp.add(hit);
        }
    }

    // java -cp /data1/tmp/LGL.jar LGL.pipeline.BindingSitesFromPETs IHH025_1r78.linker_a.1.link.multiple.pet IHH025_1r78.linker_a.1.link.multiple.peaks.txt 500 4500 5
    public static void main(String[] args) throws IOException {
        if (args.length == 7) {
            new BindingSitesFromPETs(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]), Integer.parseInt(args[5]), Integer.parseInt(args[6]));
        } else {
            System.out.println("Usage: java BindingSitesFromPETs <pet_file_before_parse> <binding_site_file> <extensionLength> <selfLigationCutoff> <minimumCoverage> <peakMode> <minimumPeakDistance>");
            System.exit(1);
        }
    }
}
