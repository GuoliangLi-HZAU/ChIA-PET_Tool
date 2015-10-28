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
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class HIT implements Comparable {

    private String chrom;
    private int loci;
    private String annotation;

    public HIT() {
        this.setChrom("");
        this.setLoci(-1);
        this.setAnnotation("---");
    }

    public HIT(String chrom, int loci) {
        this.setChrom(chrom);
        this.setLoci(loci);
        this.setAnnotation("---");
    }

    public HIT(String chrom, int loci, String annotation) {
        this.setChrom(chrom);
        this.setLoci(loci);
        this.setAnnotation(annotation);
    }

    public boolean equals(HIT hit2) {
        boolean isEqual = ((getChrom().compareTo(hit2.getChrom()) == 0) && (getLoci() == hit2.getLoci()));
        return isEqual;
    }

    @Override
    public String toString() {
        return (this.getChrom() + "\t" + this.getLoci());
    }

    public String toString2() {
        return (this.getChrom() + "\t" + this.getLoci() + "\t" + this.getAnnotation());
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
     * @return the loci
     */
    public int getLoci() {
        return loci;
    }

    /**
     * @param loci the loci to set
     */
    public void setLoci(int loci) {
        this.loci = loci;
    }

    public int compareTo(Object anotherHit) throws ClassCastException {
        if (!(anotherHit instanceof HIT)) {
            throw new ClassCastException("A HIT object expected.");
        }
        int result = this.getChrom().compareTo(((HIT) anotherHit).getChrom());
        if (result == 0) // same chromosome
        {
            result = this.getLoci() - ((HIT) anotherHit).getLoci();
        }
        return result;
    }

    public static int calculateDistance(HIT hit1, HIT hit2) {
        int distance = Integer.MAX_VALUE;
        if (hit1.getChrom().compareTo(hit2.getChrom()) == 0) // the same chromosome
        {
            distance = Math.abs(hit1.getLoci() - hit2.getLoci());
        }
        return distance;
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

    // from the current HIT, extend to the upstream with "leftExtension", and extend to the downstream with "rightExtension"
    public static Vector<REGION> hit2Region(Vector<HIT> hits, int leftExtension, int rightExtension) {
        Vector<REGION> regions = new Vector<REGION>();
        for (int i = 0; i < hits.size(); i++) {
            REGION region = new REGION(hits.get(i).getChrom(), hits.get(i).getLoci() - leftExtension, hits.get(i).getLoci() + rightExtension);
            regions.add(region);
        }
        return regions;
    }

    public static Vector<HIT> load(String HitFile) throws IOException {
        Vector<HIT> hits = new Vector<HIT>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(HitFile))));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            HIT hit = new HIT(fields[0], Integer.parseInt(fields[1]));
            if (fields.length >= 3) {
                hit.setAnnotation(fields[2]);
            }
            hits.add(hit);
        }
        fileIn.close();

        return hits;
    }

    public static void save(String HitFile, Vector<HIT> hits) throws IOException {
        PrintWriter fileOut = new PrintWriter(new BufferedWriter(new FileWriter(HitFile)));

        for (int i = 0; i < hits.size(); i++) {
            fileOut.println(hits.elementAt(i).toString2());
        }

        fileOut.close();
    }

    public static Vector<HIT> filterHit(Vector<HIT> hits, int distance) {
        Vector<HIT> hitsTemp = new Vector<HIT>();
        if (hits.size() <= 0) {
            return hitsTemp;
        } else if (hits.size() == 1) {
            hitsTemp.add(hits.elementAt(0));
            return hitsTemp;
        }

        Collections.sort(hits);

        for (int i = 1; i < hits.size() - 1; i++) {
            if ((HIT.calculateDistance(hits.elementAt(i), hits.elementAt(i - 1)) < distance) || (HIT.calculateDistance(hits.elementAt(i), hits.elementAt(i + 1)) < distance)) {
                hitsTemp.add(hits.elementAt(i));
            }
        }
        // last element
        if (HIT.calculateDistance(hits.elementAt(hits.size() - 1), hits.elementAt(hits.size() - 2)) < distance) {
            hitsTemp.add(hits.elementAt(hits.size() - 1));
        }

        return hitsTemp;
    }
}
