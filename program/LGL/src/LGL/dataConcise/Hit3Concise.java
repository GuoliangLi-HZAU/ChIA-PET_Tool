/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.dataConcise;

import LGL.util.SeqUtil;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class Hit3Concise implements Comparable {

    private int loci;
    private int startEndLabel = 0; //1: start of a tag; -1: end of a tag; 0: label to report the coverage; others: non-defined

    public Hit3Concise() {
        this.setLoci(-1);
        this.setStartEndLabel(0);
    }

    public Hit3Concise(int loci, int label) {
        this.setLoci(loci);
        this.setStartEndLabel(label);
    }

    public boolean equals(Hit3Concise hit2) {
        boolean isEqual = (this.getLoci() == hit2.getLoci());
        return isEqual;
    }

    public String toString2() {
        return (this.getLoci() + "\t" + this.getStartEndLabel());
    }

    public int compareTo(Object anotherHit) throws ClassCastException {
        if (!(anotherHit instanceof Hit3Concise)) {
            throw new ClassCastException("A Hit3Concise object expected.");
        }
        int result = this.getLoci() - ((Hit3Concise) anotherHit).getLoci();
        if (result == 0) {
            result = this.getStartEndLabel() - ((Hit3Concise) anotherHit).getStartEndLabel();
        }

        return result;
    }

    public static int calculateDistance(Hit3Concise hit1, Hit3Concise hit2) {
        return Math.abs(hit1.getLoci() - hit2.getLoci());
    }

    public static Hashtable<String, Vector<Hit3Concise>> load(String HitFile, int extensionLength) throws IOException {
        Calendar rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start loading Hit3Concise from " + HitFile);

        Hashtable<String, Vector<Hit3Concise>> hit2Hash = new Hashtable<String, Vector<Hit3Concise>>();

        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(HitFile))));
        String line;
        int nHITs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 0) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            String chromName = fields[0];
            int loci = Integer.parseInt(fields[1]);
            char strand = fields[2].charAt(0);
            Hit3Concise.addHit3s(hit2Hash, chromName, loci, strand, extensionLength);
            nHITs++;
            if (nHITs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nHITs / 1000000) + "M HITs read from " + HitFile);
            }
        }
        fileIn.close();
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] totally " + (nHITs / 1000000.0) + "M HITs read from " + HitFile);

        // sort the regions in each chromosome
        Set<String> chromNames = hit2Hash.keySet();
        if (chromNames != null) {
            Iterator<String> itr = chromNames.iterator();
            while (itr.hasNext()) {
                String chromName = itr.next();
                Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
                Collections.sort(hitsTemp);
                hit2Hash.put(chromName, hitsTemp);
            }
        }

        return hit2Hash;
    }

    public static void addHit3s(Hashtable<String, Vector<Hit3Concise>> hit2Hash, String chromName, int loci, char strand, int extensionLength) {
        Vector<Hit3Concise> hit3sTemp = hit2Hash.get(chromName);
        if (hit3sTemp == null) {
            hit3sTemp = new Vector<Hit3Concise>();
            hit2Hash.put(chromName, hit3sTemp);
        }

        Hit3Concise hit3Concise;
        if (SeqUtil.isForwardStrand(strand) == true) {
            hit3Concise = new Hit3Concise(loci, 1);
            hit3sTemp.add(hit3Concise);
            hit3Concise = new Hit3Concise(loci + extensionLength, -1);
            hit3sTemp.add(hit3Concise);
        } else {
            hit3Concise = new Hit3Concise(loci - extensionLength, 1);
            hit3sTemp.add(hit3Concise);
            hit3Concise = new Hit3Concise(loci, -1);
            hit3sTemp.add(hit3Concise);
        }
    }

    public static Hashtable<String, Vector<Hit3Concise>> loadFromRegions(String RegionFile, int extensionLength) throws IOException {
        Calendar rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start loading Hit3Concise from " + RegionFile);

        Hashtable<String, Vector<Hit3Concise>> hit2Hash = new Hashtable<String, Vector<Hit3Concise>>();

        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(RegionFile))));
        String line;
        int nHITs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 0) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            String chrom = fields[0];
            int start = Integer.parseInt(fields[1]);
            int end = Integer.parseInt(fields[2]);
            int nOccurrences = 1;
            if (fields.length >= 4) {
                nOccurrences = Integer.parseInt(fields[3]);
            }
            Hit3Concise.addHit3s(hit2Hash, chrom, start, end, extensionLength, nOccurrences);
            nHITs++;
            if (nHITs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nHITs / 1000000) + "M HITs read from " + RegionFile);
            }
        }
        fileIn.close();
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] totally " + (nHITs / 1000000.0) + "M HITs read from " + RegionFile);

        // sort the regions in each chromosome
        Set<String> chromNames = hit2Hash.keySet();
        if (chromNames != null) {
            Iterator<String> itr = chromNames.iterator();
            while (itr.hasNext()) {
                String chromName = itr.next();
                Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
                Collections.sort(hitsTemp);
                hit2Hash.put(chromName, hitsTemp);
            }
        }

        return hit2Hash;
    }

    public static void addHit3s(Hashtable<String, Vector<Hit3Concise>> hit2Hash, String chrom, int start, int end, int extensionLength, int nOccurrences) {
        Vector<Hit3Concise> hit3sTemp = hit2Hash.get(chrom);
        if (hit3sTemp == null) {
            hit3sTemp = new Vector<Hit3Concise>();
            hit2Hash.put(chrom, hit3sTemp);
        }

        Hit3Concise hit3Concise;
        hit3Concise = new Hit3Concise(start - extensionLength, nOccurrences);
        hit3sTemp.add(hit3Concise);
        hit3Concise = new Hit3Concise(end + extensionLength, -nOccurrences);
        hit3sTemp.add(hit3Concise);
    }

    public static int getNumberOfReads(Hashtable<String, Vector<Hit3Concise>> hit2Hash) {
        int nReads = 0;
        Set<String> chromNames = hit2Hash.keySet();
        if (chromNames != null) {
            Vector<String> chromNameV = new Vector<String>(chromNames);
            Collections.sort(chromNameV);
            for (int i = 0; i < chromNameV.size(); i++) {
                String chromName = chromNameV.elementAt(i);
                Vector<Hit3Concise> hitsTemp = (Vector<Hit3Concise>) hit2Hash.get(chromName);
                nReads += hitsTemp.size();
            }
        }

        return nReads;
    }

    /**
     * @return the start
     */
    public int getLoci() {
        return loci;
    }

    /**
     * @param start the start to set
     */
    public void setLoci(int loci) {
        this.loci = loci;
    }

    /**
     * @return the startEndLabel
     */
    public int getStartEndLabel() {
        return startEndLabel;
    }

    /**
     * @param startEndLabel the startEndLabel to set
     */
    public void setStartEndLabel(int startEndLabel) {
        this.startEndLabel = startEndLabel;
    }
}
