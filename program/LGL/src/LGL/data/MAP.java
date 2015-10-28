/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import LGL.util.SeqUtil;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class MAP {
    // sample input
    //iTag    chrom   strand  loci     tagLength  nMisMatches   pos1  pos2
    //1       chr17   +       7976443  20         0
    //2       chr13   +       4293528  20         2             6     16

    private String chrom;
    private char strand;
    private int loci;
    private int tagLength;
    private int nMisMatches; // number of mismatches in this hit
    private Vector<Integer> misMatchPositions;

    public MAP() {
        this.setChrom("");
        this.setStrand('+');
        this.setLoci(-1);
        this.setTagLength(Integer.MAX_VALUE);
        this.setnMisMatches(Integer.MAX_VALUE);
        this.setMisMatchPositions(null);
    }

    public MAP(String chrom, int loci, char strand) {
        this.setChrom(chrom);
        this.setLoci(loci);
        this.setStrand(strand);
        this.setTagLength(Integer.MAX_VALUE);
        this.setnMisMatches(Integer.MAX_VALUE);
        this.setMisMatchPositions(null);
    }

    public MAP(String chrom, int loci, char strand, int nMismatches) {
        this.setChrom(chrom);
        this.setLoci(loci);
        this.setStrand(strand);
        this.setTagLength(Integer.MAX_VALUE);
        this.setnMisMatches(nMismatches);
        this.setMisMatchPositions(null);
    }

    public MAP(String chrom, int loci, char strand, int nMismatches, int tagLength, Vector<Integer> misMatchPositions) {
        this.setChrom(chrom);
        this.setLoci(loci);
        this.setStrand(strand);
        this.setTagLength(tagLength);
        this.setnMisMatches(nMismatches);
        this.setMisMatchPositions(misMatchPositions);
    }

    public HIT2 toHIT2() {
        return (new HIT2(this.getChrom(), this.getLoci(), this.getStrand()));
    }

    public boolean equals(MAP map2) {
        boolean isEqual = ((this.getChrom().compareTo(map2.getChrom()) == 0) && (this.getLoci() == map2.getLoci()) && (this.getStrand() == map2.getStrand()));
        return isEqual;
    }

    @Override
    public String toString() {
        String str = (this.getChrom() + "\t" + this.getStrand() + "\t" + this.getLoci() + "\t" + this.getnMisMatches() + "\t" + this.getTagLength());
        if (misMatchPositions != null) {
            for (int i = 0; i < misMatchPositions.size(); i++) {
                str += ("\t" + misMatchPositions.elementAt(i).toString());
            }
        }
        return str;
    }

    // just for a short version of the string
    public String toString(int mode) {
        String str = (this.getChrom() + "\t" + this.getStrand() + "\t" + this.getLoci() + "\t" + this.getnMisMatches());
        return str;
    }

    public MAP reverse(int tagLength) {
        MAP map2 = new MAP();
        map2.setChrom(this.getChrom());
        map2.setnMisMatches(this.getnMisMatches());
        if (SeqUtil.isForwardStrand(this.getStrand())) {
            map2.setStrand('-');
            map2.setLoci(this.getLoci() + tagLength - 1);
        } else {
            map2.setStrand('+');
            map2.setLoci(this.getLoci());
        }
        return map2;
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
        if (!(anotherHit instanceof MAP)) {
            throw new ClassCastException("A MAP object expected.");
        }
        String anotherHitString = ((MAP) anotherHit).toString();
        return toString().compareTo(anotherHitString);
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

    /**
     * @return the tagLength
     */
    public int getTagLength() {
        return tagLength;
    }

    /**
     * @param tagLength the tagLength to set
     */
    public void setTagLength(int tagLength) {
        this.tagLength = tagLength;
    }

    /**
     * @return the nMisMatches
     */
    public int getnMisMatches() {
        return nMisMatches;
    }

    /**
     * @param nMisMatches the nMisMatches to set
     */
    public void setnMisMatches(int nMisMatches) {
        this.nMisMatches = nMisMatches;
    }

    /**
     * @return the misMatchPositions
     */
    public Vector<Integer> getMisMatchPositions() {
        return misMatchPositions;
    }

    /**
     * @param misMatchPositions the misMatchPositions to set
     */
    public void setMisMatchPositions(Vector<Integer> misMatchPositions) {
        this.misMatchPositions = misMatchPositions;
    }
}
