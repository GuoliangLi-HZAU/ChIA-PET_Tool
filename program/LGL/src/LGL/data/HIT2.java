/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import LGL.util.SeqUtil;

/**
 *
 * @author ligl
 */
public class HIT2 extends HIT {

    private char strand;

    public HIT2() {
        super();
        this.setStrand('+');
    }

    public HIT2(String chrom, int loci, char strand) {
        super(chrom, loci);
        this.setStrand(strand);
    }

    public boolean equals(HIT2 hit2) {
        boolean isEqual = ((this.getChrom().compareTo(hit2.getChrom()) == 0) && (this.getLoci() == hit2.getLoci()) && (this.getStrand() == hit2.getStrand()));
        return isEqual;
    }

    @Override
    public String toString() {
        return (this.getChrom() + "\t" + this.getLoci() + "\t" + this.getStrand());
    }

    public HIT2 reverse(int tagLength) {
        HIT2 hit2 = new HIT2();
        hit2.setChrom(this.getChrom());
        if (this.getStrand() == '+') {
            hit2.setStrand('-');
            hit2.setLoci(this.getLoci() + tagLength - 1);
        } else {
            hit2.setStrand('+');
            hit2.setLoci(this.getLoci());
        }
        return hit2;
    }

    @Override
    public int compareTo(Object anotherHit2) throws ClassCastException {
        if (!(anotherHit2 instanceof HIT2)) {
            throw new ClassCastException("A HIT2 object expected.");
        }
        int result = this.getChrom().compareTo(((HIT2) anotherHit2).getChrom());
        if (result == 0) // same chromosome
        {
            result = this.getLoci() - ((HIT2) anotherHit2).getLoci();
            if (result == 0) // same chromosome and same loci
            {
                if (this.getStrand() != ((HIT2) anotherHit2).getStrand()) { // same chromosome, different strands
                    if (this.getStrand() == '+') {
                        result = 1;
                    } else {
                        result = -1;
                    }
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

    public static int calculateDistance(HIT2 hit1, HIT2 hit2) {
        int distance = Integer.MAX_VALUE;
        //  same chromosome and same strand
        if ((hit1.getChrom().compareTo(hit2.getChrom()) == 0) && (hit1.getStrand() == hit2.getStrand())) {
            distance = Math.abs(hit1.getLoci() - hit2.getLoci());
        }
        return distance;
    }

    public MAP toMAP() {
        return (new MAP(this.getChrom(), this.getLoci(), this.getStrand(), 0, 20, null));
    }
}
