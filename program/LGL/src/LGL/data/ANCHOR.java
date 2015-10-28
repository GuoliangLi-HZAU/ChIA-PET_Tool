/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

/**
 *
 * @author ligl
 */
public class ANCHOR extends REGION {

    private int iPets_intra;
    private int iPets_inter;
    private int hitsInARegion; // as the background information

    public ANCHOR() {
        this.setChrom("");
        this.setStart(-1);
        this.setEnd(-1);
        this.setiPets_inter(0);
        this.setiPets_intra(0);
        this.setHitsInARegion(0);
    }

    public ANCHOR(String chrom, int start, int end) {
        this.setChrom(chrom);
        this.setStart(start);
        this.setEnd(end);
        this.setiPets_inter(0);
        this.setiPets_intra(0);
        this.setHitsInARegion(0);
    }

    public ANCHOR(REGION region) {
        this.setChrom(region.getChrom());
        this.setStart(region.getStart());
        this.setEnd(region.getEnd());
        this.setiPets_inter(0);
        this.setiPets_intra(0);
        this.setHitsInARegion(0);
    }

    public ANCHOR(HIT2 hit2, int extendLength) {
        this.setChrom(hit2.getChrom());
        this.setiPets_inter(0);
        this.setiPets_intra(0);
        this.setHitsInARegion(0);

        if (hit2.getStrand() == '+') {
            this.setStart(hit2.getLoci());
            this.setEnd(hit2.getLoci() + extendLength);
        } else {
            this.setStart(hit2.getLoci() - extendLength);
            this.setEnd(hit2.getLoci());
        }
    }

    @Override
    public String toString() {
        return new String(this.getChrom() + "\t" + this.getStart() + "\t" + this.getEnd() + "\t" + this.getiPets_intra() + "\t" + this.getiPets_inter() + "\t" + this.getHitsInARegion());
    }

    public static int calculateDistance(ANCHOR ANCHOR1, ANCHOR ANCHOR2) {
        int distance = Integer.MAX_VALUE;
        if (ANCHOR1.getChrom().compareTo(ANCHOR2.getChrom()) == 0) // the same chromosome
        {
            distance = Math.abs(ANCHOR1.getStart() + ANCHOR1.getEnd() - ANCHOR2.getStart() - ANCHOR2.getEnd()) / 2;
        }
        return distance;
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

    /**
     * @return the hitsInARegion
     */
    public int getHitsInARegion() {
        return hitsInARegion;
    }

    /**
     * @param hitsInARegion the hitsInARegion to set
     */
    public void setHitsInARegion(int hitsInARegion) {
        this.hitsInARegion = hitsInARegion;
    }
}
