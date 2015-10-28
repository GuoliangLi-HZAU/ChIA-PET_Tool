/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

/**
 *
 * @author ligl
 */
public class HIT3 extends HIT {

    private int startEndLabel = 0;

    public HIT3(String chrom, int loci, int label) {
        super(chrom, loci);
        this.setStartEndLabel(label);
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
