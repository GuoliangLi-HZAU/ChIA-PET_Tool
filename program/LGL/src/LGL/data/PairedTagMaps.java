/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import java.util.Vector;

/**
 *
 * @author ligl
 */
public class PairedTagMaps {

    private String label = null;
    private TagMaps headMaps = new TagMaps();
    private TagMaps tailMaps = new TagMaps();
    private int count = 1;

    public PairedTagMaps() {
    }

    public PairedTagMaps(String line) {
        String fields[] = line.split("\t");
        if (fields.length > 0) {
            this.setLabel(fields[0]);
            String subFields[] = fields[0].split(":");
            if (subFields.length >= 2) {
                this.count = Integer.parseInt(subFields[1]);
            }
        }
        if (fields.length > 1) {
            this.headMaps.setTagSeq(fields[1]);
        }
        if (fields.length > 2) {
            this.tailMaps.setTagSeq(fields[2]);
        }
    }

    public void addHeadMaps(String line) {
        this.getHeadMaps().add(line);
    }

    public void addHeadMaps(String line, int mode) {
        this.getHeadMaps().add(line, mode);
    }

    public void addTailMaps(String line) {
        this.getTailMaps().add(line);
    }

    public void addTailMaps(String line, int mode) {
        this.getTailMaps().add(line, mode);
    }

    /**
     * @return the headMaps
     */
    public TagMaps getHeadMaps() {
        return headMaps;
    }

    /**
     * @param headMaps the headMaps to set
     */
    public void setHeadMaps(TagMaps headMaps) {
        this.headMaps = headMaps;
    }

    /**
     * @return the tailMaps
     */
    public TagMaps getTailMaps() {
        return tailMaps;
    }

    /**
     * @param tailMaps the tailMaps to set
     */
    public void setTailMaps(TagMaps tailMaps) {
        this.tailMaps = tailMaps;
    }

    /**
     * @return the label
     */
    public String getLabel() {
        return label;
    }

    /**
     * @param label the label to set
     */
    public void setLabel(String label) {
        this.label = label;
    }

    public static Vector<Vector<PairedTagMaps>> parse(Vector<PairedTagMaps> pairedTagMapsVector) {
        Vector<Vector<PairedTagMaps>> pairedTagMapsVectors = new Vector<Vector<PairedTagMaps>>();

        for (int i = 0; i < pairedTagMapsVector.size(); i++) {
        }

        return pairedTagMapsVectors;
    }

    // get the paired tagmaps wiht the best unique head map and the best unique tail map
    public PairedTagMaps parseBestUU() {
        PairedTagMaps pairedTagMapsNew = null;
        MAP bestUniqueHeadMap = this.getHeadMaps().bestUniqueMap();
        MAP bestUniqueTailMap = this.getTailMaps().bestUniqueMap();
        if ((bestUniqueHeadMap != null) && (bestUniqueTailMap != null)) {
            pairedTagMapsNew = new PairedTagMaps();
            pairedTagMapsNew.setLabel(this.getLabel());
            pairedTagMapsNew.getHeadMaps().setTagSeq(this.getHeadMaps().getTagSeq());
            pairedTagMapsNew.getTailMaps().setTagSeq(this.getTailMaps().getTagSeq());
            pairedTagMapsNew.getHeadMaps().setMaps(TagMaps.map2maps(bestUniqueHeadMap));
            pairedTagMapsNew.getTailMaps().setMaps(TagMaps.map2maps(bestUniqueTailMap));
        }

        return pairedTagMapsNew;
    }

    /**
     * @return the count
     */
    public int getCount() {
        return count;
    }

    /**
     * @param count the count to set
     */
    public void setCount(int count) {
        this.count = count;
    }
}
