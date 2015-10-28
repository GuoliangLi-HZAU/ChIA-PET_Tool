/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

/**
 *
 * @author ligl
 */
public class TagMaps {
    // sample input
    //iTag    chrom   strand  loci     tagLength  nMisMatches   pos1  pos2
    //1       chr17   +       7976443  20         0
    //2       chr13   +       4293528  20         2             6     16

    private String tagSeq = null;
    // "maps" contains all the mapping information for one tag
    // each element in "maps" is the mapping information with a specific mismatch
    // for example, maps.elementAt(0) contains the mapping information for 0 mismatch
    // maps.elementAt(1) containts the mapping information for 1 mismatch
    // and so on.
    // The elements in maps.elementAt(0) is one of the mapping location with 0 mismatch
    private Vector<Vector<MAP>> maps = new Vector<Vector<MAP>>();

    public void add(String line) {
        String fields[] = line.split("\t");
        if (fields.length < 5) {
            return;
        }
        String chrom = fields[1];
        char strand = fields[2].charAt(0);
        int loci = Integer.parseInt(fields[3]);
        int nMisMatches = Integer.parseInt(fields[4]);
        int tagLength = 20;
        if (fields.length >= 6) {
            tagLength = Integer.parseInt(fields[5]);
        }
        MAP map = new MAP(chrom, loci, strand, nMisMatches, tagLength, null);
        if (this.getMaps().size() < nMisMatches + 1) {
            int currentSize = this.getMaps().size();
            for (int m = currentSize; m < nMisMatches + 1; m++) {
                this.getMaps().add(new Vector<MAP>());
            }
        }
        this.getMaps().elementAt(nMisMatches).add(map);
    }

    public void add(String line, int mode) {
        switch (mode) {
            case 1:
                add(line);
                break;
            case 2:
                add2(line);
                break;
            default:
                add(line);
                break;
        }
    }

    public void add2(String line) {
        String fields[] = line.split("\t");
        if (fields.length < 6) {
            return;
        }
        int nMisMatches = Integer.parseInt(fields[1]);
        char strand = fields[2].charAt(0);
        String chrom = fields[3];
        int loci = Integer.parseInt(fields[4]);
        int tagLength = Integer.parseInt(fields[5]);
        MAP map = new MAP(chrom, loci, strand, nMisMatches, tagLength, null);
        if (this.getMaps().size() < nMisMatches + 1) {
            int currentSize = this.getMaps().size();
            for (int m = currentSize; m < nMisMatches + 1; m++) {
                this.getMaps().add(new Vector<MAP>());
            }
        }
        this.getMaps().elementAt(nMisMatches).add(map);
    }

    public MAP bestUniqueMap() {
        MAP map = null;
        for (int i = 0; i < maps.size(); i++) {
            if (maps.elementAt(i).size() == 1) {
                map = maps.elementAt(i).elementAt(0);
                break;
            } else if (maps.elementAt(i).size() > 1) {
                break;
            }
        }
        return map;
    }

    public Vector<MAP> bestMaps() {
        Vector<MAP> mapVector = null;
        for (int i = 0; i < maps.size(); i++) {
            if (maps.elementAt(i).size() >= 1) {
                mapVector = maps.elementAt(i);
                break;
            }
        }
        return mapVector;
    }

    public static Vector<Vector<MAP>> map2maps(MAP map) {
        Vector<Vector<MAP>> tmpMaps = new Vector<Vector<MAP>>();
        int nMisMatches = map.getnMisMatches();
        for (int m = 0; m < nMisMatches + 1; m++) {
            tmpMaps.add(new Vector<MAP>());
        }
        tmpMaps.elementAt(nMisMatches).add(map);

        return tmpMaps;
    }

    public void output(String prefix, PrintWriter fileOut) throws IOException {
        Vector<MAP> mapVector = null;
        for (int i = 0; i < maps.size(); i++) {
            mapVector = maps.elementAt(i);
            for (int j = 0; j < mapVector.size(); j++) {
                fileOut.println(prefix + "\t" + mapVector.elementAt(j).toString());
            }
        }
    }

    /**
     * @return the tagSeq
     */
    public String getTagSeq() {
        return tagSeq;
    }

    /**
     * @param tagSeq the tagSeq to set
     */
    public void setTagSeq(String tagSeq) {
        this.tagSeq = tagSeq;
    }

    /**
     * @return the maps
     */
    public Vector<Vector<MAP>> getMaps() {
        return maps;
    }

    /**
     * @param maps the maps to set
     */
    public void setMaps(Vector<Vector<MAP>> maps) {
        this.maps = maps;
    }

    public int getNmaps() {
        int nMaps = 0;
        for (int i = 0; i < maps.size(); i++) {
            nMaps += maps.elementAt(i).size();
        }
        return nMaps;
    }
}
