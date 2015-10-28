/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

/**
 *
 * @author ligl
 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Hashtable;
import java.util.Vector;

public class ChromInfo {

    private String chromName;
    private int length;

    public  ChromInfo(String chromName, int length) {
        this.setChromName(chromName);
        this.setLength(length);
    }


    /**
     * @return the chromName
     */
    public String getChromName() {
        return chromName;
    }

    /**
     * @param chromName the chromName to set
     */
    public void setChromName(String chromName) {
        this.chromName = chromName;
    }

    /**
     * @return the length
     */
    public int getLength() {
        return length;
    }

    /**
     * @param length the length to set
     */
    public void setLength(int length) {
        this.length = length;
    }

    public static Vector<ChromInfo> readChromosomeInfo(String chromosomeInfoFile) throws IOException {
        Vector<ChromInfo> chromInfos = new Vector<ChromInfo>();
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(chromosomeInfoFile))));
        String line;
        int index = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String fields[] = line.split("\t");
            String chromName = fields[0];
            int length = Integer.parseInt(fields[1]);
            ChromInfo chromInfo = new ChromInfo(chromName, length);
            chromInfos.add(chromInfo);
            index++;
        }
        fileIn.close();

        return chromInfos;
    }

    public static Hashtable<String, Integer> getChromMapping(Vector<ChromInfo> chromInfos) {
        Hashtable<String, Integer> chromMapping = new Hashtable<String, Integer>();
        for (int iChrom = 0; iChrom < chromInfos.size(); iChrom++) {
            chromMapping.put(chromInfos.elementAt(iChrom).getChromName(), new Integer(iChrom));
        }
        return chromMapping;
    }
}
