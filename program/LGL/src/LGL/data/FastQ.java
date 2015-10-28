/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.data;

import java.io.BufferedReader;
import java.io.IOException;

/**
 *
 * @author ligl
 */
public class FastQ {

    private String[] fastq = null;

    public FastQ() {

    }

    public static FastQ load(BufferedReader fastqFileIn1) throws IOException {
        String[] fastq_temp = new String[4];
        String line;
        int nLines = 0;
        if ((line = fastqFileIn1.readLine()) != null) {
            fastq_temp[0] = new String(line);
            nLines++;
        }
        if ((line = fastqFileIn1.readLine()) != null) {
            fastq_temp[1] = new String(line);
            nLines++;
        }
        if ((line = fastqFileIn1.readLine()) != null) {
            fastq_temp[2] = new String(line);
            nLines++;
        }
        if ((line = fastqFileIn1.readLine()) != null) {
            fastq_temp[3] = new String(line);
            nLines++;
        }

        FastQ fastQ_1 = new FastQ();
        if (nLines == 4) {
            fastQ_1.setFastq(fastq_temp);
        }
        return fastQ_1;
    }

    /**
     * @return the fastq
     */
    public String[] getFastq() {
        return fastq;
    }

    /**
     * @param fastq the fastq to set
     */
    public void setFastq(String[] fastq) {
        this.fastq = fastq;
    }
}
