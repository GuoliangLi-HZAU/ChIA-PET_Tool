/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;

public class MappingStatistics {

    String inputFile = null;
    String outputPrefix = null;
    int scorecutoff = 20;
    int category;
    int[][] statistics;
    String[] label;

    public MappingStatistics(String inputFile, String outputPrefix, String mappingScoreCutoff) throws IOException {
        this.inputFile = inputFile;
        this.outputPrefix = outputPrefix;
        this.scorecutoff = Integer.parseInt(mappingScoreCutoff);
    }

    public void estimate() throws IOException {
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(this.inputFile)));
        PrintWriter outStatistics = new PrintWriter(new FileOutputStream(new File(this.outputPrefix + ".mapping_statistics.txt")));
        PrintWriter outSAM = new PrintWriter(new FileOutputStream(new File(this.outputPrefix + ".bedpe.selected.temp.sam")));
        PairedEndReadInSAM pairedEndReadInSAM = PairedEndReadInSAM.load(fileIn);
        //initiate arrays
        statistics = new int[3][3];
        Arrays.fill(statistics[0], 0);
        Arrays.fill(statistics[1], 0);

        while ((pairedEndReadInSAM.getSamDataForOnePairedEndRead() != null)) {
            outSAM.print(pairedEndReadInSAM.getHeader());
            if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 0 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 0) {
                statistics[0][0]++;
            } else if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 1 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 0) {
                statistics[1][0]++;
            } else if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 2 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 0) {
                statistics[2][0]++;
            } else if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 0 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 1) {
                statistics[0][1]++;
            } else if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 1 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 1) {
                statistics[1][1]++;
                outSAM.println(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0] + "\n" + pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]);
            } else if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 2 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 1) {
                statistics[2][1]++;
            } else if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 0 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 2) {
                statistics[0][2]++;
            } else if (classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[0]) == 1 && classify(pairedEndReadInSAM.getSamDataForOnePairedEndRead()[1]) == 2) {
                statistics[1][2]++;
            } else {
                statistics[2][2]++;
            }
            pairedEndReadInSAM = PairedEndReadInSAM.load(fileIn);
        }
        fileIn.close();

        label = new String[3];
        label[0] = "Non-mappable";
        label[1] = "Uniquely-mapped";
        label[2] = "Others";
        outStatistics.println("\tNon-mappable\tUniquely-mapped\tOthers");
        for (int i = 0; i < 3; i++) {
            outStatistics.print(label[i] + "\t");
            for (int j = 0; j < 2; j++) {
                outStatistics.print(statistics[i][j] + "\t");
            }
            outStatistics.print(statistics[i][2]);
            outStatistics.println();
        }
        outStatistics.close();
        outSAM.close();
    }

    // category = 0: non-mapped reads, with the criteria: label - XT:A:N, or no label
    // category = 1: uniquely-mapped reads, with the criteria: labels - XT:A:U, X0:i:1, X1:i:0, and mapping quality score is larger than or equal to the mapping score cutoff
    // category = 2: others
    public int classify(String line) {
        String[] fields = line.split("\t");
        int category_temp = 0;
        if (fields.length <= 11) {
            category_temp = 0;
        } else if (fields[11].equalsIgnoreCase("XT:A:N")) {
            category_temp = 0;
        } else if (fields[11].equalsIgnoreCase("XT:A:U") && Integer.valueOf(fields[4]) >= scorecutoff && fields[15].equalsIgnoreCase("X0:i:1") && fields[16].equalsIgnoreCase("X1:i:0")) {
            category_temp = 1;
        } else {
            category_temp = 2;
        }
        return category_temp;
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 3) {
            MappingStatistics mappingStatistics = new MappingStatistics(args[0], args[1], args[2]);
            mappingStatistics.estimate();
        } else {
            System.out.println("Usage: java MappingStatistics <inputFile> <outputPrefix> <mappingScoreCutoff>");
            System.out.println("       <inputFile>: String, short read mapping data in SAM format");
            System.out.println("       <outputPrefix>: String, prefix for output statistics and SAM files");
            System.out.println("       <mappingScoreCutoff>: Integer, mapping score cutoff for filtering the low quality reads");
        }
    }
}

class PairedEndReadInSAM {

    private String[] samDataForOnePairedEndRead = null;
    private String Header = "";

    public PairedEndReadInSAM() {
    }

    public static PairedEndReadInSAM load(BufferedReader FileIn) throws IOException {
        String[] sam_temp = new String[2];
        String line;
        int nLines = 0;
        String header = "";
        // assume that every two lines are for one paired-end read in SAM format
        while ((line = FileIn.readLine()) != null) {
            if (!line.startsWith("@")) {
                sam_temp[nLines] = line;
                nLines++;
            } else {
                header += line + "\n";
            }
            if (nLines == 2) {
                break;
            }
        }
        PairedEndReadInSAM pairedEndReadInSAM = new PairedEndReadInSAM();
        if (nLines == 2) {
            pairedEndReadInSAM.setSamDataForOnePairedEndRead(sam_temp);
            pairedEndReadInSAM.setHeader(header);
        }
        return pairedEndReadInSAM;
    }

    /**
     * @return the line
     */
    public String[] getSamDataForOnePairedEndRead() {
        return samDataForOnePairedEndRead;
    }

    /**
     * @param samDataForOnePairedEndRead the line to set
     */
    public void setSamDataForOnePairedEndRead(String[] samDataForOnePairedEndRead) {
        this.samDataForOnePairedEndRead = samDataForOnePairedEndRead;
    }

    /**
     * @return the Header
     */
    public String getHeader() {
        return Header;
    }

    /**
     * @param Header the Header to set
     */
    public void setHeader(String Header) {
        this.Header = Header;
    }

}
