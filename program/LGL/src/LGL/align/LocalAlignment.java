/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.align;

import java.io.*;
import java.util.*;

public class LocalAlignment {

    private int debugLevel = 2;
    private int[][] scoreMatrix;
    private int maxScore;
    private int maxI; // at row maxI, the score is maximum
    private int maxJ; // at column maxJ, the score is maximum
    private int minI; // at row minI, the maximum score starts from here
    private int minJ; // at column minJ, the maximum score starts from here
    private int m;
    private int n;
    private int MatchScore = 1;
    private int MismatchScore = -1;
    private int IndelScore = -1;
    private String str1;
    private String str2;
    private int nMatches;
    private int nMismatches;
    private int nIndels;
    // ****** assumptions ******
    // the insertions and deletions are represented with '-' in the aligned strings
    private StringBuilder alignedStr1 = new StringBuilder();
    private StringBuilder alignedStr2 = new StringBuilder();
    private StringBuilder alignedStatus = new StringBuilder();

    public LocalAlignment(int m, int n) {
        this.m = m;
        this.n = n;
        scoreMatrix = new int[m + 1][n + 1];
        init();
    }

    public void init() {
        initMatrix(0);
        setMaxScore(0);
        setMaxI(0);
        setMaxJ(0);
        setMinI(0);
        setMinJ(0);
        this.setStr1("");
        this.setStr2("");
        this.setAlignedStatus("");
        this.setAlignedStr1("");
        this.setAlignedStr2("");
    }

    public void initMatrix(int score) {
        for (int i = 0; i <= m; i++) {
            Arrays.fill(scoreMatrix[i], score);
        }
    }

    public void outputScoreMatrix() {
        int index1, index2;
        int length1 = this.getStr1().length();
        int length2 = this.getStr2().length();
        System.out.println("length1 = " + length1);
        System.out.println("length2 = " + length2);

        System.out.print("i = 0:");
        for (index2 = 0; index2 <= length1; index2++) {
            System.out.print("\t" + index2);
        }
        System.out.println();
        for (index1 = 0; index1 <= length1; index1++) {
            System.out.print("i = " + index1 + ":");
            for (index2 = 0; index2 <= length2; index2++) {
                System.out.print("\t" + scoreMatrix[index1][index2]);
            }
            System.out.println();
        }
    }

    public void align(String str1, String str2) {
        int length1 = str1.length();
        int length2 = str2.length();
        // adjust the dimension of the score matrix
        boolean dimensionChanged = false;
        if (m < length1) {
            m = length1;
            dimensionChanged = true;
        }
        if (n < length2) {
            n = length2;
            dimensionChanged = true;
        }
        if (dimensionChanged == true) {
            scoreMatrix = new int[m + 1][n + 1];
        }
        init();
        this.setStr1(str1);
        this.setStr2(str2);

        // alignment to generate the maximum score
        int index1, index2;
        for (index1 = 0; index1 < length1; index1++) {
            // the start position in the read is 17, not 0, in order to speed up the linker alignment
            // the setting of 17 depends on the ChIA-PET read construct structure
            for (index2 = 17; index2 < length2; index2++) {
                if (str1.charAt(index1) == str2.charAt(index2)) {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1][index2] + MatchScore;
                } else {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1][index2] + MismatchScore;
                }
                if (scoreMatrix[index1 + 1][index2 + 1] < scoreMatrix[index1 + 1][index2] + IndelScore) {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1 + 1][index2] + IndelScore;
                }
                if (scoreMatrix[index1 + 1][index2 + 1] < scoreMatrix[index1][index2 + 1] + IndelScore) {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1][index2 + 1] + IndelScore;
                }
                if (scoreMatrix[index1 + 1][index2 + 1] < 0) {
                    scoreMatrix[index1 + 1][index2 + 1] = 0;
                }
                if (maxScore < scoreMatrix[index1 + 1][index2 + 1]) {
                    setMaxScore(scoreMatrix[index1 + 1][index2 + 1]);
                    setMaxI(index1 + 1);
                    setMaxJ(index2 + 1);
                }
            }
        }

        // trace back to get the maximum aligned sub-sequences
        alignedStr1.delete(0, alignedStr1.length());
        alignedStr2.delete(0, alignedStr2.length());
        alignedStatus.delete(0, alignedStatus.length());
        index1 = maxI - 1;
        index2 = maxJ - 1;
        setnMatches(0);
        setnMismatches(0);
        setnIndels(0);

        int scoreGreaterThan0 = 1;
        while ((index1 >= 0) || (index2 >= 0)) {
            if (index1 < 0) { // deletion at the beginning of str1, insertion in str2
                alignedStr1.append('-');
                alignedStr2.append(str2.charAt(index2));
                alignedStatus.append(' ');
                index2--;
                if (scoreGreaterThan0 == 1) {
                    setnIndels(getnIndels() + 1);
                }
            } else if (index2 < 0) { // insertion in str1, deletion at the beginning of str2
                alignedStr1.append(str1.charAt(index1));
                alignedStr2.append('-');
                alignedStatus.append(' ');
                index1--;
                if (scoreGreaterThan0 == 1) {
                    setnIndels(getnIndels() + 1);
                }
            } else if ((scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1][index2] + MatchScore) && (str1.charAt(index1) == str2.charAt(index2))) {
                // match from both strs
                alignedStr1.append(str1.charAt(index1));
                alignedStr2.append(str2.charAt(index2));
                alignedStatus.append('|');
                index1--;
                index2--;
                if (scoreGreaterThan0 == 1) {
                    setnMatches(getnMatches() + 1);
                }
            } else if (scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1][index2] + MismatchScore) {
                // mismatch from both strs
                alignedStr1.append(str1.charAt(index1));
                alignedStr2.append(str2.charAt(index2));
                alignedStatus.append('X');
                index1--;
                index2--;
                if (scoreGreaterThan0 == 1) {
                    setnMismatches(getnMismatches() + 1);
                }
            } else if (scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1 + 1][index2] + IndelScore) {
                // deletion in str1, insertion in str2
                alignedStr1.append('-');
                alignedStr2.append(str2.charAt(index2));
                alignedStatus.append(' ');
                index2--;
                if (scoreGreaterThan0 == 1) {
                    setnIndels(getnIndels() + 1);
                }
            } else {
                // insertion in str1, deletion in str2
                //if (scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1][index2 + 1] + IndelScore)
                alignedStr1.append(str1.charAt(index1));
                alignedStr2.append('-');
                alignedStatus.append(' ');
                index1--;
                if (scoreGreaterThan0 == 1) {
                    setnIndels(getnIndels() + 1);
                }
            }
            if (scoreGreaterThan0 == 1) {
                if (scoreMatrix[index1 + 1][index2 + 1] <= 0) {
                    minI = index1 + 1;
                    minJ = index2 + 1;
                    scoreGreaterThan0 = 0;
                }
            }
            if (debugLevel > 2) {
                System.out.println("index1 = " + (index1 + 1) + "; index2 = " + (index2 + 1));
            }
        }
        alignedStr1.reverse();
        alignedStr2.reverse();
        alignedStatus.reverse();
        if (debugLevel > 2) {
            System.out.println("str1: " + str1);
            System.out.println("str2: " + str2);
            System.out.println("aligned score: " + maxScore);
            System.out.println("aligned str1      : " + getAlignedStr1());
            System.out.println("aligned strStatus : " + getAlignedStatus());
            System.out.println("aligned str2      : " + getAlignedStr2());
        }
    }

    public static void main(String[] args) throws IOException {
        String str1 = "GTTGGATAAGATATGGCGG";
        String str2 = "ACTGCTGCTGCCAACACTACAGTTGGATAAGATATC";
        LocalAlignment localAligner = new LocalAlignment(str1.length(), str2.length());
        localAligner.align(str1, str2);
        /*  if (args.length == 2) {
         LocalAlignment localAligner = new LocalAlignment(args[0].length(), args[1].length());
         //localAligner.align("AACCGGTT", "ACCGTATT");
         localAligner.setDebugLevel(2);
         localAligner.align(args[0], args[1]);

         } else {
         System.out.println("Usage: java LocalAlignment <sequence 1> <sequence 2>");
         System.exit(0);
         }
         */
    }

    /**
     * @return the scoreMatrix
     */
    public int[][] getScoreMatrix() {
        return scoreMatrix;
    }

    /**
     * @return the score
     */
    public int getMaxScore() {
        return maxScore;
    }

    /**
     * @param score the score to set
     */
    public void setMaxScore(int score) {
        this.maxScore = score;
    }

    /**
     * @return the maxI
     */
    public int getMaxI() {
        return maxI;
    }

    /**
     * @param maxI the maxI to set
     */
    public void setMaxI(int maxI) {
        this.maxI = maxI;
    }

    /**
     * @return the maxJ
     */
    public int getMaxJ() {
        return maxJ;
    }

    /**
     * @param maxJ the maxJ to set
     */
    public void setMaxJ(int maxJ) {
        this.maxJ = maxJ;
    }

    /**
     * @return the m
     */
    public int getM() {
        return m;
    }

    /**
     * @param m the m to set
     */
    public void setM(int m) {
        this.m = m;
    }

    /**
     * @return the n
     */
    public int getN() {
        return n;
    }

    /**
     * @param n the n to set
     */
    public void setN(int n) {
        this.n = n;
    }

    /**
     * @return the MatchScore
     */
    public int getMatchScore() {
        return MatchScore;
    }

    /**
     * @param MatchScore the MatchScore to set
     */
    public void setMatchScore(int MatchScore) {
        this.MatchScore = MatchScore;
    }

    /**
     * @return the MismatchScore
     */
    public int getMismatchScore() {
        return MismatchScore;
    }

    /**
     * @param MismatchScore the MismatchScore to set
     */
    public void setMismatchScore(int MismatchScore) {
        this.MismatchScore = MismatchScore;
    }

    /**
     * @return the IndelScore
     */
    public int getIndelScore() {
        return IndelScore;
    }

    /**
     * @param IndelScore the IndelScore to set
     */
    public void setIndelScore(int IndelScore) {
        this.IndelScore = IndelScore;
    }

    /**
     * @return the str1
     */
    public String getStr1() {
        return str1;
    }

    /**
     * @param str1 the str1 to set
     */
    public void setStr1(String str1) {
        this.str1 = str1;
    }

    /**
     * @return the str2
     */
    public String getStr2() {
        return str2;
    }

    /**
     * @param str2 the str2 to set
     */
    public void setStr2(String str2) {
        this.str2 = str2;
    }

    /**
     * @return the alignedStr1
     */
    public String getAlignedStr1() {
        return alignedStr1.toString();
    }

    /**
     * @param alignedStr1 the alignedStr1 to set
     */
    public void setAlignedStr1(String alignedStr1) {
        this.alignedStr1 = new StringBuilder(alignedStr1);
    }

    /**
     * @return the alignedStr2
     */
    public String getAlignedStr2() {
        return alignedStr2.toString();
    }

    /**
     * @param alignedStr2 the alignedStr2 to set
     */
    public void setAlignedStr2(String alignedStr2) {
        this.alignedStr2 = new StringBuilder(alignedStr2);
    }

    /**
     * @return the alignedStatus
     */
    public String getAlignedStatus() {
        return alignedStatus.toString();
    }

    /**
     * @param alignedStatus the alignedStatus to set
     */
    public void setAlignedStatus(String alignedStatus) {
        this.alignedStatus = new StringBuilder(alignedStatus);
    }

    /**
     * @return the debugLevel
     */
    public int getDebugLevel() {
        return debugLevel;
    }

    /**
     * @param debugLevel the debugLevel to set
     */
    public void setDebugLevel(int debugLevel) {
        this.debugLevel = debugLevel;
    }

    /**
     * @return the minI
     */
    public int getMinI() {
        return minI;
    }

    /**
     * @param minI the minI to set
     */
    public void setMinI(int minI) {
        this.minI = minI;
    }

    /**
     * @return the minJ
     */
    public int getMinJ() {
        return minJ;
    }

    /**
     * @param minJ the minJ to set
     */
    public void setMinJ(int minJ) {
        this.minJ = minJ;
    }

    /**
     * @return the nMatches
     */
    public int getnMatches() {
        return nMatches;
    }

    /**
     * @param nMatches the nMatches to set
     */
    public void setnMatches(int nMatches) {
        this.nMatches = nMatches;
    }

    /**
     * @return the nMismatches
     */
    public int getnMismatches() {
        return nMismatches;
    }

    /**
     * @param nMismatches the nMismatches to set
     */
    public void setnMismatches(int nMismatches) {
        this.nMismatches = nMismatches;
    }

    /**
     * @return the nIndels
     */
    public int getnIndels() {
        return nIndels;
    }

    /**
     * @param nIndels the nIndels to set
     */
    public void setnIndels(int nIndels) {
        this.nIndels = nIndels;
    }
}
