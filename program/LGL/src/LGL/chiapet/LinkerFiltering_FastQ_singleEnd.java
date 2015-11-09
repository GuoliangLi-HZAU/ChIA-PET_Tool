/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.chiapet;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 *
 * Input:
 * 1) a pair of FastQ files
 * 2) a file with linker sequences with (nLinker) linkers
 * 3) configurations
 *     a) minimum_linker_alignment_score: default >= 14
 *     b) minimum_tag_length: 20 bp
 *     c) flip_head_tag: default: 1: Yes (default), 0: No
 *     d) flip_tail_tag: default: 1: Yes, 0: No (default)
 *     e) debug_level: integer (default: 1)
 *
 * Description:
 * FastQ description: Refer to Wikipedia http://en.wikipedia.org/wiki/FASTQ_format
 * There are two formats for linker file
 * 1) just with one column: the column contains the half-linkers
 * 2) there are 3 columns:
 *    1st: the first column is the half-linkers,
 *    2nd: the second column is the bar-code start position (the first index is 1)
 *    3rd: the third column is the length of the bar code
 *
 * Assumptions:
 * 1) there are nLinker possible linker sequences
 * 2) all the linkers have the same length
 * 3) all the linkers are 5bps or more
 *
 * Output:
 * all the output files are in one folder
 * 1) nLinker*nLinker pairs of FastQ files
 * 2) debug file
 *
 * command line:
 * java LinkerFiltering_FastQ_PET <Input information file>
 *
 */
import LGL.align.LocalAlignment;
import LGL.data.FastQ;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;

/**
 *
 * @author ligl
 */
public class LinkerFiltering_FastQ_singleEnd {

    int debug_level = 1;
    String fastQFile = null;
    String linkerFile = null;
    String inputInfoFile = null;
    String outputFolder = null;
    String outputPrefix = "output";
    String[] linkers = null;
    int nLinkers = 0;
    int[] barCodeStart = null;
    int[] barCodeLength = null; // default: 4
    int maxLinkerLength;
    int minLinkerLength;
    int minimum_linker_alignment_score = 14; // minimum_linker_alignment score: default 14
    int minimum_tag_length = 20; // default: 20
    int minSecondBestScoreDiff = 3;
    int output_data_with_ambiguous_linker_info = 1; // default: output the data with ambiguous linker info
    int flip_head_tag = 0; // default: flip head tag
    int flip_tail_tag = 0; // default: Not flip tail tag
    int linker_length = 16;
    int read1_end = 20;

    int[][] scoreDistribution;
    int[][] secondBestScoreDiffDistribution;
    int maxSecondBestScoreDiff = 0;

    int maximum_length_tag1 = 1000; // default: 1000bp
    int maximum_length_tag2 = 21;
    int maxRealTagLength = 0; // used to control the output length for the tag length distribution
    int[][] tagLengthDistribution;

    int[][] linkerCompositionDistribution;
    int nAmbiguousLinkerComposition = 0;

    Calendar rightNow = Calendar.getInstance();
    LocalAlignment localAligner; // for linker alignment
    PrintWriter debug_output = null;

    String[] letter = {"A", "B", "C", "D", "E", "F"};

    public LinkerFiltering_FastQ_singleEnd(String inputInfoFile, String outputFolder, String outputPrefix) throws IOException {
        this.inputInfoFile = inputInfoFile;
        this.outputFolder = outputFolder;
        this.outputPrefix = outputPrefix;
        System.out.println("[" + rightNow.getTime().toString() + "] start LinkerFiltering_FastQ_singleEnd ... ");
        // read input information
        readInputInformation();
        // test input information
        testInputInformation();
        //if(false)
        {
            // read linkers
            readLinkers();
            // output running information
            outputRunningInfo();
            localAligner = new LocalAlignment(maxLinkerLength, maxLinkerLength);

            // initiate the distribution arrays
            // tag length distribution
            tagLengthDistribution = new int[2][maximum_length_tag1];
            Arrays.fill(tagLengthDistribution[0], 0);
            Arrays.fill(tagLengthDistribution[1], 0);

            scoreDistribution = new int[2][this.maxLinkerLength + 1];
            Arrays.fill(scoreDistribution[0], 0);
            Arrays.fill(scoreDistribution[1], 0);
            secondBestScoreDiffDistribution = new int[2][this.maxLinkerLength * 2 + 1];
            Arrays.fill(secondBestScoreDiffDistribution[0], 0);
            Arrays.fill(secondBestScoreDiffDistribution[1], 0);
            linker_length = linkers[0].length();
            read1_end = linker_length + maximum_length_tag1;
            filterSequenceByLinker();
            printDistribution();
        }
    }

    private void printDistribution() throws IOException {
        PrintWriter fileOut = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".linker_alignment_score_distribution.txt")));
        //fileOut.println("Distribution of linker alignment scores");
        for (int i = 0; i < scoreDistribution[0].length; i++) {
            fileOut.println(i + "\t" + scoreDistribution[0][i] + "\t" + scoreDistribution[1][i]);
        }
        fileOut.close();

        fileOut = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".linker_alignment_score_difference_distribution.txt")));
        //fileOut.println("\nDistribution of alignment score differences between best alignment and second best alignment");
        for (int i = 0; i < secondBestScoreDiffDistribution[0].length && i <= maxSecondBestScoreDiff; i++) {
            fileOut.println(i + "\t" + secondBestScoreDiffDistribution[0][i] + "\t" + secondBestScoreDiffDistribution[1][i]);
        }
        fileOut.close();

        fileOut = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".tag_length_distribution.txt")));
        //fileOut.println("\nDistribution of tag lengths");
        for (int i = 0; i < tagLengthDistribution[0].length && i <= maxRealTagLength; i++) {
            fileOut.println(i + "\t" + tagLengthDistribution[0][i] + "\t" + tagLengthDistribution[1][i]);
        }
        fileOut.close();

        fileOut = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".linker_composition_distribution.txt")));
        //fileOut.println("\nDistribution of linker compositions");
        String lin1 = null;
        String lin2 = null;
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                if (i == 0) {
                    lin1 = "A";
                } else if (i == 1) {
                    lin1 = "B";
                }
                if (j == 0) {
                    lin2 = "A";
                } else if (j == 1) {
                    lin2 = "B";
                }

                //int newi = i+1;
                //int newj = j+1;
                fileOut.print(lin1 + "_" + lin2 + "\t");
            }
        }
        fileOut.println("Ambiguous\tTotal");
        int total = 0;
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut.print(linkerCompositionDistribution[i][j] + "\t");
                total += linkerCompositionDistribution[i][j];
            }
        }
        total += nAmbiguousLinkerComposition;
        fileOut.println(nAmbiguousLinkerComposition + "\t" + total);
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut.print(String.format("%.2f", (100.0 * linkerCompositionDistribution[i][j] / total)) + "%" + "\t");
            }
        }
        fileOut.println(String.format("%.2f", (100.0 * nAmbiguousLinkerComposition / total)) + "%" + "\t" + "100%");

        fileOut.close();

    }

    public void filterSequenceByLinker() throws IOException {
        BufferedReader fastqFileIn = new BufferedReader(new InputStreamReader(new FileInputStream(this.fastQFile)));
        // for every linker pair, output a pair of FastQ files
        PrintWriter[] fileOut = new PrintWriter[nLinkers * nLinkers * 2 + 2];
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut[(i * nLinkers + j) * 2] = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + "." + (i + 1) + "_" + (j + 1) + ".R1.fastq")));
                fileOut[(i * nLinkers + j) * 2 + 1] = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + "." + (i + 1) + "_" + (j + 1) + ".R2.fastq")));
            }
        }
        if (this.output_data_with_ambiguous_linker_info == 1) {
            fileOut[nLinkers * nLinkers * 2] = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".ambiguous.R1.fastq")));
            fileOut[nLinkers * nLinkers * 2 + 1] = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".ambiguous.R2.fastq")));
        }
        if (debug_level >= 2) {
            debug_output = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + "debug.LinkerFiltering_FastQ_PET.txt")));
        }

        long nPETs = 0;
        FastQ fastQ1 = FastQ.load(fastqFileIn);

        while (fastQ1.getFastq() != null) {
            processOnePET(fastQ1, fileOut);
            nPETs++;
            if (nPETs % 1000000 == 0) {
                for (int i = 0; i < nLinkers; i++) {
                    for (int j = 0; j < nLinkers; j++) {
                        fileOut[(i * nLinkers + j) * 2].flush();
                        fileOut[(i * nLinkers + j) * 2 + 1].flush();
                    }
                }
                if (this.output_data_with_ambiguous_linker_info == 1) {
                    fileOut[nLinkers * nLinkers * 2].flush();
                    fileOut[nLinkers * nLinkers * 2 + 1].flush();
                }
                System.out.println((nPETs / 1000000) + " Million PETs processed");
                if (debug_level >= 2) {
                    debug_output.flush();
                }
            }

            fastQ1 = FastQ.load(fastqFileIn);
        }
        fastqFileIn.close();

        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut[(i * nLinkers + j) * 2].close();
                fileOut[(i * nLinkers + j) * 2 + 1].close();
            }
        }
        if (this.output_data_with_ambiguous_linker_info == 1) {
            fileOut[nLinkers * nLinkers * 2].close();
            fileOut[nLinkers * nLinkers * 2 + 1].close();
        }
        PrintWriter basicOut = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".basic_statistics.txt")));
        basicOut.write("Total PETs\t" + nPETs + "\n");
        basicOut.close();
        if (debug_level >= 2) {
            debug_output.close();
        }
    }

    public void processOnePET(FastQ fastQ, PrintWriter[] fileOut) throws IOException {
        int[] results_1;
        if (read1_end < fastQ.getFastq()[1].length()) {
            results_1 = processOneSequence(fastQ.getFastq()[1].substring(0, read1_end), 0);
        } else {
            results_1 = processOneSequence(fastQ.getFastq()[1], 0);
        }
        String revComRead = LGL.util.SeqUtil.revComplement(fastQ.getFastq()[1]);
        StringBuilder revComQuality = new StringBuilder(fastQ.getFastq()[3]).reverse();
        int[] results_2 = processOneSequence(revComRead, 1);
        /* result format illustration
         results[0] = tag_Start;
         results[1] = tag_End;
         results[2] = bestLinkerIndex;
         results[3] = bestScore; // best linker alignment score
         results[4] = secondBestScoreDiff; // difference between the best alignment score and the second-best alignment score
         results[5] = minJ; // start index in sequence for the best local alignment
         results[6] = minI; // start index in linker for the best local alignment
         results[7] = barcodeStatus; //1: the barcode in the reads is the same as in the designed sequence
         *
         */

        PrintWriter output1;
        PrintWriter output2;
        int tag_length_1 = results_1[1] - results_1[0];
        int tag_length_2 = results_2[1] - results_2[0];
        // criteria to accept the linker filtering results
        // 1. linker sequence has alignment score larger than or equal to the minimum linker alignment score
        // 2. tag length is larger than or equal to the minimum tag length
        // 3. tag length is smaller than or equal to the max tag length
        // 4. the difference of alignment score between the best alignment and the second-best alignment to the linker sequences is larger than the required cutoff
        // 5. the barcode in the reads must be the same as the designed linker barcode
        if ((results_1[3] >= minimum_linker_alignment_score)
                && ((results_2[3]) >= minimum_linker_alignment_score)
                && (tag_length_1 >= minimum_tag_length)
                && (tag_length_2 >= minimum_tag_length)
                && (tag_length_1 <= maximum_length_tag1)
                && (tag_length_2 <= maximum_length_tag2)
                && (results_1[4] >= minSecondBestScoreDiff)
                && (results_2[4] >= minSecondBestScoreDiff)
                && (results_1[7] == 1)
                && (results_2[7] == 1)) {
            output1 = fileOut[(results_1[2] * nLinkers + results_2[2]) * 2];
            output2 = fileOut[(results_1[2] * nLinkers + results_2[2]) * 2 + 1];
            // update linker composition distribution
            linkerCompositionDistribution[results_1[2]][results_2[2]]++;

            // output head
            if (flip_head_tag == 1) {
                output1.println(fastQ.getFastq()[0]);
                output1.println(LGL.util.SeqUtil.revComplement(fastQ.getFastq()[1].substring(results_1[0], results_1[1])));
                output1.println(fastQ.getFastq()[2]);
                output1.println(new StringBuilder(fastQ.getFastq()[3].substring(results_1[0], results_1[1])).reverse());
            } else {
                output1.println(fastQ.getFastq()[0]);
                output1.println(fastQ.getFastq()[1].substring(results_1[0], results_1[1]));
                output1.println(fastQ.getFastq()[2]);
                output1.println(fastQ.getFastq()[3].substring(results_1[0], results_1[1]));
            }

            // output tail
            if (flip_tail_tag == 0) {
                output2.println(fastQ.getFastq()[0]);
                output2.println(revComRead.substring(results_2[0], results_2[1]));
                output2.println(fastQ.getFastq()[2]);
                output2.println(revComQuality.substring(results_2[0], results_2[1]));
            } else {
                output2.println(fastQ.getFastq()[0]);
                output2.println(LGL.util.SeqUtil.revComplement(revComRead.substring(results_2[0], results_2[1])));
                output2.println(fastQ.getFastq()[2]);
                output2.println(new StringBuilder(revComQuality.substring(results_2[0], results_2[1])).reverse());
            }
        } else { // PETs with ambiguous linkers
            nAmbiguousLinkerComposition++;
            if (this.output_data_with_ambiguous_linker_info == 1) {
                output1 = fileOut[nLinkers * nLinkers * 2];
                output2 = fileOut[nLinkers * nLinkers * 2 + 1];

                output1.println(fastQ.getFastq()[0]);
                output1.println(fastQ.getFastq()[1]);
                output1.println(fastQ.getFastq()[2]);
                output1.println(fastQ.getFastq()[3]);

                output2.println(fastQ.getFastq()[0]);
                output2.println(LGL.util.SeqUtil.revComplement(fastQ.getFastq()[1]));
                output2.println(fastQ.getFastq()[2]);
                output2.println(LGL.util.SeqUtil.revComplement(fastQ.getFastq()[3]));
            }
        }

        if (debug_level >= 2) {
            // output the sequences and linker alignment for debugging purpose
            debug_output.print(fastQ.getFastq()[1].substring(results_1[0], results_1[1]) + "\t" + fastQ.getFastq()[1].substring(results_1[1]) + "\t" + results_1[2] + "\t" + results_1[3] + "\t" + results_1[4] + "\t" + results_1[5] + "\t" + results_1[6] + "\t");
            debug_output.print(fastQ.getFastq()[0]);
        }
    }

    int nOutput = 10;
    int iOutput = 0;

    // 1: barcode test passed; 0: barcode test failed
    private int testBarcode(LocalAlignment a, int offset, int len) {
        if ((offset < 0) || (len < 0)) {
            // either offset or len is less than 0, barcode test is not applicable
            return 1;
        }
        // ****** assumptions ******
        // 1) the aligned strings are from the beginning of the original sequences, and to the end of the best local alignment
        // 2) the aligned str1 is the linker sequence
        // 3) the insertions and deletions are represented with '-'
        // 4) matched alignment status is represented with '|'
        String alignedStr1 = a.getAlignedStr1();
        String alignedStatus = a.getAlignedStatus();
        StringBuilder trimmedAlignedStatus = new StringBuilder();
        for (int i = 0; i < alignedStr1.length(); i++) {
            if (alignedStr1.charAt(i) != '-') {
                trimmedAlignedStatus.append(alignedStatus.charAt(i));
            }
        }
        if (debug_level >= 2) {
            if (iOutput < nOutput) {
                System.out.println("originalStr1: " + a.getStr1());
                System.out.println("alignedStr1:  " + alignedStr1);
                System.out.println("status:       " + alignedStatus);
                System.out.println("alignedStr2:  " + a.getAlignedStr2());
                System.out.println("originalStr2: " + a.getStr2());
                System.out.println("trimmedAlignedStatus: " + trimmedAlignedStatus);
                iOutput++;
            }
        }

        int barcodePassed = 1; // 1: barcode test passed; 0: barcode test failed
        if (trimmedAlignedStatus.length() < offset + len - 1) {
            // the barcode is not completely covered
            barcodePassed = 0;
        } else {
            int i = offset - 1;
            for (int j = 0; j < len; j++) {
                if (trimmedAlignedStatus.charAt(i) != '|') {
                    barcodePassed = 0;
                    break;
                }
                i++;
            }
        }

        return barcodePassed;
    }

    // iRead: 0 for read1, 1 for read2
    public int[] processOneSequence(String seq, int iRead) throws IOException {
        // align different linkers
        int bestScore = -1;
        int secondBestScore = -1;
        int bestLinkerIndex = -1;
        int minI = -1; // index in linker
        int minJ = -1; // index in sequence
        int barcodeStatus = -1;

        for (int i = 0; i < this.linkers.length; i++) {
            localAligner.align(this.linkers[i], seq);
            int score = localAligner.getMaxScore();

            if (bestScore < score) {
                secondBestScore = bestScore;
                bestScore = score;
                bestLinkerIndex = i;
                minI = localAligner.getMinI(); // index in linker
                minJ = localAligner.getMinJ(); // index in sequence
                barcodeStatus = testBarcode(localAligner, this.barCodeStart[i], this.barCodeLength[i]);
            } else if (secondBestScore < score) {
                secondBestScore = score;
            }
        }
        if ((bestScore >= 0) && (bestScore <= this.maxLinkerLength)) {
            scoreDistribution[iRead][bestScore]++;
        }
        int secondBestScoreDiff = bestScore - secondBestScore;

        if ((secondBestScoreDiff >= 0) && (secondBestScoreDiff <= 2 * this.maxLinkerLength)) {
            secondBestScoreDiffDistribution[iRead][secondBestScoreDiff]++;
        }
        if (maxSecondBestScoreDiff < secondBestScoreDiff) {
            maxSecondBestScoreDiff = secondBestScoreDiff;
        }

        int tag_Start = 0;

        int tag_End = minJ - minI;
        if (tag_End < 0) {
            tag_End = 0;
        }
        if (iRead == 0) {
            if (tag_End < 0) {
                tagLengthDistribution[iRead][0]++;
            } else if (tag_End >= maximum_length_tag1) {
                tagLengthDistribution[iRead][maximum_length_tag1 - 1]++;
            } else {
                tagLengthDistribution[iRead][tag_End]++;
            }
        } else { // iRead == 1
            if (tag_End < 0) {
                tagLengthDistribution[iRead][0]++;
            } else if (tag_End >= maximum_length_tag2) {
                tagLengthDistribution[iRead][maximum_length_tag2 - 1]++;
                tag_Start = tag_End - maximum_length_tag2;
            } else {
                tagLengthDistribution[iRead][tag_End]++;
            }

        }
// record the real tag length, for output the tag length distribution purpose
        if (maxRealTagLength < tag_End - tag_Start) {
            maxRealTagLength = tag_End - tag_Start;
        }

        // sequence index starts with 0
        int[] results = new int[8];
        results[0] = tag_Start;
        results[1] = tag_End;
        results[2] = bestLinkerIndex;
        results[3] = bestScore;
        results[4] = secondBestScoreDiff;
        results[5] = minJ;
        results[6] = minI;
        results[7] = barcodeStatus;

        return results;
    }

    public void readLinkers() throws IOException {
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(this.linkerFile)));
        ArrayList<String> tempLinkers = new ArrayList<String>();
        String line;

        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            tempLinkers.add(line);
        }
        fileIn.close();

        this.linkers = new String[tempLinkers.size()];
        this.maxLinkerLength = 0;
        this.minLinkerLength = Integer.MAX_VALUE;

        this.barCodeStart = new int[tempLinkers.size()];
        this.barCodeLength = new int[tempLinkers.size()];

        for (int i = 0; i < linkers.length; i++) {
            String[] fields = tempLinkers.get(i).split("\t");
            this.linkers[i] = fields[0];

            // shortest linker length
            if (this.minLinkerLength > this.linkers[i].length()) {
                this.minLinkerLength = this.linkers[i].length();
            }
            // longest linker length
            if (this.maxLinkerLength < this.linkers[i].length()) {
                this.maxLinkerLength = this.linkers[i].length();
            }

            //bar code start
            if (fields.length >= 2) {
                this.barCodeStart[i] = Integer.parseInt(fields[1]);
            } else {
                this.barCodeStart[i] = -1;
            }
            //bar code length
            if (fields.length >= 3) {
                this.barCodeLength[i] = Integer.parseInt(fields[2]);
            } else {
                this.barCodeLength[i] = -1;
            }
        }
        nLinkers = linkers.length;
        if (nLinkers <= 0) {
            System.out.println("No linker sequence information.\tStop!!!");
            System.exit(0);
        }
        if (nLinkers > 100) {
            System.out.println("Too many linkers. Please check...\tStop!!!");
            System.exit(0);
        }
        if (this.minLinkerLength < 5) {
            System.out.println("the shortest linker is less than 5bp.\tStop!!!");
            System.exit(0);
        }
        // create linker composition distribution matrix, and initiate the values to 0
        linkerCompositionDistribution = new int[nLinkers][nLinkers];
        for (int i = 0; i < nLinkers; i++) {
            Arrays.fill(linkerCompositionDistribution[i], 0);
        }
    }

    public void outputRunningInfo() throws IOException {
        PrintWriter fileOut = new PrintWriter(new FileOutputStream(new File(this.outputFolder, this.outputPrefix + ".runningInformation.LinkerFiltering_FastQ_PET.txt")));
        fileOut.println("Fastq file\t" + new File(this.fastQFile).getAbsolutePath());
        fileOut.println("Linker file\t" + new File(this.linkerFile).getAbsolutePath());
        for (int i = 0; i < linkers.length; i++) {
            fileOut.println("Linker" + letter[i] + "\t" + this.linkers[i] + " " + this.barCodeStart[i] + " " + this.barCodeLength[i]);
        }
        fileOut.println("Output folder\t" + new File(this.outputFolder).getAbsolutePath());
        fileOut.println("Output prefix\t" + this.outputPrefix);
        fileOut.println("Minimum linker alignment score\t" + this.minimum_linker_alignment_score);
        fileOut.println("Minimum tag length\t" + this.minimum_tag_length);
        fileOut.println("Maximum tag1 length\t" + this.maximum_length_tag1);
        fileOut.println("Maximum tag2 length\t" + this.maximum_length_tag2);
        fileOut.println("Minimum SecondBestScore difference\t" + this.minSecondBestScoreDiff);
        fileOut.println("Output data with ambiguous linker info\t" + this.output_data_with_ambiguous_linker_info);
        // fileOut.println("Flip head tag\t" + this.flip_head_tag);
        // fileOut.println("Flip tail tag\t" + this.flip_tail_tag);
        //fileOut.println("Debug level\t"   + this.debug_level);
        fileOut.close();
    }

    public void testInputInformation() throws IOException {
        boolean requiredInputMissed = false;
        // test required input parameters
        if (this.fastQFile == null) {
            System.out.println("No FASTQ file input!");
            requiredInputMissed = true;
        }
        if (this.linkerFile == null) {
            System.out.println("No linker file input!");
            requiredInputMissed = true;
        }
        // if required inputs are missed, exit
        if (requiredInputMissed == true) {
            System.exit(0);
        }

        // compare the maximum_tag_length and minimum_tag_length
        if (maximum_length_tag1 < minimum_tag_length || maximum_length_tag2 < minimum_tag_length) {
            System.out.println("maximum_tag_length is smaller than the minimum_tag_length! Stop ...");
            System.exit(0);
        }

        // test optional input parameters
        if (this.outputFolder == null) {
            System.out.println("No output folder! Set it as \"out\"");
            this.outputFolder = "out";
        }
        // test whether outputFolder exists
        File fOutputFolder = new File(this.outputFolder);
        if (!fOutputFolder.exists()) { // no this folder
            fOutputFolder.mkdir(); // create this folder
        } else if (!fOutputFolder.isDirectory()) { // a file exists with the same name, but not a folder
            System.out.println(this.outputFolder + " exists, but it is not a directory!!! please check...");
            System.exit(0);
        }
        this.outputFolder = fOutputFolder.getPath();
    }

    public void readInputInformation() throws IOException {
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(this.inputInfoFile)));
        String line;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String[] fields = line.split("[ \t][ \t]*");
            //fastq_file
            if (fields[0].equalsIgnoreCase("fastq_file")) {
                this.fastQFile = fields[1];
            }
            //linker_file
            if (fields[0].equalsIgnoreCase("linker_file")) {
                this.linkerFile = fields[1];
            }
            //minimum_linker_alignment_score
            if (fields[0].equalsIgnoreCase("minimum_linker_alignment_score")) {
                this.minimum_linker_alignment_score = Integer.parseInt(fields[1]);
            } //minimum_tag_length
            if (fields[0].equalsIgnoreCase("minimum_tag_length")) {
                this.minimum_tag_length = Integer.parseInt(fields[1]);
            } // maximum_length_tag1
            if (fields[0].equalsIgnoreCase("maximum_length_tag1")) {
                this.maximum_length_tag1 = Integer.parseInt(fields[1]);
            } // maximum_length_tag2
            if (fields[0].equalsIgnoreCase("maximum_length_tag2")) {
                this.maximum_length_tag2 = Integer.parseInt(fields[1]);
            } //minSecondBestScoreDiff
            if (fields[0].equalsIgnoreCase("minSecondBestScoreDiff")) {
                this.minSecondBestScoreDiff = Integer.parseInt(fields[1]);
            } //output_data_with_ambiguous_linker_info
            if (fields[0].equalsIgnoreCase("output_data_with_ambiguous_linker_info")) {
                this.output_data_with_ambiguous_linker_info = Integer.parseInt(fields[1]);
            } //flip_head_tag
            if (fields[0].equalsIgnoreCase("flip_head_tag")) {
                this.flip_head_tag = Integer.parseInt(fields[1]);
            } //flip_tail_tag
            if (fields[0].equalsIgnoreCase("flip_tail_tag")) {
                this.flip_tail_tag = Integer.parseInt(fields[1]);
            } //debug_level
            if (fields[0].equalsIgnoreCase("debug_level")) {
                this.debug_level = Integer.parseInt(fields[1]);
            }
        }
        fileIn.close();
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 3) {
            new LinkerFiltering_FastQ_singleEnd(args[0], args[1], args[2]);
        } else {
            System.out.println("Usage: java LinkerFiltering_FastQ_singleEnd <Input Information file> <outputFolder> <outputPrefix>");
            System.out.println("Configurations");
            System.out.println(" * a) minimum_linker_alignment_score: default 14");
            System.out.println(" * b) minimum_tag_length: 20 bp");
            System.out.println(" * c) flip_head_tag: default: 1: Yes (default), 0: No");
            System.out.println(" * d) flip_tail_tag: default: 1: Yes, 0: No (default)");
            System.exit(0);
        }
    }
}
