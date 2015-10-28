/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.chiapet;

import LGL.data.HIT;
import LGL.data.HIT2;
import LGL.data.PET;
import LGL.util.SeqUtil;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Calendar;

/**
 *
 * @author ligl
 */
public class PetClassification {

    int debugLevel = 3;
    Calendar rightNow = Calendar.getInstance();
    int selfLigationCutoff = 3000;

    // iPET: for inter-ligation PETs
    // sPET: for self-ligation PETs
    // oPET: for other PETs in short distance
    public PetClassification(String inFile, String iPEToutFile, String sPEToutFile, String oPEToutFile, int cutoff) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start PetClassification ... ");

        this.selfLigationCutoff = cutoff;

        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] " + "start loading PETs from " + inFile);

        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inFile))));
        PrintWriter iPETfileOut = new PrintWriter(new BufferedWriter(new FileWriter(iPEToutFile)));
        PrintWriter sPETfileOut = new PrintWriter(new BufferedWriter(new FileWriter(sPEToutFile)));
        PrintWriter oPETfileOut = new PrintWriter(new BufferedWriter(new FileWriter(oPEToutFile)));

        String line;
        int nPETs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            String[] fields = line.split("\t");
            HIT2 hit1 = new HIT2(fields[0], Integer.parseInt(fields[1]), fields[2].charAt(0));
            HIT2 hit2 = new HIT2(fields[3], Integer.parseInt(fields[4]), fields[5].charAt(0));
            double weight = 1.0;
            if(fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            PET pet = new PET(hit1, hit2, weight);
            classify(pet, iPETfileOut, sPETfileOut, oPETfileOut);
            nPETs++;
            if (nPETs % 1000000 == 0) {
                iPETfileOut.flush();
                sPETfileOut.flush();
                oPETfileOut.flush();

                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000) + "M PETs read from " + inFile);
            }
        }
        fileIn.close();
        iPETfileOut.close();
        sPETfileOut.close();
        oPETfileOut.close();
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000.0) + "M PETs read from " + inFile);
    }

    // this applies to the ChIA-PET data with read1 fliped （reverse complements)
    void classify(PET pet, PrintWriter iPETfileOut, PrintWriter sPETfileOut, PrintWriter oPETfileOut, int flip_head) {
        if (pet.getHead().getChrom().compareToIgnoreCase(pet.getTail().getChrom()) != 0) {
            // different chromosomes, as inter-chromosomal inter-ligations
            iPETfileOut.println(pet.toString());
        } else if (HIT.calculateDistance((HIT) pet.getHead(), (HIT) pet.getTail()) > selfLigationCutoff) {
            // same chromosome, long distance, as intra-chromosomal inter-ligations
            iPETfileOut.println(pet.toString());
        } else if (pet.getHead().getStrand() != pet.getTail().getStrand()) {
            // same chromosome, short distance, different strands, as other PETs in short distance
            oPETfileOut.println(pet.toString());
        } else if (SeqUtil.isForwardStrand(pet.getHead().getStrand()) == true) {
            if (pet.getHead().getLoci() < pet.getTail().getLoci()) {
                // same chromosome, same positive strands, short distance, head is before tail, as self-ligations
                sPETfileOut.println(pet.toString());
            } else {
                // same chromosome, same positive strands, short distance, head is after tail, as other inter-ligations
                oPETfileOut.println(pet.toString());
            }
        } else {
            if (pet.getHead().getLoci() > pet.getTail().getLoci()) {
                // same chromosome, same negative strands, short distance, head is before tail, as self-ligations
                sPETfileOut.println(pet.toString());
            } else {
                // same chromosome, same negative strands, short distance, head is after tail, as other inter-ligations
                oPETfileOut.println(pet.toString());
            }
        }
    }

    // self-ligation PETs are in the mode: <-- -->
    void classify(PET pet, PrintWriter iPETfileOut, PrintWriter sPETfileOut, PrintWriter oPETfileOut) {
        if (pet.getHead().getChrom().compareToIgnoreCase(pet.getTail().getChrom()) != 0) {
            // different chromosomes, as inter-chromosomal inter-ligations
            iPETfileOut.println(pet.toString());
        } else if (HIT.calculateDistance((HIT) pet.getHead(), (HIT) pet.getTail()) > selfLigationCutoff) {
            // same chromosome, long distance, as intra-chromosomal inter-ligations
            iPETfileOut.println(pet.toString());
        } else if (pet.getHead().getStrand() == pet.getTail().getStrand()) {
            // same chromosome, short distance, same strands, as other PETs in short distance
            oPETfileOut.println(pet.toString());
        } else if (SeqUtil.isForwardStrand(pet.getHead().getStrand()) == true) {
            if (pet.getHead().getLoci() < pet.getTail().getLoci()) {
                // same chromosome, short distance, read1:plus:small, read2:minus:large, as other PETs in short distance
                oPETfileOut.println(pet.toString());
            } else {
                // same chromosome, short distance, read1:plus:large, read2:minus:small, as self-ligations
                sPETfileOut.println(pet.toString());
            }
        } else {
            if (pet.getHead().getLoci() > pet.getTail().getLoci()) {
                // same chromosome, short distance, read1:minus:large, read2:plus:small, as other PETs in short distance
                oPETfileOut.println(pet.toString());
            } else {
                // same chromosome, short distance, read1:minus:small, read2:plus:large, as self-ligations
                sPETfileOut.println(pet.toString());
            }
        }
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 5) {
            new PetClassification(args[0], args[1], args[2], args[3], Integer.parseInt(args[4]));
        } else {
            System.out.println("Usage: java PetClassification <pet_file_before_parse> <iPet_file_after_parse> <sPet_file_after_parse> <oPet_file_after_parse> <selfLigationCutoff>");
            System.out.println("       <pet_file_before_parse>: input file in pet format");
            System.out.println("       <iPet_file_after_parse>: output file with inter-ligation PETs");
            System.out.println("       <sPet_file_after_parse>: output file with self-ligation PETs");
            System.out.println("       <oPet_file_after_parse>: output file with other PETs in short distance");
            System.exit(1);
        }
    }
}
