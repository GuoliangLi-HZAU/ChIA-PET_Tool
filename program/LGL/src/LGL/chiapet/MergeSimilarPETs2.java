/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package LGL.chiapet;

/**
 *
 * @author gli
 */

import LGL.data.HIT2;
import LGL.data.PET;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.util.Vector;

/*
 *
 * 1) the input file needs to be sorted by
 *    a) chromosome in the read1,
 *    b) chromosome in the read2,
 *    c) strand in the read1,
 *    d) strand in the read2
 *
 */
public class MergeSimilarPETs2 {

    int debugLevel = 4;
    Calendar rightNow = Calendar.getInstance();

    int sortingLabel = 0; // no sorting on the head and tail inside a single PET [default]
    // The different PETs are sorted by chromosomes and strands.
    Vector<PET> pets = new Vector<PET>();
    int distance = 2;

    public MergeSimilarPETs2(String inputPetFile, String outputFile, int distance, int sortingLabel) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start MergeSimilarPETs2 ... ");

        this.distance     = distance;
        this.sortingLabel = sortingLabel;

        // load inter-ligation PETs
        String chrom1 = "";
        String chrom2 = "";
        char strand1 = '!';
        char strand2 = '!';
        int petIndex = 0; // artificial index
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputPetFile))));
        PET.save(outputFile, pets, false);
        String line;
        int nPETs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 0) // skip the short lines
            {
                continue;
            }
            String[] fields = line.split("\t");

            HIT2 hit1 = new HIT2(fields[0], Integer.parseInt(fields[1]), fields[2].charAt(0));
            HIT2 hit2 = new HIT2(fields[3], Integer.parseInt(fields[4]), fields[5].charAt(0));
            double weight = 1.0;
            if (fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            PET newPet = new PET(hit1, hit2, weight, petIndex, this.sortingLabel);
            if (fields.length >= 8) {
                newPet.setIndex(Integer.parseInt(fields[7]));
            }
            if ((!chrom1.equalsIgnoreCase(newPet.getHead().getChrom())) 
                    || (!chrom2.equalsIgnoreCase(newPet.getTail().getChrom()))
                    || (strand1 != newPet.getHead().getStrand())
                    || (strand2 != newPet.getTail().getStrand())) {
                // merge the same and similar PETs
                pets = PET.mergeSimilarPets(pets, distance);
                PET.save(outputFile, pets, true);
                pets.clear();
                chrom1 = newPet.getHead().getChrom();
                chrom2 = newPet.getTail().getChrom();
                strand1 = newPet.getHead().getStrand();
                strand2 = newPet.getTail().getStrand();
            }

            pets.add(newPet);
            nPETs++;
            if (nPETs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000) + "M PETs read from " + inputPetFile);
            }
        }
        fileIn.close();
        pets = PET.mergeSimilarPets(pets, distance);
        PET.save(outputFile, pets, true);
        System.out.println("number of pets = " + nPETs);
    }


    public static void main(String[] args) throws IOException {
        if (args.length == 3) {
            int distance = Integer.parseInt(args[2]);
            int sortingLabel = 0; // no sorting [default]

            new MergeSimilarPETs2(args[0], args[1], distance, sortingLabel);
        }
        else if(args.length==4) {
            int distance = Integer.parseInt(args[2]);
            int sortingLabel = Integer.parseInt(args[3]);
            new MergeSimilarPETs2(args[0], args[1],distance, sortingLabel);
        }
        else {
            System.out.println("Usage: java MergeSimilarPETs2 <PET file> <merged PET file> <distance> [<sortingLabel>]");
            System.out.println("Input PET file is sorted by the chromosome pairs");
            System.out.println("<distance> is the required maximum distance between different PETs");
            System.out.println("[<sortingLabel>] is optional. The default is no sorting for the head and tail");
            System.out.println("sortingLabel:  1: ascending order");
            System.out.println("              -1: descending order");
            System.out.println("              0 or other: no sorting  [default]");
            System.exit(1);
        }
    }
}
