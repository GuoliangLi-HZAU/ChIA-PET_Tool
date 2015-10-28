/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.file;

import LGL.data.CLUSTER;
import LGL.data.HIT2;
import LGL.data.PET;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Calendar;
import java.util.Hashtable;

public class Pet2Cluster1 {

    public static void main(String[] args) throws IOException {

        if (args.length < 4) {
            System.out.println("Usage: java Pet2Cluster1 <input_pet_file> <output_cluster_file> <extensionLength> <chromLengthInfo> [<sortingLabel>]");
            System.out.println("[sortingLabel] is optional. The default is sorting the head and tail in ascending order");
            System.out.println("sortingLabel:  1: ascending order [default]");
            System.out.println("              -1: descending order");
            System.out.println("              otherwise: no sorting");
            System.exit(1);
        }

        String petFile = args[0];
        String clusterFile = args[1];
        int extensionLength = Integer.parseInt(args[2]);

        String chromLengthInfo = args[3];
        Hashtable<String, Integer> chromLength = new Hashtable<String, Integer>();
        BufferedReader fileChromLength = new BufferedReader(new InputStreamReader(new FileInputStream(new File(chromLengthInfo))));
        String lineTemp;
        while ((lineTemp = fileChromLength.readLine()) != null) {
            String[] temp = lineTemp.split("[ \t][ \t]*");
            if (temp.length > 1) {
                chromLength.put(temp[0], Integer.valueOf(temp[1]));
            }
        }
        int sortingLabel = 1; // ascending order [default]
        if (args.length >= 5) {
            sortingLabel = Integer.parseInt(args[4]);
        }
        Calendar rightNow = Calendar.getInstance();
        PrintWriter clusterFileOut = new PrintWriter(new BufferedWriter(new FileWriter(clusterFile)));
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream(new File(petFile))));
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
            if (fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            PET newPet = new PET(hit1, hit2, weight);

            if (fields.length >= 8) {
                newPet.setIndex(Integer.parseInt(fields[7]));
            }
            CLUSTER cluster = new CLUSTER(newPet, extensionLength, sortingLabel, chromLength);
            clusterFileOut.println(cluster.toString());

            nPETs++;
            if (nPETs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000) + "M PETs read from " + petFile);
                clusterFileOut.flush();
            }
        }

        fileChromLength.close();
        fileIn.close();
        clusterFileOut.close();

    }
}
